function wf = getWaveForms(gwfparams, refilter)
% function wf = getWaveForms(gwfparams)
%
% Extracts individual spike waveforms from the raw datafile, for multiple
% clusters. Returns the waveforms and their means within clusters.
%
% Contributed by C. Schoonover and A. Fink
% Patched by M Macedo-Lima 10/08/20
%
% % EXAMPLE INPUT
% gwfparams.dataDir = '/path/to/data/';    % KiloSort/Phy output folder
% gwfparams.fileName = 'data.dat';         % .dat file containing the raw 
% gwfparams.dataType = 'int16';            % Data type of .dat file (this should be BP filtered)
% gwfparams.nCh = 32;                      % Number of channels that were streamed to disk in .dat file
% gwfparams.wfWin = [-40 41];              % Number of samples before and after spiketime to include in waveform
% gwfparams.nWf = 2000;                    % Number of waveforms per unit to pull out
% gwfparams.spikeTimes =    [2,3,5,7,8,9]; % Vector of cluster spike times (in samples) same length as .spikeClusters
% gwfparams.spikeClusters = [1,2,1,1,1,2]; % Vector of cluster IDs (Phy nomenclature)   same length as .spikeTimes
%
% % OUTPUT
% wf.unitIDs                               % [nClu,1]            List of cluster IDs; defines order used in all wf.* variables
% wf.spikeTimeKeeps                        % [nClu,nWf]          Which spike times were used for the waveforms
% wf.waveForms                             % [nClu,nWf,nCh,nSWf] Individual waveforms
% wf.waveFormsMean                         % [nClu,nCh,nSWf]     Average of all waveforms (per channel)
%                                          % nClu: number of different clusters in .spikeClusters
%                                          % nSWf: number of samples per waveform
%
% % USAGE
% wf = getWaveForms(gwfparams);

% If refilter flag is on, filter again with a bandpass
if refilter
    ops = gwfparams.ops;
    ops.fshigh = 300;
    ops.fslow = 6000;
    ops.comb = 0;
    ops.CAR = 0;
    
    ops.fbinary = ops.fclean;
    
    if ~isfield(ops, 'igood')
        ops.igood=true(size(1:ops.NchanTOT));
    end
    
    if isfield(ops, 'badchannels')
        ops.igood(ops.badchannels) = 0;
    end
    % Create new filtered file
    temp_dir = dir(ops.fclean);
    ops.fclean = fullfile(temp_dir.folder, [temp_dir.name(1:end-4) '300hz.dat']);
    
    caraslab_gpufilter(temp_dir.folder, ops)
    
    fileName = ops.fclean;

else
    ops = gwfparams.ops;
    temp_dir = dir(ops.fclean);
    if isfile(fullfile(temp_dir.folder, [temp_dir.name(1:end-4) '300hz.dat']))  % if no filtering is selected, check if filtered file already exists and use it
        ops.fclean = fullfile(temp_dir.folder, [temp_dir.name(1:end-4) '300hz.dat']);
        fileName = ops.fclean;
    else
        % Load .dat and KiloSort/Phy output
        fileName = fullfile(gwfparams.rawDir,gwfparams.fileName);  
    end
end

filenamestruct = dir(fileName);
dataTypeNBytes = numel(typecast(cast(0, gwfparams.dataType), 'uint8')); % determine number of bytes per sample
nSamp = filenamestruct.bytes/(gwfparams.nCh*dataTypeNBytes);  % Number of samples per channel
wfNSamples = length(gwfparams.wfWin(1):gwfparams.wfWin(end));
mmf = memmapfile(fileName, 'Format', {gwfparams.dataType, [gwfparams.nCh nSamp], 'x'});


chMap = readNPY(fullfile(gwfparams.dataDir, 'channel_map.npy'))+1;               % Order in which data was streamed to disk; must be 1-indexed for Matlab
nChInMap = numel(chMap);

%% phy2
unitIDs = gwfparams.good_clusters;
numUnits = size(gwfparams.good_clusters,1);
spikeTimeKeeps = nan(numUnits,gwfparams.nWf);
allSpikeTimes = cell(numUnits, 1);

%%
waveForms = nan(numUnits,gwfparams.nWf,nChInMap,wfNSamples);
waveFormsMean = nan(numUnits,nChInMap,wfNSamples);

% Added by MML to speed up measurements. Will measure only good clusters
%% phy1
% good_clusters = gwfparams.cluster_quality(strcmp(gwfparams.cluster_quality.group(:,1), 'g'),1);
% good_clusters = table2array(good_clusters);
% [ , good_cluster_idx] = intersect(unitIDs, good_clusters);
%% phy2
good_clusters = gwfparams.good_clusters;
%%

counter = 1;
for curUnitInd=1:numUnits

% for curUnitInd=1:length(good_cluster_idx)
    %% phy1
%     curUnitID = unitIDs(curUnitInd);
    
    %% phy2
    curUnitID = good_clusters(curUnitInd);
    %%
    curSpikeTimes = gwfparams.spikeTimes(gwfparams.spikeClusters==curUnitID);
    allSpikeTimes{curUnitInd} = curSpikeTimes;
    curUnitnSpikes = size(curSpikeTimes,1);
    if ismember(curUnitID, good_clusters)
        spikeTimesRP = curSpikeTimes(randperm(curUnitnSpikes));
        spikeTimeKeeps(curUnitInd,1:min([gwfparams.nWf curUnitnSpikes])) = sort(spikeTimesRP(1:min([gwfparams.nWf curUnitnSpikes])));
        for curSpikeTime = 1:min([gwfparams.nWf curUnitnSpikes])
            try
            tmpWf = mmf.Data.x(1:gwfparams.nCh, ...
                round(spikeTimeKeeps(curUnitInd,curSpikeTime)+ gwfparams.wfWin(1):...
                spikeTimeKeeps(curUnitInd,curSpikeTime)+ gwfparams.wfWin(end)));
            catch ME
                if strcmp(ME.identifier, 'MATLAB:badsubscript')
                    % This catches the error of trying to fetch the last
                    % spike in the recording that spans larger than the
                    % recording itself
                    continue
                else
                    throw(ME)
                end
            end
            waveForms(curUnitInd,curSpikeTime,:,:) = tmpWf(chMap,:);
        end
        waveFormsMean(curUnitInd,:,:) = squeeze(nanmean(waveForms(curUnitInd,:,:,:),2));
    
        % Edited by MML
        disp(['Completed ' int2str(counter) ' units of ' int2str(length(good_clusters)) '.']);
        counter = counter + 1;
    else
        continue
    end
end

% Package in wf struct
wf.unitIDs = unitIDs;
wf.spikeTimeKeeps = spikeTimeKeeps;
wf.waveForms = waveForms;
wf.waveFormsMean = waveFormsMean;

%% MML: phy2
wf.allSpikeTimePoints = allSpikeTimes;
end