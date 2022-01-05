function caraslab_createconfig(Savedir,chanMap, badchannels, fetch_tstart_from_behav, recording_type, single_dir)
%
% This function sets configuration parameters for kilosort.
% 
% Input variables:
%   Savedir:
%       path to directory where -mat and -csv files will be saved
%
%   chanMap:
%       path to appropriate probe type's channel map
%
%   badchannels:
%       known bad/disconnected channels. These will still be sorted but will not
%       be used for CAR filter. For some mysterious reason to me, Kilosort
%       performs better when bad channels are included in the sorting
%
%   fetch_tstart_from_behav:
%       if 1, program will attempt to find a behavioral file within
%           Savedir/CSV files to grab the first relevant timestamp. This
%           helps eliminate noise from when the animal is plugged in
%           outside of the booth.
%           If Aversive is in the name of the folder, the first spout onset is
%           selected; if Aversive is not in the name of the folder, the first
%           stimulus onset is selected.
%       if 0, the entire recording ([0 Inf]) will be signalded to kilosort
%
%   recording_type: Placeholder used to adjust sampling rate. Not
%       super-useful
%       'synapse' or 'intan'

%Written by ML Caras Mar 26 2019
% Patched by M Macedo-Lima 9/8/20


%% %%%%% Don't edit this; keep scrolling %%%%% %%
%Check that channel map exists
if ~exist(chanMap,'file')
    fprintf('\nCannot find channel map!\n')
    return
end

if nargin > 5
    datafolders = {};
    [~, datafolders{end+1}, ~] = fileparts(single_dir);
else
    %Prompt user to select folders
    datafolders_names = uigetfile_n_dir(Savedir,'Select data directory');
    datafolders = {};
    for i=1:length(datafolders_names)
        [~, datafolders{end+1}, ~] = fileparts(datafolders_names{i});
    end
end
%Load in the channel map and identify bad channels
chandata = load(chanMap);
% badchannels = chandata.chanMap(chandata.connected == 0);

%% Loop through files
for i = 1:numel(datafolders)
    clear temp ops 
    
    cur_path.name = datafolders{i};
    cur_savedir = [Savedir filesep cur_path.name];
    matfilename = fullfile(cur_savedir, [cur_path.name '.mat']);
    infofilename = fullfile(cur_savedir, [cur_path.name '.info']);
    
    % Catch error if -mat file is not found
    try
        temp = load(infofilename, '-mat');
        %Get the sampling rate and number of channels
        ops.fs = temp.epData.streams.RSn1.fs;
        ops.NchanTOT = numel(temp.epData.streams.RSn1.channels);    %both active and dead
    catch ME
        if strcmp(ME.identifier, 'MATLAB:load:couldNotReadFile') || strcmp(ME.identifier, 'MATLAB:nonExistentField')
%             fprintf('\n-mat file not found')  % For concatenated recordings, this will fail. So just hard-code for now
            switch recording_type
                case 'synapse'
                    ops.fs = 24414.0625;
                case 'intan'
                    ops.fs = 30000;
            end
            ops.NchanTOT = length(chandata.chanMap);
        else
            fprintf(ME.identifier)
            fprintf(ME.message)
            continue
        end
    end
    
    %Save the path to the -mat data file (contains raw voltage data)
    ops.rawdata = matfilename;
    
    
    
    %% %%%%% Edit these if necessary %%%%% %%
    % nt0 is the number of points for template matching
    % The default on Kilosort is 61, equivalent to ~2 ms for a 30 kHz sampling
    % If you start seeing double spikes in phy, consider reducing this
    ops.nt0 = round(ops.fs * 0.002);
    
    if ~fetch_tstart_from_behav
        ops.trange = [0 Inf]; % time range to sort
    else
        tstart = fetch_tstart_from_behavior(fullfile(cur_savedir, 'CSV files'));
        ops.trange = [tstart Inf];
    end

    ops.CAR = 0;  % CAR after highpass
    
    ops.kilosort_filter = 0; % filter during kilosort, turned off since we're prefiltering
    
    ops.deline = 0;  % Fails too often; use comb for now
    
    ops.comb = 1;  % Comb filter before highpass
    
    ops.rm_artifacts = 1;  % Remove super high amplitude events
    ops.std_threshold = 25;  % Threshold for artifact rejection (50)
    
    ops.Nchan = ops.NchanTOT - numel(badchannels);              %number of active channels

    ops.fbinary = fullfile(cur_savedir, [cur_path.name '.dat']);
    ops.fclean = fullfile(cur_savedir, [cur_path.name '_CLEAN.dat']);
    
    %Define the channel map and associated parameters
    ops.chanMap  = chanMap;
    ops.badchannels = badchannels;

    % frequency for high pass filtering (150)
    ops.fshigh = 150;   

    % minimum firing rate on a "good" channel (0 to skip)
    ops.minfr_goodchannels = 0;

    % threshold on projections (like in Kilosort1, can be different for last pass like [10 4])
    ops.Th = [10 4];  

    % how important is the amplitude penalty (like in Kilosort1, 0 means not used, 10 is average, 50 is a lot) 
    ops.lam = 10;  

    % splitting a cluster at the end requires at least this much isolation for each sub-cluster (max = 1)
    ops.AUCsplit = 0.9; 

    % minimum spike rate (Hz), if a cluster falls below this for too long it gets removed (1/50)
    ops.minFR = 0; 
    
    % number of samples to average over (annealed from first to second value) 
    ops.momentum = [20 400];

    % spatial constant in um for computing residual variance of spike
    ops.sigmaMask = 30;

    % threshold crossings for pre-clustering (in PCA projection space) (8)
    ops.ThPre = 4;  % This at 4 helps recordings not crash for finding too few spikes
        
    %% danger, changing these settings can lead to fatal errors; Edit  if you know what you're doing
    % options for determining PCs
    ops.spkTh           = -2.5;      % spike threshold in standard deviations (-6)
    ops.reorder         = 1;       % whether to reorder batches for drift correction. 
    ops.nskip           = 25;  % how many batches to skip for determining spike PCs

    ops.GPU                 = 1; % has to be 1, no CPU version yet, sorry
    % ops.Nfilt               = 1024; % max number of clusters
    ops.nfilt_factor        = 8; % max number of clusters per good channel (even temporary ones) (4)
    ops.ntbuff              = 64;    % samples of symmetrical buffer for whitening and spike detection
    
    % MML comment: I'm guessing that the first term 64*1024 has to be a 
    % multiple of 32. Not the whole of ops.NT
    ops.NT                  = 32*1024 + ops.ntbuff; % must be multiple of 32 + ntbuff. This is the batch size (try decreasing if out of memory). 
    
    ops.whiteningRange      = 32; % number of channels to use for whitening each channel
    ops.nSkipCov            = 25; % compute whitening matrix from every N-th batch
    ops.scaleproc           = 200;   % int16 scaling of whitened data
    ops.nPCs                = 3; % how many PCs to project the spikes into
    ops.useRAM              = 0; % not yet available

    %Save configuration file
    configfilename  = fullfile(cur_savedir,'config.mat');
    save(configfilename,'ops')
    fprintf('Saved configuration file: %s\n', configfilename)
end



