function [wf, gwfparams] = get_waveforms_from_folder(cur_savedir, bp_filter, only_good)
    %Load in configuration file (contains ops struct)
    % Catch error if -mat file is not found and skip folder
    try
        temp = load(fullfile(cur_savedir, 'config.mat'));
        ops = temp.ops;
    catch ME
        if strcmp(ME.identifier, 'MATLAB:load:couldNotReadFile')
            fprintf('\nFile not found\n')
            return
        else
            fprintf(ME.identifier)  % file not found has no identifier?? C'mon MatLab...
            fprintf(ME.message)
            return  % Continue here instead of break because I don't know how to catch 'file not found' exception; maybe using ME.message?
        end
    end

    % Define I/O and waveform parameters
    gwfparams.dataDir = cur_savedir;    % KiloSort/Phy output folder 
    gwfparams.ops = ops;

    % Define output name for cluster; this is specific for my
    % file-naming convention and should be tweaked
    split_dir = split(cur_savedir, '/'); 
    subj_id = split(split_dir{end-1}, '-');
    subj_id = join(subj_id(1:3), '-'); 
    subj_id = subj_id{1}; 
    recording_id = split_dir{end};
    prename = [subj_id '_' recording_id];  % this is what goes into the .txt file name

    gwfparams.rawDir = cur_savedir;
    gwfparams.sr = ops.fs;
    gwfparams.nCh = ops.NchanTOT; % Number of channels that were streamed to disk in .dat file
    gwfparams.dataType = 'int16'; % Data type of .dat file (this should be BP filtered)

    gwfparams.wfWin = [round(-(0.001*gwfparams.sr)) round(0.003*gwfparams.sr)]; % Number of samples before and after spiketime to include in waveform
    gwfparams.nWf = 2000; % Max number of waveforms per unit to pull out for averaging
    gwfparams.spikeTimes = readNPY(fullfile(gwfparams.dataDir, 'spike_times.npy')); % Vector of cluster spike times (in samples) same length as .spikeClusters
    gwfparams.spikeClusters = readNPY(fullfile(gwfparams.dataDir, 'spike_clusters.npy')); % Vector of cluster IDs (Phy nomenclature)   same length as .spikeTimes
    gwfparams.channelShanks = readNPY(fullfile(gwfparams.dataDir, 'channel_shanks.npy')); % Vector of cluster shanks
    gwfparams.channelPositions = readNPY(fullfile(gwfparams.dataDir, 'channel_positions.npy')); % Vector of cluster shanks

    gwfparams.chanMap = readNPY(fullfile(gwfparams.dataDir, 'channel_map.npy')); % this is important in esp if you got rid of files. 
    try
        gwfparams.cluster_quality = tdfread(fullfile(gwfparams.dataDir, 'cluster_info.tsv'));
    catch ME
        if strcmp(ME.identifier, 'MATLAB:load:couldNotReadFile')
            fprintf('\nFile not found\n')
            return
        else
            fprintf(ME.identifier)  % file not found has no identifier?? C'mon MatLab...
            fprintf([ME.message '\n'])
            return 
        end
    end

    try
        gwfparams.fileName = dir(ops.fclean).name; % .dat file containing the raw used for sorting
    catch ME
        if strcmp(ME.identifier, 'MATLAB:needMoreRhsOutputs') % this error might happen when file has been moved
            % try to find the fclean file name within the cur_savedir
            split_fclean_path = split(ops.fclean, '/');
            fclean = split_fclean_path{end};
            gwfparams.fileName = fullfile(cur_savedir, fclean); % .dat file containing the raw used for sorting
            % Update ops
            ops.fclean = fullfile(cur_savedir, fclean);
        end
    end
    % Save new config
    save(fullfile(cur_savedir, 'config.mat'), 'ops');

    % Store (redundant) copy in gwfparams
    gwfparams.ops = ops;

    % Get good only for measuring
    if only_good
        gwfparams.good_clusters = gwfparams.cluster_quality.cluster_id(gwfparams.cluster_quality.group(:,1)=='g');
    else
        gwfparams.good_clusters = gwfparams.cluster_quality.cluster_id(gwfparams.cluster_quality.group(:,1)=='g' | gwfparams.cluster_quality.group(:,1)=='m');
    end


    % Get waveforms from .dat
    wf = getWaveForms(gwfparams, bp_filter);  
end