function id_noise_templates_wrapper(Savedir, sel, restore_original_npys_first, py_code_folder)

% IN DEVELOPMENT %

%     From the Allen Brain Institude https://github.com/AllenInstitute/ecephys_spike_sorting
%     Adapted to run from MatLab by M Macedo-Lima, January, 2021

%     
%     """
%     Uses a set of heuristics to identify noise units based on waveform shape
% 
%     Inputs:
%     -------
%     cluster_ids : all unique cluster ids
%     templates : template for each unit output by Kilosort
%     channel_map : mapping between template channels and actual probe channels
% 
%     Outputs:
%     -------
%     cluster_ids : same as input
%     is_noise : boolean array, True at index of noise templates
% 
%     """

    % Load python version. >3.6 does not seem to work with this version of
    % MatLab
    try
        pe = pyenv('Version', '/usr/bin/python3.7');
    catch ME
        if strcmp(ME.identifier, 'MATLAB:Pyenv:PythonLoadedInProcess')
            % Python already loaded; skip
        else
            throw(ME)
        end
    end

    py.sys.setdlopenflags(int32(10));        % Set RTLD_NOW and RTLD_DEEPBIND
    
    % Add folder to python system path.
    insert(py.sys.path, int64(0), py_code_folder);
    % Load module
    func = py.importlib.import_module('id_noise_templates');
    py.importlib.reload(func);
    
    % Select folders
    if ~sel
        datafolders = caraslab_lsdir(Savedir);
        datafolders = {datafolders.name};

    elseif sel  
        %Prompt user to select folder
        datafolders_names = uigetfile_n_dir(Savedir,'Select data directory');
        datafolders = {};
        for i=1:length(datafolders_names)
            [~, datafolders{end+1}, ~] = fileparts(datafolders_names{i});
        end
    end

    %For each data folder...
    for i = 1:numel(datafolders)
        clear ops rez params
        close all
        
        
        t0 = tic;
        
        cur_path.name = datafolders{i};
        cur_savedir = [Savedir filesep cur_path.name];
        
        fprintf('Removing noise from: %s\n', cur_savedir);

        % Load in rez file
        % Catch error if -mat file is not found and skips folder
        try
            rez_file = load(fullfile(cur_savedir, 'rez.mat'));
            rez = rez_file.rez;
        catch ME
            if strcmp(ME.identifier, 'MATLAB:load:couldNotReadFile')
                fprintf('\n-mat file not found\n')
                continue
            else
                fprintf(ME.identifier)
                fprintf(ME.message)
                break
            end
        end
                
        if restore_original_npys_first && isfield(rez, 'stored_ks_originals_flag')
            restore_original_npys(cur_savedir, rez);
        end
        
        % Read in npy files
        spike_times = readNPY(fullfile(cur_savedir, 'spike_times.npy'));
        spike_clusters = readNPY(fullfile(cur_savedir, 'spike_clusters.npy'));
        try
            cluster_ids = tdfread(fullfile(cur_savedir, 'cluster_info.tsv'));
            cluster_ids = int16(cluster_ids.cluster_id);
        catch ME
            if strcmp(ME.identifier, 'stats:tdfread:OpenFailed')
                cluster_ids = unique(spike_clusters);
            end
        end

        spike_templates = readNPY(fullfile(cur_savedir, 'spike_templates.npy'));
        unwhitening_mat = readNPY(fullfile(cur_savedir, 'whitening_mat_inv.npy'));
        
        amplitudes = readNPY(fullfile(cur_savedir, 'amplitudes.npy'));
        channel_map = readNPY(fullfile(cur_savedir, 'channel_map.npy'));
        templates = readNPY(fullfile(cur_savedir, 'templates.npy')); 
        template_zero_padding = 18;
        templates = templates(:, template_zero_padding:end, :);
        pc_features = readNPY(fullfile(cur_savedir, 'pc_features.npy'));
        channel_shanks = readNPY(fullfile(cur_savedir, 'channel_shanks.npy'));
        
        if ~isfield(rez, 'stored_ks_originals_flag')
            rez.spike_times = spike_times;
            rez.spike_clusters = spike_clusters;
            rez.spike_templates = spike_templates;
            rez.amplitudes = amplitudes;
            rez.pc_features = pc_features;
            rez.stored_ks_originals_flag = 1;
        end
        
        unwhitened_temps = zeros(size(templates));
        for temp_idx= 1:size(templates, 1)
            unwhitened_temps(temp_idx,:,:) = squeeze(templates(temp_idx,:,:)) * unwhitening_mat;
        end
        
        % Load in parameters
        params.smoothed_template_amplitude_threshold = 0.2;  % , help='Fraction of max amplitude for calculating spread')
        params.template_amplitude_threshold = 0.2;  % , help='Fraction of max amplitude for calculating spread')
        params.smoothed_template_filter_width = int16(2);  % , help='Smoothing window for calculating spread')
        
        params.min_spread_threshold = int16(1);  % , help='Minimum number of channels for a waveform to be considered good')
        params.mid_spread_threshold = int16(5);  % , help='Over this channel spread, waveform shape must be considered')
        params.max_spread_threshold = int16(20);  % , help='Maximum channel spread for a good unit')
        
        params.channel_amplitude_thresh = 0.25;  % , help='Fraction of max amplitude for considering channels in spatial peak detection')
        params.peak_height_thresh = 0.2;  % , help='Minimum height for spatial peak detection')
        params.peak_prominence_thresh = 0.2;  % , help='Minimum prominence for spatial peak detection')
        params.peak_channel_range = int16(24);  % , help='Range of channels for detecting spatial peaks')
        params.peak_locs_std_thresh = 5;  % , help='Maximum standard deviation of peak locations for good units')

        params.min_temporal_peak_location = int16(8);  % , help='Minimum peak index for good unit')
        params.max_temporal_peak_location = int16(24);  % , help='Maximum peak index for good unit')

        params.template_shape_channel_range = int16(12);  % , help='Range of channels for checking template shape')
        params.wavelet_index = int16(2);  % , help='Wavelet index for noise template shape detection')
        params.min_wavelet_peak_height = int16(0);  % , help='Minimum wavelet peak height for good units')
        params.min_wavelet_peak_loc = int16(15);  % , help='Minimum wavelet peak location for good units')
        params.max_wavelet_peak_loc = int16(25);  % , help='Maximum wavelet peak location for good units')
        
        params.multiprocessing_worker_count = int16(4);  % not currently in use; makes MatLab crash
        
        % shuttle MatLab variables to Python
        cluster_ids = py.numpy.array(unique(cluster_ids));
        unwhitened_temps = py.numpy.array(unwhitened_temps);
        channel_map = py.numpy.array(int16(channel_map));
        params = py.dict(params);
        
        % Run noise detection        
        is_noise = py.id_noise_templates.id_noise_templates(cluster_ids, unwhitened_temps, channel_map, params);
        
        % Shuttle Python variable back to MatLab
        is_noise = boolean(is_noise);
        
        % Remove noise templates from all npys
        [~, noise_templates] = find(is_noise);
        noise_templates = noise_templates-1;
        [spikes_to_remove, ~] = ismember(spike_templates, noise_templates);
        
        if ~isempty(spikes_to_remove)
            [spike_times, spike_clusters, spike_templates, amplitudes, pc_features] = ...
                remove_spikes(spike_times, spike_clusters, spike_templates, ...
                amplitudes, pc_features, spikes_to_remove);
        end
        
        % Write new npys
        writeNPY(spike_times, fullfile(cur_savedir, 'spike_times.npy'));
        writeNPY(spike_clusters, fullfile(cur_savedir, 'spike_clusters.npy'));
        writeNPY(spike_templates, fullfile(cur_savedir, 'spike_templates.npy'));
        writeNPY(amplitudes, fullfile(cur_savedir, 'amplitudes.npy'));
        writeNPY(pc_features, fullfile(cur_savedir, 'pc_features.npy'));
        
        % Write new rez
        save(fullfile(cur_savedir, 'rez.mat'), 'rez');
        
        tEnd = toc(t0);
        fprintf('Done in: %d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));
    
    end
    
