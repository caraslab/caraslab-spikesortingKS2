function remove_double_counted_spikes(Savedir, restore_original_npys_first)
    %     From the Allen Brain Institude https://github.com/AllenInstitute/ecephys_spike_sorting
    %     Translated from Python and adapted by M Macedo-Lima, December, 2020

    %     """ Remove putative double-counted spikes from Kilosort outputs
    % 
    %     Inputs:
    %     ------
    %     spike_times : numpy.ndarray (num_spikes x 0)
    %         Spike times in samples 
    %     spike_clusters : numpy.ndarray (num_spikes x 0)
    %         Cluster IDs for each spike time
    %     spike_templates : numpy.ndarray (num_spikes x 0)
    %         Template IDs for each spike time
    %     amplitudes : numpy.ndarray (num_spikes x 0)
    %         Amplitude value for each spike time
    %     channel_map : numpy.ndarray (num_units x 0)
    %         Original data channel for pc_feature_ind array
    %     templates : numpy.ndarray (num_units x num_channels x num_samples)
    %         Spike templates for each unit
    %     pc_features : numpy.ndarray (num_spikes x num_pcs x num_channels)
    %         Pre-computed PCs for blocks of channels around each spike
    %     pc_feature_ind : numpy.ndarray (num_units x num_channels)
    %         Channel indices of PCs for each unit
    %     sample_rate : Float
    %         Sample rate of spike times
    %     params : dict of parameters
    %         'within_unit_overlap_window' : time window for removing overlapping spikes
    %         'between_unit_overlap_window' : time window for removing overlapping spikes
    %         'between_unit_channel_distance' : number of channels over which to search for overlapping spikes
    %     epochs : list of Epoch objects
    %         contains information on Epoch start and stop times
    % 
    %     
    %     Outputs:
    %     --------
    %     spike_times : numpy.ndarray (num_spikes x 0)
    %         Spike times in seconds (same timebase as epochs)
    %     spike_clusters : numpy.ndarray (num_spikes x 0)
    %         Cluster IDs for each spike time
    %     spike_templates : numpy.ndarray (num_spikes x 0)
    %         Template IDs for each spike time
    %     amplitudes : numpy.ndarray (num_spikes x 0)
    %         Amplitude value for each spike time
    %     pc_features : numpy.ndarray (num_spikes x num_pcs x num_channels)
    %         Pre-computed PCs for blocks of channels around each spike
    %     overlap_matrix : numpy.ndarray (num_clusters x num_clusters)
    %         Matrix indicating number of spikes removed for each pair of clusters
    % 
    %     """

    %Prompt user to select folder
    datafolders_names = uigetfile_n_dir(Savedir,'Select data directory');
    datafolders = {};
    for i=1:length(datafolders_names)
        [~, datafolders{end+1}, ~] = fileparts(datafolders_names{i});
    end


    %For each data folder...
    for i = 1:numel(datafolders)
        clear ops rez
        close all
        
        t0 = tic;
        
        cur_path.name = datafolders{i};
        cur_savedir = [Savedir filesep cur_path.name];
        
        % Load in configuration file (contains ops struct) just to get the
        % sampling rate
        % Catch error if -mat file is not found and skips folder
        try
            config_file = load(fullfile(cur_savedir, 'config.mat'));
            ops = config_file.ops;
            
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
        
        fprintf('Removing double-counted spikes from: %s\n', cur_savedir)
        
        if restore_original_npys_first && isfield(rez, 'stored_ks_originals_flag')
            restore_original_npys(cur_savedir, rez);
        end
        
        
        % Read in npy files
        spike_times = readNPY(fullfile(cur_savedir, 'spike_times.npy'));
        spike_clusters = readNPY(fullfile(cur_savedir, 'spike_clusters.npy'));
        spike_templates = readNPY(fullfile(cur_savedir, 'spike_templates.npy'));
        amplitudes = readNPY(fullfile(cur_savedir, 'amplitudes.npy'));
        channel_map = readNPY(fullfile(cur_savedir, 'channel_map.npy'));
        templates = readNPY(fullfile(cur_savedir, 'templates.npy')); 
        pc_features = readNPY(fullfile(cur_savedir, 'pc_features.npy'));
        channel_shanks = readNPY(fullfile(cur_savedir, 'channel_shanks.npy'));

                
        % Save originals in rez.mat if not done before
        if ~isfield(rez, 'stored_ks_originals_flag')
            rez.spike_times = spike_times;
            rez.spike_clusters = spike_clusters;
            rez.spike_templates = spike_templates;
            rez.amplitudes = amplitudes;
            rez.pc_features = pc_features;
            rez.stored_ks_originals_flag = 1;
        end
        
        sample_rate = ops.fs;
        
        unit_list = unique(spike_clusters);
        
        % Use the templates file to grab best channels
        [~, temp_max_idx] = max(abs(max(templates, [], 2) - min(templates,[], 2)), [], 3);
        peak_channels = channel_map(temp_max_idx);

        overlap_matrix = zeros(size(peak_channels));
        
        % Delete one of the spikes that fall within 0.15 ms of each other (from AllenInstitute);
        % Might be extreme but does a pretty good job deleting double-counted spikes
        % MML note: 0.15 ms wasn't doing much for me so I increased to 0.3
        within_unit_overlap_window = 0.0003;
        between_unit_overlap_window = 0.0003;
        
        within_unit_overlap_samples = round(within_unit_overlap_window * sample_rate);
        between_unit_overlap_samples = round(between_unit_overlap_window * sample_rate);

        fprintf('Removing within-unit overlapping spikes...\n')

        spikes_to_remove = [];  % Store indices of overlapping spikes
        for idx1=1:length(unit_list)
            unit_id1 = unit_list(idx1);

            for_unit1 = find(spike_clusters == unit_id1);

            to_remove = find_within_unit_overlap(spike_times(for_unit1), within_unit_overlap_samples);

            overlap_matrix(idx1, idx1) = length(to_remove);

            spikes_to_remove = [spikes_to_remove; for_unit1(to_remove)];
        
        end
        
        if ~isempty(spikes_to_remove)
            [spike_times, spike_clusters, spike_templates, amplitudes, pc_features] = ...
                remove_spikes(spike_times, spike_clusters, spike_templates, ...
                amplitudes, pc_features, spikes_to_remove);
        end
        
        fprintf('Removing between-unit overlapping spikes...\n')

        spikes_to_remove = [];  % Store indices of overlapping spikes
        for idx1=1:(length(unit_list))
            unit_id1 = unit_list(idx1);
            
            % For debugging
%             if unit_list(idx1) ~= 496
%                 continue
%             end
            
            for_unit1 = find(spike_clusters == unit_id1);
            if isempty(for_unit1)
                continue
            end
            template_unit1 = spike_templates(for_unit1);
            template_unit1 = template_unit1(1);
            
            for idx2=(idx1+1):length(unit_list)
                unit_id2 = unit_list(idx2);
                
                % For debugging
%                 if unit_list(idx2) ~= 512
%                     continue
%                 end

                % MML edit: 
                % Avoid repetition and only look at them if they're on the
                % same shank
                for_unit2 = find(spike_clusters == unit_id2);
                if isempty(for_unit2)
                    continue
                end
                template_unit2 = spike_templates(for_unit2);
                template_unit2 = template_unit2(1);
                if channel_shanks(channel_map == peak_channels(template_unit1+1)) == channel_shanks(channel_map == peak_channels(template_unit2+1))
                    
                    [to_remove1, to_remove2] = find_between_unit_overlap(spike_times(for_unit1), spike_times(for_unit2), between_unit_overlap_samples);

                    overlap_matrix(idx1, idx2) = length(to_remove1) + length(to_remove2);

                    spikes_to_remove = [spikes_to_remove; [for_unit1(to_remove1); for_unit2(to_remove2)]];
                end
            end
        end
        
        if ~isempty(spikes_to_remove)
            [spike_times, spike_clusters, spike_templates, amplitudes, pc_features] = ...
                remove_spikes(spike_times, spike_clusters, spike_templates, ...
                amplitudes, pc_features, spikes_to_remove);
        end
        
        writeNPY(spike_times, fullfile(cur_savedir, 'spike_times.npy'));
        writeNPY(spike_clusters, fullfile(cur_savedir, 'spike_clusters.npy'));
        writeNPY(spike_templates, fullfile(cur_savedir, 'spike_templates.npy'));
        writeNPY(amplitudes, fullfile(cur_savedir, 'amplitudes.npy'));
        writeNPY(pc_features, fullfile(cur_savedir, 'pc_features.npy'));
        writeNPY(overlap_matrix, fullfile(cur_savedir, 'overlap_matrix.npy'));
        
        save(fullfile(cur_savedir, 'rez.mat'), 'rez');
        
        tEnd = toc(t0);
        fprintf('Done in: %d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));
    
    end
    
        
    function [to_remove] = find_within_unit_overlap(spike_train, overlap_window)
%         """
%         Finds overlapping spikes within a single spike train.
%     
%         Parameters
%         ----------
%         spike_train : numpy.ndarray
%             Spike times (in samples)
%         overlap_window : int
%             Number of samples to search for overlapping spikes
%     
%     
%         Outputs
%         -------
%         spikes_to_remove : numpy.ndarray
%             Indices of overlapping spikes in spike_train
%     
%         """
        to_remove = find(diff(spike_train) < overlap_window);
    end
    
    
    function [to_remove1, to_remove2] = find_between_unit_overlap(spike_train1, spike_train2, overlap_window)
% 
%         """
%         Finds overlapping spikes between two spike trains
% 
%         Parameters
%         ----------
%         spike_train1 : numpy.ndarray
%             Spike times (in samples)
%         spike_train2 : numpy.ndarray
%             Spike times (in samples)
%         overlap_window : int
%             Number of samples to search for overlapping spikes
% 
% 
%         Outputs
%         -------
%         spikes_to_remove1 : numpy.ndarray
%             Indices of overlapping spikes in spike_train1
%         spikes_to_remove2 : numpy.ndarray
%             Indices of overlapping spikes in spike_train2
% 
%         """

        spike_train = [spike_train1; spike_train2];
        original_inds = [1:length(spike_train1) 1:length(spike_train2)];
        cluster_ids = [ zeros(length(spike_train1), 1); ones(length(spike_train2),1) ];

        [~, spike_train_order] = sort(spike_train);
        sorted_train = spike_train(spike_train_order);
        sorted_inds = original_inds(spike_train_order);
        sorted_inds = sorted_inds(2:end);

        sorted_cluster_ids = cluster_ids(spike_train_order);
        sorted_cluster_ids  = sorted_cluster_ids(2:end);

        to_remove = diff(sorted_train) < overlap_window;

        to_remove1 = sorted_inds(to_remove & (sorted_cluster_ids == 0));
        to_remove2 = sorted_inds(to_remove & (sorted_cluster_ids == 1));
    end

end