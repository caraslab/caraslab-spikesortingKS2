function plot_unit_shanks(Savedir, show_plots, bp_filter, load_previous_gwfparams, plot_mean, plot_std, plot_wf_samples)
    % This function reads the probe geometry in channel map and outputs the
    % spike means and SEM organized in space. If filter_300hz==0, it will
    % search for the 300hz bandpass filtered file. Otherwise, it will filter
    % again
    % Adapted parts from Cortexlab repository to read waveforms
    % By M Macedo-Lima, November 2020
    
    % KNOWN BUG:
    %   It's not perfect. If waveforms are too large, the scaling factor
    %   clobbers them together. Low priority fix right now...

    % Select folders to run
    %Prompt user to select folder
    datafolders_names = uigetfile_n_dir(Savedir,'Select data directory');
    datafolders = {};
    for dummy_idx=1:length(datafolders_names)
        [~, datafolders{end+1}, ~] = fileparts(datafolders_names{dummy_idx});
    end


    %For each data folder...
    for folder_idx = 1:numel(datafolders)
        clear ops rez
        close all
        
        cur_path.name = datafolders{folder_idx};
        cur_savedir = [Savedir filesep cur_path.name];

        %Load in configuration file (contains ops struct)
        % Catch error if -mat file is not found and skip folder
        try
            load(fullfile(cur_savedir, 'config.mat'));
            readNPY(fullfile(cur_savedir, 'spike_times.npy'));
        catch ME
            if strcmp(ME.identifier, 'MATLAB:load:couldNotReadFile')
                fprintf('\nFile not found\n')
                continue
            else
                fprintf(ME.identifier)  % file not found has no identifier?? C'mon MatLab...
                fprintf(ME.message)
                continue  % Continue here instead of break because I don't know how to catch 'file not found' exception; maybe using ME.message?
            end
        end

        %Start timer
        t0 = tic;
        
        %% Define I/O and waveform parameters
        
        % Define output name for cluster; this is specific for my
        % file-naming convention and should be tweaked
        split_dir = split(cur_savedir, '/'); 
        subj_id = split(split_dir{end-1}, '-');
        subj_id = join(subj_id(1:3), '-'); 
        subj_id = subj_id{1}; 
        recording_id = split_dir{end};
        prename = [subj_id '_' recording_id];  % this is what goes into the .txt file name

        if load_previous_gwfparams
            [wf, gwfparams] = try_load_previous_gwfparams(cur_savedir, bp_filter, 1);
        else
            fprintf('Running wf extraction...\n')
            [wf, gwfparams] = get_waveforms_from_folder(cur_savedir, bp_filter, 1);
            % Save gwfparams and wf for future use
            fprintf('Saving gwfparams and wf structs to mat file\n')
            save(fullfile(day_dir, 'extracted_wfs.mat'), 'gwfparams', 'wf', '-v7.3');
        end
        
        % In an older version of get_waveforms_from_folder I forgot to add
        % this field to gwfparams. Will become obsolete
        gwfparams.channelPositions = readNPY(fullfile(gwfparams.dataDir, 'channel_positions.npy')); % Vector of cluster shanks

        good_cluster_idx = wf.unitIDs; % store their Phy IDs

        % Delete old plots in folder
        old_files = dir([cur_savedir '/*shankWaveforms*pdf']);
        for dummy_idx=1:numel(old_files)
            delete(fullfile(old_files(dummy_idx).folder, old_files(dummy_idx).name));
        end
        
        %% Plot
        for wf_idx=1:length(good_cluster_idx)
            if show_plots
                figure('Name', prename);
            else
                figure('Name', prename, 'visible', 'off');
            end
            
            cluster_phy_id = good_cluster_idx(wf_idx);

            cur_shank = gwfparams.cluster_quality.sh(gwfparams.cluster_quality.cluster_id == cluster_phy_id);
            shank_channels_0ind = gwfparams.chanMap(gwfparams.channelShanks == cur_shank);
            [~, shank_channels_idx] = intersect(gwfparams.chanMap, shank_channels_0ind);
            
            % Rescale x and y to 0-1
            scaled_channel_positions = gwfparams.channelPositions;
            scaled_channel_positions(:,1) = rescale(scaled_channel_positions(:,1));
            scaled_channel_positions(:,2) = rescale(scaled_channel_positions(:,2));
         
            % Normalize all channels by the best channel
            best_channel_0ind = gwfparams.cluster_quality.ch(gwfparams.cluster_quality.cluster_id == cluster_phy_id);
            best_chanMap_idx = find(gwfparams.chanMap == best_channel_0ind);
            best_wf_mean = squeeze(wf.waveFormsMean(wf_idx, best_chanMap_idx,:));
            abs_best_wf_mean = abs(best_wf_mean);
            cur_wf_mean_peak_idx = abs_best_wf_mean == max(abs_best_wf_mean);
            norm_peak_value = abs(best_wf_mean(cur_wf_mean_peak_idx)) / 100;

            % If more than 14 channels plot 13 closest channels to best channel
            % Use an euclidean distances matrix to figure this out
            if length(shank_channels_idx) > 14
                euclidean_distances = zeros(length(scaled_channel_positions), length(scaled_channel_positions));

                for ch_i=1:length(shank_channels_idx)
                    for ch_j=1:length(shank_channels_idx)
                        euclidean_distances(shank_channels_idx(ch_i), shank_channels_idx(ch_j)) = ...
                            sqrt( (scaled_channel_positions(shank_channels_idx(ch_i),1)-scaled_channel_positions(shank_channels_idx(ch_j),1))^2 +...
                            (scaled_channel_positions(shank_channels_idx(ch_i),2)-scaled_channel_positions(shank_channels_idx(ch_j),2))^2 );
                    end
                end

                closest_distances = nonzeros(sort(euclidean_distances(:,best_chanMap_idx)));
                closest_distances = closest_distances(1:min(14, length(closest_distances)));  % grab 5 non-zero firsts
                [~, closest_channels_idx] = intersect(euclidean_distances(:,best_chanMap_idx), closest_distances);

                closest_channels_idx = [closest_channels_idx; best_chanMap_idx];
                % Exclude channels not part of the current shank
                [shank_channels_idx, ~] = intersect(shank_channels_idx, closest_channels_idx);
            end
            
            
            for shank_channel_idx=1:numel(shank_channels_idx)
                cur_chanMap_idx = shank_channels_idx(shank_channel_idx);
                
                % Grab x and y positions; rescale them to look closer
                % together; might need to be tweaked for each probe type
                x_offset = scaled_channel_positions(cur_chanMap_idx, 1)*120;
                y_offset = scaled_channel_positions(cur_chanMap_idx, 2)*5000;
                
                % Squeeze out and store raw waveforms and averages
                snip_points = 20;  % snip some points at the end of template
                cur_wfs = wf.waveForms(wf_idx, :, cur_chanMap_idx,1:end-snip_points);
                cur_wfs = squeeze(cur_wfs);
                cur_wfs = rmmissing(cur_wfs, 1);
                
                cur_wf_std = std(cur_wfs./norm_peak_value, 0, 1);

                cur_wf_mean = wf.waveFormsMean(wf_idx, cur_chanMap_idx,1:end-snip_points);
                cur_wf_mean = squeeze(cur_wf_mean);

                % 10x upsample mean waveform with spline interpolation
                samplingRateIncrease = 10;
                newXSamplePoints = linspace(1, length(cur_wf_mean), length(cur_wf_mean) * samplingRateIncrease);
                cur_wf_mean_upsample = spline(1:length(cur_wf_mean), cur_wf_mean, newXSamplePoints);
                cur_wf_std_upsample = spline(1:length(cur_wf_std), cur_wf_std, newXSamplePoints);
                
                x_time = linspace(0, length(cur_wf_mean) / gwfparams.sr, length(cur_wf_mean));
                x_time = x_time * 1000; % in ms

                upsampled_x_time = linspace(0, length(cur_wf_mean_upsample) / gwfparams.sr / samplingRateIncrease, length(cur_wf_mean_upsample));
                upsampled_x_time = upsampled_x_time*1000;
                
                % Normalize amplitude so that peak is at -1 or 1
                cur_wf_mean_upsample = cur_wf_mean_upsample ./ norm_peak_value;
                
                title(['Cluster' num2str(cluster_phy_id)]);

                hold on
                
                if plot_wf_samples
                  % If you want to output individual spikes around the mean instead of SEM
                    for spike_idx=1:min(size(cur_wfs, 1), 50)  % plot a max of 50 spikes
                        cur_wf = cur_wfs(spike_idx,:);
                        if isnan(mean(cur_wfs(spike_idx, 1:5)))
                            continue  % sometimes waveforms are NaN which is weird...
                        end
                        % Normalize amplitude of individual spike using the average
                        cur_wf = cur_wf ./ norm_peak_value;

                        cax = plot(x_time + x_offset, cur_wf + y_offset, 'black'); % Plot 0-centered traces
                        cax.Color(4)=0.1;
                    end
                end
                
                if plot_std
                    std_patch_top = (cur_wf_mean_upsample + cur_wf_std_upsample) + y_offset;
                    std_patch_bottom = (cur_wf_mean_upsample - cur_wf_std_upsample) + y_offset ;
                    patch_cax = patch([upsampled_x_time + x_offset fliplr(upsampled_x_time + x_offset)],...
                        [std_patch_bottom fliplr(std_patch_top)], 'black', 'EdgeColor','none');
                    patch_cax.FaceAlpha=0.3;
                end
                
                if plot_mean
                    plot(upsampled_x_time + x_offset, cur_wf_mean_upsample + y_offset, 'black', 'linewidth', 0.5);
                end
            end

            axis off;
           print([gwfparams.dataDir '\' prename '_shankWaveforms' '_cluster' num2str(cluster_phy_id)], '-dpdf', '-bestfit', '-painters');
        end
    end
end
