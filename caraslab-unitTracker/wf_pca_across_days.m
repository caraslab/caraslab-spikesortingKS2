function wf_pca_across_days(Savedir, show_plots, bp_filter, load_previous_gwfparams)
    % This function loops through recording folders and compares waveforms.
    % Option to compare firing rates and autocorrelograms between
    % consecutive days if they occured on the same shank.
    % Created by M Macedo-Lima December, 2020

    %Prompt user to select folder
    datafolders_names = uigetfile_n_dir(Savedir,'Select data directory');
    datafolders = {};
    for i=1:length(datafolders_names)
        [~, datafolders{end+1}, ~] = fileparts(datafolders_names{i});
    end

    % It is important that folder names contain at least a date as a second
    % item after an underscore (_). Adjust this part for your own naming
    % conventions if different
    
    % Sort according to dates and times in folder names
    date_time = regexp(datafolders, '\d+', 'match');
    recording_dates = [];
    recording_times = [];
    try
        for i=1:length(date_time)
           recording_dates = [recording_dates str2num(date_time{i}{1})];
           if length(date_time{i}) > 1
              recording_times = [recording_times str2num(date_time{i}{2})];
           else
              recording_times = [recording_times 0];
           end
        end
    catch ME
        if strcmp(ME.identifier, 'MATLAB:badsubscript')
            fprintf('\nfile not found\n')
            return
        else
            fprintf(ME.identifier)
            fprintf(ME.message)
            return
        end
    end

    % now sort hierarchically first date, then time
    temp_cell = horzcat(datafolders', num2cell([recording_dates' recording_times']) );
    temp_cell = sortrows(temp_cell, [2 3]);
    datafolders = temp_cell(:,1)';

    unique_days = unique(recording_dates);
    if length(unique_days) < 2
        fprintf('Only one day was found... Nothing to compare\n')
        return
    end
    
    % Cycle through folders
    for day1_idx=1:length(unique_days)-1
        close all
        
        day2_idx=day1_idx+1;
        
        day1_dir = fullfile(Savedir, datafolders{day1_idx});
        day2_dir = fullfile(Savedir, datafolders{day2_idx});

        output_dir = [num2str(unique_days(day1_idx)) 'vs' num2str(unique_days(day2_idx)) '_pca_comparison'];
        output_file_name = [num2str(unique_days(day1_idx)) 'vs' num2str(unique_days(day2_idx))];
        mkdir(fullfile(Savedir, output_dir));
        
        % Load waveforms if gwfparams cannot be found in folder or if
        % load_previous_gwfparams flag is 0
        if load_previous_gwfparams
            [day1_wf, day1_gwfparams] = try_load_previous_gwfparams(day1_dir, bp_filter, 1);
            [day2_wf, day2_gwfparams] = try_load_previous_gwfparams(day2_dir, bp_filter, 1);
        else
            fprintf('Running wf extraction...\n')
            [day1_wf, day1_gwfparams] = get_good_waveforms_from_folder(day1_dir, bp_filter, 1);
            [day2_wf, day2_gwfparams] = get_good_waveforms_from_folder(day2_dir, bp_filter, 1);
            
            % Save gwfparams and wf for future use
            fprintf('Saving gwfparams and wf structs to mat file\n')
            wf = day1_wf;
            gwfparams = day1_gwfparams;
            save(fullfile(day_dir, 'extracted_wfs.mat'), 'gwfparams', 'wf', '-v7.3');
            wf = day2_wf;
            gwfparams = day2_gwfparams;
            save(fullfile(day_dir, 'extracted_wfs.mat'), 'gwfparams', 'wf', '-v7.3');
            
            clear wf gwfparams
        end
        
        % Snip waveforms from 0.5 to 1.5 ms to increase PCA discriminability
        % Tweak this if needed
        start_wf_sample = round(0.0005 * day1_gwfparams.ops.fs);
        end_wf_sample = start_wf_sample + round(0.0015 * day1_gwfparams.ops.fs);
       
        %% Loop through waveforms, check if same shank, and compare
        for day1_wf_idx=1:length(day1_wf.unitIDs)
            % Normalize wfs around 0 and so that peak is at -1 or 1
            [day1_cur_wfs, day1_cur_wf_mean, day1_shank] = normalize_wfs(day1_wf, day1_gwfparams, day1_wf_idx);
            
            % Snip around peak-trough for PCA
            day1_cur_wfs = day1_cur_wfs(:,start_wf_sample:end_wf_sample);
            day1_cur_wf_mean = day1_cur_wf_mean(start_wf_sample:end_wf_sample);
            
            for day2_wf_idx=1:length(day2_wf.unitIDs)
                [day2_cur_wfs, day2_cur_wf_mean, day2_shank] = normalize_wfs(day2_wf, day2_gwfparams, day2_wf_idx);
                
                % Skip if units belong to different shanks
                if day1_shank ~= day2_shank
                    continue
                else   
                    day2_cur_wfs = day2_cur_wfs(:,start_wf_sample:end_wf_sample);
                    day2_cur_wf_mean = day2_cur_wf_mean(start_wf_sample:end_wf_sample);
                    
                    %% ignore these metrics for now
%                     %% Waveform similarity
%                     % As in Fraser & Schwartz, 2012 J. Neurosci. apud
%                     % Jackson & Fetz, 2007
%                     % "Peak value of the cross-correlogram between the 
%                     % average waveform shape in session 1 vs. session 2. 
%                     % [...] The resulting coefficient is Fisher transformed
%                     % (the arc tangent of the hypotenuse function) to make
%                     % it more normally distributed
%                     wf_corr = xcorr(day1_cur_wf_mean, day2_cur_wf_mean, 'normalized');
%                     wf_zcorr = atanh(max(wf_corr));  % Fisher transform
%                     
%                     %% Autocorrelogram similarity
%                     % Autocorrelogram parameters in s
%                     bin_size = 0.005;
%                     max_lag = 0.1;
%                     p_corr = ...
%                         autocorrelogram_similarity(day1_wf.allSpikeTimePoints{day1_wf_idx}, ...
%                         day2_wf.allSpikeTimePoints{day2_wf_idx}, ...
%                         day1_gwfparams.ops, day2_gwfparams.ops, ...
%                         max_lag, bin_size);
%                     
%                     %% Firing rate similarity
%                     % "Difference between the log of the mean rates"
%                     fr_diff = firing_rate_similarity(day1_wf.allSpikeTimePoints{day1_wf_idx}, ...
%                                               day2_wf.allSpikeTimePoints{day2_wf_idx}, ...
%                                               day1_gwfparams.ops, day2_gwfparams.ops);
%                     
                    
                    %% Plots
                    if show_plots
                        figure;
                    else
                        figure('visible', 'off');
                    end
                    
                    [coeffs, scores] = pca_plots();
                    
                    % Save figure
                    print(fullfile(Savedir, output_dir, [output_file_name '_' num2str(day1_wf.unitIDs(day1_wf_idx)) 'vs' ...
                        num2str(day2_wf.unitIDs(day2_wf_idx))]), '-dpdf', '-bestfit', '-painters');
                    
                    % Save PCA parameters
                    TT = array2table([ [repmat(day1_wf.unitIDs(day1_wf_idx), size(day1_cur_wfs, 1), 1); ...
                                        repmat(day2_wf.unitIDs(day2_wf_idx), size(day2_cur_wfs, 1), 1)],...
                                      scores(:,1), scores(:,2), scores(:,3)],...
                        'VariableNames', ...
                        {'Cluster' 'PC1' 'PC2' 'PC3'});
                    writetable(TT, fullfile(Savedir, output_dir, [output_file_name '_' num2str(day1_wf.unitIDs(day1_wf_idx)) 'vs' ...
                        num2str(day2_wf.unitIDs(day2_wf_idx)) '_PCA.csv']));
                    
                    close all
                end
            end
        
        end
        %% Clear RAM
        clear day1_wf day2_wf day1_gwfparams day2_gwfparams
    end

    
    function [coeffs, scores] = pca_plots()
        % Concatenate wfs in a new variable for PCA
        all_wfs = [day1_cur_wfs; day2_cur_wfs];
        % Keep track of indices
        original_inds = [1:size(day1_cur_wfs, 1) 1:size(day2_cur_wfs, 1)];
        % Keep track of clusters
        cluster_ids = [ zeros(size(day1_cur_wfs, 1), 1); ones(size(day2_cur_wfs, 1),1) ];

        % X axis in seconds for waveform o=plots
        x_time = linspace(0, length(day1_cur_wfs(1,:)) / day1_gwfparams.sr, length(day1_cur_wfs(1,:)));

        % Plot waveforms
        w1 = subplot(4, 3, 1);
        w2 = subplot(4, 3, 2);
        w1vsw2 = subplot(4, 3, 3);

        hold(w1, 'on')
        hold(w2, 'on')
        hold(w1vsw2, 'on')

        ylims = 2;
        ylim(w1, [-ylims ylims])  
        ylim(w2, [-ylims ylims])  
        ylim(w1vsw2, [-ylims ylims])  

        wf1_colors = {'DarkTurquoise', 'DarkCyan'};
        wf2_colors = {'Amethyst', 'DarkMagenta'};
        for spike_idx=1:min(size(day1_cur_wfs, 1), 500)  % plot a max of 500 spikes
            cax = plot(w1, x_time, day1_cur_wfs(spike_idx,:), 'color', rgb(wf1_colors{1})); % Plot 0-centered traces
            cax.Color(4)=0.1;

            cax = plot(w1vsw2, x_time, day1_cur_wfs(spike_idx,:), 'color', rgb(wf1_colors{1})); % Plot 0-centered traces
            cax.Color(4)=0.1;
        end

        for spike_idx=1:min(size(day2_cur_wfs, 1), 500) % plot a max of 500 spikes
            cax = plot(w2, x_time, day2_cur_wfs(spike_idx,:), 'color', rgb(wf2_colors{1})); % Plot 0-centered traces
            cax.Color(4)=0.1;

            cax = plot(w1vsw2, x_time, day2_cur_wfs(spike_idx,:), 'color', rgb(wf2_colors{1})); % Plot 0-centered traces
            cax.Color(4)=0.1;
        end

        % Plot means
        plot(w1, x_time, day1_cur_wf_mean, 'color', rgb(wf1_colors{2}));
        plot(w2, x_time, day2_cur_wf_mean, 'color', rgb(wf2_colors{2}));
        plot(w1vsw2, x_time, day1_cur_wf_mean, 'color', rgb(wf1_colors{2}));
        plot(w1vsw2, x_time, day2_cur_wf_mean, 'color', rgb(wf2_colors{2}));

        % Turn off axes
        set(w1.XAxis,'visible','off'); set(w1.YAxis,'visible','off')
        set(w2.XAxis,'visible','off'); set(w2.YAxis,'visible','off')
        set(w1vsw2.XAxis,'visible','off'); set(w1vsw2.YAxis,'visible','off')

        % Titles
        title(w1, ['Cluster ' num2str(day1_wf.unitIDs(day1_wf_idx))])
        title(w2, ['Cluster ' num2str(day2_wf.unitIDs(day2_wf_idx))])
        title(w1vsw2, 'Overlay')

        % Run PCA
        [coeffs, scores] = pca(all_wfs, 'Centered', true);
        
        % PC1 vs PC2
        pc_x = 1;
        pc_y = 2;
        pca_ax = subplot(4, 3, [4 7 10]);
        % Remove outliers (e.g. noisy waveforms)
        euc_distances = sqrt(scores(:,pc_x).^2 + scores(:,pc_y).^2);
        is_outlier = isoutlier(euc_distances, 'grubbs');
        
        h = biplot(pca_ax, coeffs(:,[pc_x pc_y]),'Scores', scores(~is_outlier,[pc_x pc_y]));
        xlabel(pca_ax, 'PC1');
        ylabel(pca_ax, 'PC2');
        % Remove loading vectors
        for k = 1:size(all_wfs, 2)*2
            h(k).LineStyle = 'none';
            h(k).MarkerEdgeColor = 'none';
        end
        % Recolor markers
        for k = size(all_wfs, 2)*2+1:size(all_wfs, 2)*2+size(all_wfs(~is_outlier), 1)
            reset_idx = k - size(all_wfs, 2)*2;
            if cluster_ids(reset_idx) == 0
                h(k).MarkerEdgeColor = rgb(wf1_colors{1});
            else
                h(k).MarkerEdgeColor = rgb(wf2_colors{1});
            end
            h(k).MarkerSize = 5;
        end
        
        % PC1 vs PC3
        pca_ax = subplot(4, 3, [5 8 11]);
        pc_x = 1;
        pc_y = 3;
        % Remove outliers (e.g. noisy waveforms)
        euc_distances = sqrt(scores(:,pc_x).^2 + scores(:,pc_y).^2);
        is_outlier = isoutlier(euc_distances, 'grubbs');
        
        h = biplot(pca_ax, coeffs(:,[pc_x pc_y]),'Scores', scores(~is_outlier,[pc_x pc_y]));
        xlabel(pca_ax, 'PC1');
        ylabel(pca_ax, 'PC3');
        % Remove loading vectors
        for k = 1:size(all_wfs, 2)*2
            h(k).LineStyle = 'none';
            h(k).MarkerEdgeColor = 'none';
        end
        % Recolor markers
        for k = size(all_wfs, 2)*2+1:size(all_wfs, 2)*2+size(all_wfs(~is_outlier), 1)
            reset_idx = k - size(all_wfs, 2)*2;
            if cluster_ids(reset_idx) == 0
                h(k).MarkerEdgeColor = rgb(wf1_colors{1});
            else
                h(k).MarkerEdgeColor = rgb(wf2_colors{1});
            end
            h(k).MarkerSize = 5;
        end

        
        % PC2 vs PC3
        pca_ax = subplot(4, 3, [6 9 12]);
        pc_x = 2;
        pc_y = 3;
        % Remove outliers (e.g. noisy waveforms)
        euc_distances = sqrt(scores(:,pc_x).^2 + scores(:,pc_y).^2);
        is_outlier = isoutlier(euc_distances, 'grubbs');
        
        h = biplot(pca_ax, coeffs(:,[pc_x pc_y]),'Scores', scores(~is_outlier,[pc_x pc_y]));
        xlabel(pca_ax, 'PC2');
        ylabel(pca_ax, 'PC3');
        % Remove loading vectors
        for k = 1:size(all_wfs, 2)*2
            h(k).LineStyle = 'none';
            h(k).MarkerEdgeColor = 'none';
        end
        % Recolor markers
        for k = size(all_wfs, 2)*2+1:size(all_wfs, 2)*2+size(all_wfs(~is_outlier), 1)
            reset_idx = k - size(all_wfs, 2)*2;
            if cluster_ids(reset_idx) == 0
                h(k).MarkerEdgeColor = rgb(wf1_colors{1});
            else
                h(k).MarkerEdgeColor = rgb(wf2_colors{1});
            end
            h(k).MarkerSize = 5;
        end
    end

%% Ignore these metrics for now
% 
%     function p_corr = autocorrelogram_similarity(unit1_allSpikeTimePoints, unit2_allSpikeTimePoints, ...
%                                                  unit1_ops, unit2_ops, ...
%                                                  max_lag, bin_size)
%         % Day 1 autocorrelogram
%         bytes       = get_file_size(unit1_ops.fclean); % size in bytes of raw binary
%         nTimepoints = floor(bytes/unit1_ops.NchanTOT/2); % number of total timepoints
%         day1_tend    = double(nTimepoints) / unit1_ops.fs; % ending timepoint
%         hist_edges = 0:bin_size:day1_tend;    % 1 ms bins
%         timestamps = (double(unit1_allSpikeTimePoints) / unit1_ops.fs);
%         binned = histcounts(timestamps, hist_edges); % timestamps to binned
%         [unit1_autocor, lags] = xcorr(binned, max_lag/bin_size);
%         unit1_autocor = unit1_autocor(lags > 0);
% 
%         % Day 2 autocorrelogram
%         bytes       = get_file_size(unit2_ops.fclean); % size in bytes of raw binary
%         nTimepoints = floor(bytes/unit2_ops.NchanTOT/2); % number of total timepoints
%         day2_tend    = double(nTimepoints) / unit2_ops.fs; % ending timepoint
%         hist_edges = 0:bin_size:day2_tend;    % 1 ms bins
%         timestamps = (double(unit2_allSpikeTimePoints) / unit2_ops.fs);
%         binned = histcounts(timestamps, hist_edges); % timestamps to binned
%         [unit2_autocor, lags] = xcorr(binned, max_lag/bin_size);
%         unit2_autocor = unit2_autocor(lags > 0);
% 
%         % Correlation
%         p_corr = corrcoef(unit1_autocor, unit2_autocor);
%         p_corr = atanh(p_corr(2,1)); % Fisher transform
%     end
% 
%     function fr_diff = firing_rate_similarity(unit1_allSpikeTimePoints, ...
%                                               unit2_allSpikeTimePoints, ...
%                                               unit1_ops, unit2_ops)
%         % Unit 1
%         bytes       = get_file_size(unit1_ops.fclean); % size in bytes of raw binary
%         nTimepoints = floor(bytes/unit1_ops.NchanTOT/2); % number of total timepoints
%         max_time    = double(nTimepoints) / unit1_ops.fs; % ending timepoint
%         % Subtract periods of noise removal
%         if isfield(unit1_ops, 'concat_tranges')
%             % Read the concat_tranges to subtract the max_time from the offsets
%             trange_or_tranges = unit1_ops.concat_tranges;
%         else
%             % for non-concatenated files, just read the offset
%             trange_or_tranges = unit1_ops.trange;
%         end
% 
%         for row=1:size(trange_or_tranges, 1)
%             max_time = max_time - trange_or_tranges(row, 1);
%         end
% 
%         day1_fr = length(unit1_allSpikeTimePoints) / max_time;
% 
%         % Unit 2
%         bytes       = get_file_size(unit2_ops.fclean); % size in bytes of raw binary
%         nTimepoints = floor(bytes/unit2_ops.NchanTOT/2); % number of total timepoints
%         max_time    = double(nTimepoints) / unit2_ops.fs; % ending timepoint
%         % Subtract periods of noise removal
%         if isfield(unit2_ops, 'concat_tranges')
%             % Read the concat_tranges to subtract the max_time from the offsets
%             trange_or_tranges = unit2_ops.concat_tranges;
%         else
%             % for non-concatenated files, just read the offset
%             trange_or_tranges = unit2_ops.trange;
%         end
% 
%         for row=1:size(trange_or_tranges, 1)
%             max_time = max_time - trange_or_tranges(row, 1);
%         end
% 
%         day2_fr = length(unit2_allSpikeTimePoints) / max_time;
% 
%         fr_diff = abs(log(day1_fr) - log(day2_fr));
%     end
end