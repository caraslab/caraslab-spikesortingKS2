function pca_compare_waveforms(Savedir, show_plots, bp_filter)
    % This function loops through recording folders and compares waveforms between
    % consecutive days if they occured on the same shank. The output is a plot 
    % showing the compared waveforms and the correlations between the PC scores 
    % up to the 3rd PC. It also outputs a csv file per comparison with the PC scores.
    % getWaveforms is from the Cortexlab github repository
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
        clear ops
        close all
        
        day2_idx=day1_idx+1;
        
        day1_dir = fullfile(Savedir, datafolders{day1_idx});
        day2_dir = fullfile(Savedir, datafolders{day2_idx});

        output_dir = [num2str(unique_days(day1_idx)) 'vs' num2str(unique_days(day2_idx)) '_pca_comparison'];
        output_file_name = [num2str(unique_days(day1_idx)) 'vs' num2str(unique_days(day2_idx))];
        mkdir(fullfile(Savedir, output_dir));
        
        % Load waveforms
        [day1_wf, day1_gwfparams] = get_waveforms_from_folder(day1_dir, bp_filter);
        [day2_wf, day2_gwfparams] = get_waveforms_from_folder(day2_dir, bp_filter);
        
        % Snip waveforms from 0.5 to 2 ms to increase PCA discriminability
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
                    
                    % Concatenate wfs in a new variable for PCA
                    all_wfs = [day1_cur_wfs; day2_cur_wfs];
                    % Keep track of indices
                    original_inds = [1:size(day1_cur_wfs, 1) 1:size(day2_cur_wfs, 1)];
                    % Keep track of clusters
                    cluster_ids = [ zeros(size(day1_cur_wfs, 1), 1); ones(size(day2_cur_wfs, 1),1) ];
                    
                    % Start plotting
                    figure;
                    
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
                    
                    wf1_colors = {'CornflowerBlue', 'MediumBlue'};
                    wf2_colors = {'LightCoral', 'FireBrick'};
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
                    pca_ax = subplot(4, 3, [4 7 10]);
                    h = biplot(pca_ax, coeffs(:,1:2),'Scores', scores(:,1:2));
                    xlabel(pca_ax, 'PC1');
                    ylabel(pca_ax, 'PC2');
                    % Remove loading vectors
                    for k = 1:size(all_wfs, 2)*2
                        h(k).LineStyle = 'none';
                        h(k).MarkerEdgeColor = 'none';
                    end
                    % Recolor markers
                    for k = size(all_wfs, 2)*2+1:size(all_wfs, 2)*2+size(all_wfs, 1)
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
                    h = biplot(pca_ax, coeffs(:,[1 3]),'Scores', scores(:,[1 3]));
                    xlabel(pca_ax, 'PC1');
                    ylabel(pca_ax, 'PC3');
                    % Remove loading vectors
                    for k = 1:size(all_wfs, 2)*2
                        h(k).LineStyle = 'none';
                        h(k).MarkerEdgeColor = 'none';
                    end
                    % Recolor markers
                    for k = size(all_wfs, 2)*2+1:size(all_wfs, 2)*2+size(all_wfs, 1)
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
                    h = biplot(pca_ax, coeffs(:,[2 3]),'Scores', scores(:,[2 3]));
                    xlabel(pca_ax, 'PC2');
                    ylabel(pca_ax, 'PC3');
                    % Remove loading vectors
                    for k = 1:size(all_wfs, 2)*2
                        h(k).LineStyle = 'none';
                        h(k).MarkerEdgeColor = 'none';
                    end
                    % Recolor markers
                    for k = size(all_wfs, 2)*2+1:size(all_wfs, 2)*2+size(all_wfs, 1)
                        reset_idx = k - size(all_wfs, 2)*2;
                        if cluster_ids(reset_idx) == 0
                            h(k).MarkerEdgeColor = rgb(wf1_colors{1});
                        else
                            h(k).MarkerEdgeColor = rgb(wf2_colors{1});
                        end
                        h(k).MarkerSize = 5;
                    end
                    
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
                end
            end
        
        end
    end
    
    function [cur_wfs, cur_wf_mean, shank] = normalize_wfs(wf, gwfparams, wf_idx)
        % This function levels wfs to be around 0V and rescales so that
        % their peak is at -1 or 1
        
        cluster_phy_id = wf.unitIDs(wf_idx);
        % Get channel with highest amplitude for this unit
        
        % This usually works but sometimes phy fails at detecting the best
        % channel, so let's do it manually
%         best_channel = gwfparams.cluster_quality.ch(...
%             gwfparams.cluster_quality.cluster_id == cluster_phy_id);

        temp_wfs = wf.waveFormsMean(wf_idx, :,:);
        temp_wfs = abs(temp_wfs);
        
        [max_amplitude, ~] = max(temp_wfs, [], 3);
        clear temp_wfs
        
        [~, best_channel] = max(max_amplitude);
        best_channel_0in = best_channel - 1;

        % Grab best channel index
        best_channel_idx = find(gwfparams.chanMap == best_channel_0in);
        % Get shank
        shank = gwfparams.cluster_quality.sh(...
            gwfparams.cluster_quality.cluster_id == cluster_phy_id);  

        % Get waveforms at best channel
        cur_wfs = wf.waveForms(wf_idx, :, best_channel_idx,:);
        cur_wfs = squeeze(cur_wfs);
        % Find and remove nan wfs
        cur_wfs(any(isnan(cur_wfs), 2),:) = [];	

        % Normalize amplitude so that peak is at -1 or 1
        cur_wf_mean = wf.waveFormsMean(wf_idx, best_channel_idx,:);
        cur_wf_mean = squeeze(cur_wf_mean);
        
        % Recenter and rescale wfs
        cur_wf_mean = cur_wf_mean - mean(cur_wf_mean(1:5));
        abs_wf = abs(cur_wf_mean);
        peak_idx = find(abs_wf == max(abs_wf));
        peak_value = abs(cur_wf_mean(peak_idx));
        
        if length(peak_value) > 1
            peak_value = peak_value(1);
        end
       
        
        cur_wf_mean = cur_wf_mean ./ peak_value;

        cur_wfs = cur_wfs ./ peak_value;
    end
    
    function [wf, gwfparams, ops] = get_waveforms_from_folder(cur_savedir, bp_filter)
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
        gwfparams.fileName = dir(ops.fclean).name; % .dat file containing the raw used for sorting
        gwfparams.dataType = 'int16'; % Data type of .dat file (this should be BP filtered)

        gwfparams.wfWin = [round(-(0.001*gwfparams.sr)) round(0.003*gwfparams.sr)]; % Number of samples before and after spiketime to include in waveform
        gwfparams.nWf = 500; % Max number of waveforms per unit to pull out for averaging
        gwfparams.spikeTimes = readNPY(fullfile(gwfparams.dataDir, 'spike_times.npy')); % Vector of cluster spike times (in samples) same length as .spikeClusters
        gwfparams.spikeClusters = readNPY(fullfile(gwfparams.dataDir, 'spike_clusters.npy')); % Vector of cluster IDs (Phy nomenclature)   same length as .spikeTimes
        gwfparams.channelShanks = readNPY(fullfile(gwfparams.dataDir, 'channel_shanks.npy')); % Vector of cluster shanks

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

        % Get good only for measuring
%         gwfparams.good_clusters = gwfparams.cluster_quality.cluster_id(gwfparams.cluster_quality.group(:,1)=='g' | gwfparams.cluster_quality.group(:,1)=='m');
        gwfparams.good_clusters = gwfparams.cluster_quality.cluster_id(gwfparams.cluster_quality.group(:,1)=='g');

        % Get waveforms from .dat
        wf = getWaveForms(gwfparams, bp_filter);  
    end
end