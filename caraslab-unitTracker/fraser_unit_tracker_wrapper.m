function fraser_unit_tracker_wrapper(Savedir, show_plots, bp_filter, load_previous_gwfparams, chanMap)
    % This function implements an algorithm to detect unit survival across
    % recording days from Fraser et al., 2011: unitIdentification.m
    %
    % The gist of the approach: 
    % 1. Compare units within and between consecutive days using 4 metrics: 
    %   a) autocorrelograms: computeDailyAutocorrelations.m
    %   b) baseline firing: computeDailyBaserates.m
    %   c) cross-correlograms: computeDailyCorrelations.m
    %   d) waveform similarity: computeWaveScore.m
    % 2. Define potential-same and definite-different unit groups by using a
    %   physical criterion such as different channel (original publication) or
    %   different shank (here)
    % 3. Use the 4 scores to iteratively train a quadratic discriminant classifier using
    %   partially supervised expectation-maximization to fit a mixture of Gaussian models
    %   Essentially, we use the definite-different units (hard label) to dictate
    %   the shape of one of the Gaussian distributions and the classifier
    %   computes a 4-dimensional decision boundary to segregate the other
    %
    % Many modifications were done to their code to adapt the algorithm
    % originally developed for tetrode recordings.
    % Modifications are the following:
    %   1. Waveforms on the same probe shank are treated as potential-same.
    %   2. Recording duration adapts to noise removal information at the
    %       beginning of recordings from the TDT recording system. Important in
    %       computeDailyBaserates.m
    %   3. Probe geometry is a key factor in computeWaveScore. Waveforms are
    %       compared between and across all channels (e.g. 16x16 comparisons in a 16ch shank)
    %       In these, I implemented a weighted Euclidean distance to confer more
    %       weight to comparisons between closer channels and maximum weight to
    %       same-channel comparisons

    % IMPORTANT: This code cannot be relied on 100% yet; After running, make sure to
    % inspect every waveform survival plot and make manual changes to the
    % unitSurvival.csv file (i.e., change a unit's new survival code) if you find
    % spurious grouping

    % Created by M Macedo-Lima December, 2021
    
    close all;
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
    

    channel_map = load(chanMap);
    % Create containers for the Fraser function format
    shank = {};
    unit = {};
    unitnames = {};
    spiketimes = {};
    wmean = {};
    gwfparamses = {};
    allShankChan_wfmeans = {};
    shank_channelses = {};
    % Cycle through folders
    for day_idx=1:length(unique_days)
               
        day_dir = fullfile(Savedir, datafolders{day_idx});
        
        % Load waveforms if gwfparams cannot be found in folder or if
        % load_previous_gwfparams flag is 0
        if load_previous_gwfparams
            [wf, gwfparams] = try_load_previous_gwfparams(day_dir, bp_filter, 1);
        else
            fprintf('Running wf extraction...\n')
            [wf, gwfparams] = get_waveforms_from_folder(day_dir, bp_filter, 1);
            % Save gwfparams and wf for future use
            fprintf('Saving gwfparams and wf structs to mat file\n')
            save(fullfile(day_dir, 'extracted_wfs.mat'), 'gwfparams', 'wf', '-v7.3');
        end
        
        %  To get unit ids
        split_dir = split(day_dir, filesep); 
        subj_id = split(split_dir{end-1}, '-');
        subj_id = join(subj_id(1:3), '-'); 
        subj_id = subj_id{1}; 
        recording_id = split_dir{end};
        prename = [subj_id '_' recording_id];  % this is what goes into the .txt file name
        
        day_shanks = [];
        day_units = [];
        day_unitIDs = {};
        day_spiketimes = {};
        day_wfmeans = {};
        day_gwfparams = {};
        day_allShankChan_wfmeans = {};
        day_shank_channels = {};
        for wf_idx=1:length(wf.unitIDs)
            % Normalize wfs around 0 and so that peak is at -1 or 1
            [~, day_cur_wf_mean, day_shank, allShankChan_wfmean, shank_channels] = normalize_wfs(wf, gwfparams, wf_idx);
            day_shanks(end+1) = day_shank;
            day_units(end+1) = wf_idx;
            day_unitIDs{end+1} = [prename '_cluster' int2str(wf.unitIDs(wf_idx))];
            day_spiketimes{end+1} = double(wf.allSpikeTimePoints{wf_idx}) / gwfparams.sr;
            day_wfmeans{end+1} = day_cur_wf_mean;
            day_gwfparams{end+1} = gwfparams;
            day_allShankChan_wfmeans{end+1} = allShankChan_wfmean;
            day_shank_channels{end+1} = shank_channels;
        end
        
        shank{end+1} = day_shanks;
        unit{end+1} = day_units;
        unitnames{end+1} = day_unitIDs;
        spiketimes{end+1} = day_spiketimes;
        wmean{end+1} = day_wfmeans;
        gwfparamses{end+1} = day_gwfparams;
        allShankChan_wfmeans{end+1} = day_allShankChan_wfmeans;
        shank_channelses{end+1} = day_shank_channels;
        
        clear wf gwfparams
    end
    
    % Run survival tracker
    if ~show_plots
        figure('visible', 'off');
    else
        figure;
    end
    
    %% Calculate survival
    [survival, score, corrscore, wavescore, autoscore, basescore, correlations] ...
        = unitIdentification(shank, unit, spiketimes, allShankChan_wfmeans, gwfparamses, shank_channelses, channel_map, 'plot');
    
    %% Continue to saving and plotting
    % Save variables
    mkdir(fullfile(Savedir, 'Unit tracking'));
    save(fullfile(Savedir, 'Unit tracking', 'unitTracking_output.mat'), 'survival', ...
        'score', 'corrscore', 'wavescore',  'autoscore', 'basescore', 'correlations', '-v7.3');
    
    % Save plots
    screen_size = get(0, 'ScreenSize');
    origSize = get(gcf, 'Position'); % grab original on screen size
    set(gcf, 'Position', [0 0 screen_size(3) screen_size(4)] ); %set to scren size
    set(gcf,'PaperPositionMode','auto') %set paper pos for printing
    print(fullfile(Savedir, 'Unit tracking', 'classifier_plots.pdf'), '-dpdf', '-bestfit', '-painters');
    
    % Output a survival csv with new unit IDs
    % Create unique IDs then change them when unit survives
    all_unit_n = sum(cellfun('length', unit));
    unique_unitID = reroll(unit, 1:all_unit_n);
    flatten_unitNames = {};
    flatten_unitIDs = [zeros(all_unit_n, 1)];
    TT = array2table([(1:all_unit_n)' (1:all_unit_n)'],...
            'VariableNames',{'Cluster' 'Survival_ID'});
    TT.Cluster = num2cell(TT.Cluster);

    table_row = 0;
    first_entry_flag = 1;
    for day_idx=1:(length(survival))
        
        % Kip days without units
        if isempty(unitnames{day_idx})
            continue
        end
        
        day_survival = survival{day_idx};
        day_unitIDs = unique_unitID{day_idx};
        

        [day1_surviving_idx, day2_surviving_idx] = find(day_survival == 1);
        
        day_survivingUnitIDs = day_unitIDs(day1_surviving_idx);
        
        unique_unitID{day_idx+1}(day2_surviving_idx) = day_survivingUnitIDs;
        
        % Fill first day once
        if first_entry_flag == 1
            for unit_idx=1:length(unitnames{day_idx})
                TT.Cluster(table_row + unit_idx) = unitnames{day_idx}(unit_idx);
                TT.Survival_ID(table_row + unit_idx) = unique_unitID{day_idx}(unit_idx);
            end
            table_row = table_row + unit_idx;
            first_entry_flag = 0;
        end
        
        for unit_idx=1:length(unitnames{day_idx+1})
            TT.Cluster(table_row + unit_idx) = unitnames{day_idx+1}(unit_idx);
            TT.Survival_ID(table_row + unit_idx) = unique_unitID{day_idx+1}(unit_idx);
        end
        table_row = table_row + unit_idx;
    end
    
    % create a table with unit names and IDs
    writetable(TT, fullfile(Savedir, 'Unit tracking', [subj_id '_unitSurvival.csv']));        
    
    % Plot waveforms overlain; probe geometry is respected
    unique_IDs = unique(TT.Survival_ID);
    rerolled_IDs = reroll(unit, TT.Survival_ID);
    for id_idx=1:length(unique_IDs)
        cur_ID = unique_IDs(id_idx);
        if sum(TT.Survival_ID == cur_ID) > 1  % unit survived past 1 day
            figure;
            hold on
            day_colormap = lines;
            cax_list=[];
            surviving_indices = find(TT.Survival_ID == cur_ID);
            for surviving_indices_iid=1:length(surviving_indices)
                plotted_days = [];
                for day_iid=1:length(rerolled_IDs)
                    cur_wf_mask = rerolled_IDs{day_iid} == surviving_indices(surviving_indices_iid);
                    if sum(cur_wf_mask) == 0
                        continue
                    end
                    plotted_days = [plotted_days; day_iid];
                    gwfparams = gwfparamses{day_iid}{cur_wf_mask};
                    shank_channels = shank_channelses{day_iid}{cur_wf_mask};

                    % Rescale x and y to 0-1
                    scaled_channel_positions = [channel_map.xcoords channel_map.ycoords];
                    scaled_channel_positions(:,1) = rescale(scaled_channel_positions(:,1));
                    scaled_channel_positions(:,2) = rescale(scaled_channel_positions(:,2));

                    for shank_channel_idx=1:numel(shank_channels)
                        cur_chanMap_ch = shank_channels(shank_channel_idx);

                        % Grab x and y positions; rescale them to look closer
                        % together; might need to be tweaked for each probe type
                        x_offset = scaled_channel_positions(cur_chanMap_ch, 1)*120;
                        y_offset = scaled_channel_positions(cur_chanMap_ch, 2)*5000;

                        % Squeeze out and store raw waveforms and averages
                        snip_points = 20;  % snip some points at the end of template

                        cur_wf_mean = allShankChan_wfmeans{day_iid}{cur_wf_mask}(1:end-snip_points, shank_channel_idx);
                        
                        if any(isnan(cur_wf_mean))
                            continue
                        end
                        
                        % 10x upsample mean waveform with spline interpolation
                        samplingRateIncrease = 10;
                        amplitude_gain = 50;
                        newXSamplePoints = linspace(1, length(cur_wf_mean), length(cur_wf_mean) * samplingRateIncrease);
                        cur_wf_mean_upsample = spline(1:length(cur_wf_mean), cur_wf_mean, newXSamplePoints) * amplitude_gain;

                        x_time = linspace(0, length(cur_wf_mean) / gwfparams.sr, length(cur_wf_mean));
                        x_time = x_time * 1000; % in ms

                        upsampled_x_time = linspace(0, length(cur_wf_mean_upsample) / gwfparams.sr / samplingRateIncrease, length(cur_wf_mean_upsample));
                        upsampled_x_time = upsampled_x_time*1000;

                        cax = plot(upsampled_x_time + x_offset, cur_wf_mean_upsample + y_offset, 'linewidth', 2, 'Color', day_colormap(length(plotted_days), :));
                        if shank_channel_idx == 1
                            cax_list = [cax_list; cax];
                        end
                    end
                end

                axis off;
                legend(cax_list, TT(TT.Survival_ID == cur_ID,:).Cluster, 'Interpreter', 'none', 'Position', [0.1 0.1 0.15 0.1], 'Box', 'off');
                % Save plot
                screen_size = get(0, 'ScreenSize');
                set(gcf, 'Position', [0 0 screen_size(3) screen_size(4)] ); %set to scren size
                set(gcf,'PaperPositionMode','auto') %set paper pos for printing
                print(fullfile(Savedir, 'Unit tracking', ...
                    ['survival_wfs_' TT(TT.Survival_ID == cur_ID,:).Cluster{1} '_ID' ...
                    int2str(surviving_indices(surviving_indices_iid))]), '-dpdf', '-bestfit', '-painters');
            end
        end
    end
end