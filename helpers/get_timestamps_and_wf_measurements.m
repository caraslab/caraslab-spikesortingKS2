function get_timestamps_and_wf_measurements(Savedir, show_plots, bp_filter)
    % This function retrieves timestamps and waveforms from phy files
    % Outputs are .txt files with timestamps, .csv and .pdf files with waveform
    % measurements and plots    
    % Adapted from Cortexlab repository
    % Patched by M Macedo-Lima 9/8/20
    % New in this patch:
    %   - Call mode (selecting directories)
    %   - show_plots (0 or 1) to visualize waveforms
    %   - bp_filter (0 or 1) to refilter before measuring units. Since KS2
    %       uses a 150 high-pass, using this only could give you noisy units.
    %       The new filter is a 300-6000 Hz bandpass. A new file will be
    %       output, called *CLEAN_300Hz.dat; can be deleted after this
    
    % Select folders to run
 
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
        

        cur_path.name = datafolders{i};
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
        fprintf('Outputting timestamps and measuring waveforms for: %s\n', cur_savedir)

        % Delete old timestamp files in folder
        old_files = dir([cur_savedir '/*cluster*txt']);
        for i=1:numel(old_files)
            delete(fullfile(old_files(i).folder, old_files(i).name));
        end
        
        old_files = dir(fullfile(cur_savedir, 'CSV files', '*_waveforms.csv'));
        for dummy_idx=1:numel(old_files)
            delete(fullfile(old_files(dummy_idx).folder, old_files(dummy_idx).name));
        end
        
        %% Define I/O and waveform parameters
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
        
        % Store in ops for future use
        ops.subj_id = subj_id;
        ops.recording_id = recording_id;
        ops.subj_recording_id = prename;
        
        
        gwfparams.rawDir = cur_savedir;
        gwfparams.sr = ops.fs;
        gwfparams.nCh = ops.NchanTOT; % Number of channels that were streamed to disk in .dat file
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
        
        gwfparams.dataType = 'int16'; % Data type of .dat file (this should be BP filtered)

        gwfparams.wfWin = [round(-(0.001*gwfparams.sr)) round(0.003*gwfparams.sr)]; % Number of samples before and after spiketime to include in waveform
        gwfparams.nWf = 2000; % Max number of waveforms per unit to pull out for averaging
        gwfparams.spikeTimes = readNPY(fullfile(gwfparams.dataDir, 'spike_times.npy')); % Vector of cluster spike times (in samples) same length as .spikeClusters
        gwfparams.spikeClusters = readNPY(fullfile(gwfparams.dataDir, 'spike_clusters.npy')); % Vector of cluster IDs (Phy nomenclature)   same length as .spikeTimes
        gwfparams.channelShanks = readNPY(fullfile(gwfparams.dataDir, 'channel_shanks.npy')); % Vector of cluster shanks

        gwfparams.chanMap = readNPY(fullfile(gwfparams.dataDir, 'channel_map.npy')); % this is important in esp if you got rid of files. 
        try
            gwfparams.cluster_quality = tdfread(fullfile(gwfparams.dataDir, 'cluster_info.tsv'));
        catch ME
            if strcmp(ME.identifier, 'MATLAB:load:couldNotReadFile')
                fprintf('\nFile not found\n')
                continue
            else
                fprintf(ME.identifier)  % file not found has no identifier?? C'mon MatLab...
                fprintf([ME.message '\n'])
                continue  % Continue here instead of break because I don't know how to catch 'file not found' exception; maybe using ME.message?
            end
        end

        % Get good and mua for measuring
        gwfparams.good_clusters = gwfparams.cluster_quality.cluster_id(gwfparams.cluster_quality.group(:,1)=='g' | gwfparams.cluster_quality.group(:,1)=='m');

        %% Get waveforms from .dat
        wf = getWaveForms(gwfparams, bp_filter);  

        
        %% Save gwfparams and wf for future use
        fprintf('Saving gwfparams and wf structs to mat file\n')
        save(fullfile(gwfparams.dataDir, 'extracted_wfs.mat'), 'gwfparams', 'wf', '-v7.3');
        
        %% Store some other stuff
        good_cluster_idx = wf.unitIDs; % store their Phy IDs

        % Create a folder called CSV files within saved directory, because
        % older versions of Phy crash if random .csv files are present
        % Contains other file types too, but oh well...
        if (exist(fullfile(gwfparams.dataDir, 'CSV files')) == 0)
            mkdir(fullfile(gwfparams.dataDir, 'CSV files')); 
        end

        % Save averages to .mat file
        % Useful in case you want to replot using BETTER software, i.e. Python
        save(fullfile(gwfparams.dataDir, 'CSV files', [prename '_waveform_averages.mat']), 'wf', '-v7.3');

        %% OUTPUT from getWaveforms looks like this:
        % wf.unitIDs                               % [nClu,1]            List of cluster IDs; defines order used in all wf.* variables
        % wf.spikeTimeKeeps                        % [nClu,nWf]          Which spike times were used for the waveforms
        % wf.waveForms                             % [nClu,nWf,nCh,nSWf] Individual waveforms
        % wf.waveFormsMean                         % [nClu,nCh,nSWf]     Average of all waveforms (per channel)
        %                                          % nClu: number of different clusters in .spikeClusters
        %                                          % nSWf: number of samples per waveform

        %% Measure waveform averages
        % Prealocate some variables
        ptp_durations = zeros(length(gwfparams.good_clusters), 1);
        ptp_ratios = zeros(length(gwfparams.good_clusters), 1);
        repolarization_durations = zeros(length(gwfparams.good_clusters), 1);
        best_channels_csv = zeros(length(gwfparams.good_clusters), 1);
        shanks = zeros(length(gwfparams.good_clusters), 1);
        cluster_quality = zeros(length(gwfparams.good_clusters), 1);  % zeros will be MU, ones will be SU
        
        if show_plots
            figure('Name', prename);
        else
            figure('Name', prename, 'visible', 'off');
        end
        for wf_idx=1:length(good_cluster_idx)
            cluster_phy_id = good_cluster_idx(wf_idx);

            % Get channel with highest amplitude for this unit
            % This usually works but sometimes phy fails at detecting the best
            % channel, so let's do it manually
    %         best_channel = gwfparams.cluster_quality.ch(...
    %             gwfparams.cluster_quality.cluster_id == cluster_phy_id);

            temp_wfs = wf.waveFormsMean(wf_idx, :,:);
            temp_wfs = squeeze(abs(temp_wfs));

            [max_amplitude, ~] = max(temp_wfs, [], 2);
            
            clear temp_wfs

            [~, best_channel_order] = max(max_amplitude);

            % Grab best channel index
            best_channel_0in = gwfparams.chanMap(best_channel_order);
            
            % Store channel and shank info
            best_channels_csv(wf_idx) = best_channel_0in + 1;
            shanks(wf_idx) = gwfparams.cluster_quality.sh(gwfparams.cluster_quality.cluster_id == cluster_phy_id);  
            
            % Grab best channel index
            best_channel_idx = find(gwfparams.chanMap == best_channel_0in);
            
            % Store cluster quality
            cluster_quality_char = gwfparams.cluster_quality.group(gwfparams.cluster_quality.cluster_id == cluster_phy_id);
            if (cluster_quality_char == 'g')
            	cluster_quality(wf_idx) = 1;
            end
            
            % Squeeze out and store raw waveforms and averages
            cur_wfs = wf.waveForms(wf_idx, :, best_channel_idx,:);
            cur_wfs = squeeze(cur_wfs);
            cur_wf_mean = wf.waveFormsMean(wf_idx, best_channel_idx,:);
            cur_wf_mean = squeeze(cur_wf_mean);
            
            if isnan(mean(cur_wf_mean))
                continue  % sometimes waveforms are NaN which is weird...
            end

            % Recenter average at zero
            cur_wf_mean = cur_wf_mean - mean(cur_wf_mean(1:5));

            %% Write spike times to txt file
            fprintf('Outputting spike times to txt file\n')
            fileID = fopen(fullfile(gwfparams.dataDir, [prename '_cluster' int2str(cluster_phy_id) '.txt']), 'w');
            fprintf(fileID, '%.6f\n', (double(wf.allSpikeTimePoints{wf_idx}) / gwfparams.sr));
            fclose(fileID);
            
            %% Write spike waveforms to csv file
            % Delete old files first
            fprintf('Outputting spike waveforms to csv file\n')
            writetable(array2table(cur_wfs), fullfile(gwfparams.dataDir, 'CSV files', ...
                [prename '_cluster' int2str(cluster_phy_id) '_waveforms.csv']));
            
            %% Continue calculations
            % 10x upsample mean waveform with spline interpolation
            samplingRateIncrease = 10;
            newXSamplePoints = linspace(1, length(cur_wf_mean), length(cur_wf_mean) * samplingRateIncrease);
            cur_wf_mean_upsample = spline(1:length(cur_wf_mean), cur_wf_mean, newXSamplePoints)';

            % peak to peak
            abs_wf = abs(cur_wf_mean_upsample);
            peak_idx = find(abs_wf == max(abs_wf));

            % Normalize amplitude so that peak is at -1 or 1
            cur_wf_mean_upsample = cur_wf_mean_upsample ./ abs(cur_wf_mean_upsample(peak_idx));

            % Find trough
            % Could be replaced by  trough_idx = find(cur_wf_mean_upsample == max(abs(cur_wf_mean_upsample(peak_idx:end))));
            % But if average has a third peak of different polarity, the above might
            % behave strangely...
            if cur_wf_mean_upsample(peak_idx) < 0 % downward deflection
                trough_idx = find(cur_wf_mean_upsample == max(cur_wf_mean_upsample(peak_idx:end)));
            else   % upward deflection
                trough_idx = find(cur_wf_mean_upsample == min(cur_wf_mean_upsample(peak_idx:end)));
            end

            % Grab PTP duration
            ptp_duration = (trough_idx - peak_idx) / gwfparams.sr * 1000 / samplingRateIncrease; %in ms

            assert(ptp_duration > 0)  % break if something goes wrong

            % Trainito et al., 2019 says: 
            % we excluded waveforms that
            % satisfied any of three criteria for atypical shape: 
            % (1) the amplitude of the main trough was smaller than the subsequent positive peak (n = 41), 
            % (2) the trace was noisy, defined as > = 6 local maxima of magnitude > = 0.01 (n = 38), 
            % (3) there was one or more local maxima in the period between the main trough and the subsequent peak (n = 35).
            % I'll implement number 1 here as a safety mechanism, replace
            % below with an  'if true then skip unit' statement if you keep getting this error
            assert(abs(cur_wf_mean_upsample(peak_idx)) > abs(cur_wf_mean_upsample(trough_idx)));

            % PTP ratio
            ptp_ratio = abs(cur_wf_mean_upsample(peak_idx) / cur_wf_mean_upsample(trough_idx));

            % Repolarization duration, a measure used by Trainito et al. 2019
            % to cluster cortical neurons in monkey
            % Defined as the time between the trough and the inflection point
            % (when second derivative crosses 0)
            % Derivative is too jittery around 0. Smoothing it by ~1 ms (i.e. 0.001 * samplingRateIncrease * ops.fs)
            % seems to correct for this. Might need to be adjusted for
            % different sampling rates
            second_dv_dt = smooth(gradient(gradient(cur_wf_mean_upsample)), round(0.001 * samplingRateIncrease * ops.fs), 'loess');

            % Below debug code is for checking the derivative plots; make sure to use
            % breakpoints
    %         figure
    %         hold on
    %         plot(cur_wf_mean_upsample)
    %         plot(second_dv_dt*1000)
    %         plot(gradient(gradient(cur_wf_mean_upsample))*500)
    %         yline(0)        

            % find zero crossings
            upcross = find(second_dv_dt(1:end-1) <= 0 & second_dv_dt(2:end) > 0);
            downcross = find(second_dv_dt(1:end-1) >= 0 & second_dv_dt(2:end) < 0);

            all_crosses = sort([upcross; downcross]);

            % find inflection after the trough
            try
                dummy_next_points = all_crosses(trough_idx - all_crosses < 0);
                inflection_point = dummy_next_points(1);
                repolarization_duration = (inflection_point - trough_idx) / gwfparams.sr * 1000 / samplingRateIncrease;
                repolarization_durations(wf_idx) = repolarization_duration;
            catch ME
                fprintf('Repolarization calculation failed.')
            end

            % Store measurements
            ptp_durations(wf_idx) = ptp_duration;
            ptp_ratios(wf_idx) = ptp_ratio;

            %% PLOT
            subplot(ceil(length(good_cluster_idx) / 4), 4, wf_idx);
            x_time = linspace(0, length(cur_wf_mean) / gwfparams.sr, length(cur_wf_mean));
            x_time = x_time * 1000; % in ms

            upsampled_x_time = linspace(0, length(cur_wf_mean_upsample) / gwfparams.sr / samplingRateIncrease, length(cur_wf_mean_upsample));
            upsampled_x_time = upsampled_x_time*1000;

            abs_cur_wf_mean = abs(cur_wf_mean);
            cur_wf_mean_peak_idx = abs_cur_wf_mean == max(abs_cur_wf_mean);
            norm_peak_value = abs(cur_wf_mean(cur_wf_mean_peak_idx));
            % Normalize amplitude so that peak is at -1
            cur_wf_mean = cur_wf_mean ./ norm_peak_value;
            
            title(['Cluster' num2str(cluster_phy_id)]);
            
            hold on

            for spike_idx=1:min(size(cur_wfs, 1), 500)  % plot a max of 500 spikes
                cur_wf = cur_wfs(spike_idx,:) - mean(cur_wfs(spike_idx,1:5));
                if isnan(mean(cur_wfs(spike_idx, 1:5)))
                    continue  % sometimes waveforms are NaN which is weird...
                end
                % Normalize amplitude of individual spike using the average
                cur_wf = cur_wf ./ norm_peak_value;

                cax = plot(x_time, cur_wf, 'black'); % Plot 0-centered traces
                cax.Color(4)=0.1;
            end


            plot(upsampled_x_time, cur_wf_mean_upsample, 'red', 'linewidth', 1.5);
            xlim([min(x_time) max(x_time)])
            ylims = 5; 
            ylim([-ylims ylims])  

            % Vertical peak lines
            plot([peak_idx / gwfparams.sr * 1000 / samplingRateIncrease  peak_idx / gwfparams.sr * 1000 / samplingRateIncrease], ...
                 [-ylims, ylims], '--')

            plot([trough_idx / gwfparams.sr * 1000 / samplingRateIncrease trough_idx / gwfparams.sr * 1000 / samplingRateIncrease], ...
                 [-ylims, ylims], '--')

            % Horizontal measuring line
            plot([-ylims, ylims], ...
                 [cur_wf_mean_upsample(peak_idx) cur_wf_mean_upsample(peak_idx)], '-')

             plot([-ylims, ylims], ...
                 [cur_wf_mean_upsample(trough_idx) cur_wf_mean_upsample(trough_idx)], '-')

            % Inflection point
            plot([inflection_point / gwfparams.sr * 1000 / samplingRateIncrease inflection_point / gwfparams.sr * 1000 / samplingRateIncrease],...
            [-ylims, ylims], '-')

            duration_txt = ['PTP = ' num2str(ptp_duration, 2) ' ms'];
            ratio_txt = ['Ratio = ' num2str(ptp_ratio, 2)];
            inflection_txt = ['Inflection = ' num2str(repolarization_duration, 3) ' ms'];
            text(2, double(-ylims/2), duration_txt, 'fontsize',10);
            text(2, double(-ylims/2-ylims/5), ratio_txt, 'fontsize',10);
            text(2, double(-ylims/2-2*ylims/5), inflection_txt, 'fontsize',10);

        end

        screen_size = get(0, 'ScreenSize');

        origSize = get(gcf, 'Position'); % grab original on screen size
        set(gcf, 'Position', [0 0 screen_size(3) screen_size(4)] ); %set to scren size
        set(gcf,'PaperPositionMode','auto') %set paper pos for printing
        print([gwfparams.dataDir '\CSV files\' prename '_waveform_measurements'], '-dpdf', '-bestfit', '-painters');

        %% write wf measurements
        TT = array2table([gwfparams.good_clusters best_channels_csv shanks cluster_quality ptp_durations ptp_ratios repolarization_durations],...
            'VariableNames',{'Cluster' 'Best_channel' 'Shank' 'Cluster_quality' 'PTP_duration_ms' 'PTP_ratio' 'Repolarization_duration_ms'});
        % Change code into words
        TT.Cluster_quality = num2cell(TT.Cluster_quality);
        
        for i=1:length(TT.Cluster_quality)
            if TT.Cluster_quality{i} == 1
                TT.Cluster_quality(i) = {'good'};
            else
                TT.Cluster_quality(i) = {'mua'};
            end
        end
        writetable(TT, fullfile(gwfparams.dataDir, 'CSV files', [prename '_waveform_measurements.csv']));
                
        tEnd = toc(t0);
        fprintf('Done in: %d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));
    end
end
