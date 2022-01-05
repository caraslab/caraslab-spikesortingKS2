function cluster_quality_metrics(Savedir, show_plots, bp_filter, load_previous_gwfparams)
% This function runs 3 quality control metrics on the curated clusters:
% 1. ISI violation false positive rate: how many false positive spikes in a
%   cluster.
% 2. Fraction of spikes missing: based on the probability distribution of 
%   spikes detected for a unit, how many are estimated to be missing? 
% 3. Presence ratio: for how much of the recording is a unit present? The
%   recording time is divided in 100 bins and the fraction of bins with at
%   least one spike present is calculated.
% Adapted from the Allen Brain Institute github

% By default, the AllenSDK applies filters so only units above a set of thresholds are returned.
% 
% The default filter values are as follows:
% 
% isi_fprate < 0.5
% amplitude_cutoff < 0.1
% presence_ratio > 0.9

% Inspired by the AllenInstitude github code. 
% Parts of this code are adapted from Cortexlab repository
% Implemented by M Macedo-Lima, December 2020


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
        
        fprintf('Calculating quality control metrics for: %s\n', cur_savedir)
        
        % Define output name for cluster; this is specific for my
        % file-naming convention and should be tweaked
        split_dir = split(cur_savedir, '/'); 
        subj_id = split(split_dir{end-1}, '-');
        subj_id = join(subj_id(1:3), '-'); 
        subj_id = subj_id{1}; 
        recording_id = split_dir{end};
        prename = [subj_id '_' recording_id];  % this is what goes into the .txt file name
        
        if load_previous_gwfparams
            [wf, gwfparams] = try_load_previous_gwfparams(cur_savedir, bp_filter, 0);
        else
            fprintf('Running wf extraction...\n')
            [wf, gwfparams] = get_waveforms_from_folder(cur_savedir, bp_filter, 0);
            % Save gwfparams and wf for future use
            fprintf('Saving gwfparams and wf structs to mat file\n')
            save(fullfile(day_dir, 'extracted_wfs.mat'), 'gwfparams', 'wf', '-v7.3');
        end
        
        %% Get waveforms from .dat
        fprintf('Sampling waveforms and calculating averages...')

        good_cluster_idx = wf.unitIDs; % store their Phy IDs

        % Create a folder called CSV files within saved directory, because
        % older versions of Phy crash if random .csv files are present
        % Contains other file types too, but oh well...
        if (exist(fullfile(gwfparams.dataDir, 'CSV files')) == 0)
            mkdir(fullfile(gwfparams.dataDir, 'CSV files')); 
        end

        %% OUTPUT from getWaveforms looks like this:
        % wf.unitIDs                               % [nClu,1]            List of cluster IDs; defines order used in all wf.* variables
        % wf.spikeTimeKeeps                        % [nClu,nWf]          Which spike times were used for the waveforms
        % wf.waveForms                             % [nClu,nWf,nCh,nSWf] Individual waveforms
        % wf.waveFormsMean                         % [nClu,nCh,nSWf]     Average of all waveforms (per channel)
        %                                          % nClu: number of different clusters in .spikeClusters
        %                                          % nSWf: number of samples per waveform

        %% Measure waveform averages
        % Prealocate some variables
        best_channels_csv = zeros(length(gwfparams.good_clusters), 1);
        shanks = zeros(length(gwfparams.good_clusters), 1);
        cluster_quality = zeros(length(gwfparams.good_clusters), 1);  % zeros will be MU, ones will be SU
        fpRate_list = zeros(length(gwfparams.good_clusters), 1);
        fraction_missing_list = zeros(length(gwfparams.good_clusters), 1);
        presence_ratio_list = zeros(length(gwfparams.good_clusters), 1);
        if show_plots
            figure('Name', prename);
        else
            figure('Name', prename, 'visible', 'off');
        end
        for wf_idx=1:length(good_cluster_idx)
            cluster_phy_id = good_cluster_idx(wf_idx);

            % Get channel with highest amplitude for this unit
            % best_channel = gwfparams.cluster_quality.ch(...
            %   gwfparams.cluster_quality.cluster_id == cluster_phy_id);
            
            % The above usually works but sometimes phy fails at detecting the best
            % channel, so let's do it manually
            temp_wfs = wf.waveFormsMean(wf_idx, :,:);
            temp_wfs = abs(temp_wfs);

            [max_amplitude, ~] = max(temp_wfs, [], 3);
            clear temp_wfs

            [~, best_channel] = max(max_amplitude);
            best_channel_0in = best_channel - 1;

            % Grab best channel index
            best_channel_idx = find(gwfparams.chanMap == best_channel_0in);
            
            % Store channel and shank info
            best_channels_csv(wf_idx) = best_channel_0in + 1;
            shanks(wf_idx) = gwfparams.cluster_quality.sh(gwfparams.cluster_quality.cluster_id == cluster_phy_id);  
            
            % Store cluster quality
            cluster_quality_char = gwfparams.cluster_quality.group(gwfparams.cluster_quality.cluster_id == cluster_phy_id);
            if (cluster_quality_char == 'g')
            	cluster_quality(wf_idx) = 1;
            end
            
            % Squeeze out and store raw waveforms and averages
            cur_wfs = wf.waveForms(wf_idx, :, best_channel_idx,:);
            cur_wfs = squeeze(cur_wfs);
            
            %% ISI violations
%             """Calculate ISI violations for a spike train.
% 
%             Based on metric described in Hill et al. (2011) J Neurosci 31: 8699-8705
            
            % Declare some parameters first
            all_spike_times = double(wf.allSpikeTimePoints{wf_idx}) / gwfparams.sr;
            min_time = ops.trange(1);
            
            filenamestruct = dir(ops.fclean); % .dat file containing the raw used for sorting
            if isempty(filenamestruct)
                % try to find the fclean file name within the cur_savedir
                split_fclean_path = split(ops.fclean, '/');
                fclean = split_fclean_path{end};
                filenamestruct = dir(fullfile(cur_savedir, fclean)); % .dat file containing the raw used for sorting
                
                % If still not found, try adding the 300hz.dat suffix, else
                % break
                if isempty(filenamestruct)
                    fclean = split_fclean_path{end}(1:end-4);
                    filenamestruct = dir(fullfile(cur_savedir, [fclean '300hz.dat']));
                else
                    fprintf('DAT file not found. Check if it is in folder')
                    break
                end
            end
            
            % This is just to get the recording duration. You can also find
            % this from rez.mat outputted from kilosort. Might be more
            % reliable...
            dataTypeNBytes = numel(typecast(cast(0, gwfparams.dataType), 'uint8')); % determine number of bytes per sample
            nSamp = filenamestruct.bytes/(gwfparams.nCh*dataTypeNBytes);  % Number of samples per channel
            max_time = nSamp / gwfparams.sr;
            
            % Modify these if you'd like
            %     isi_threshold : threshold for isi violation (i.e. refractory period)
            %     min_isi : threshold for duplicate spikes (if run after the postprocessing step, it shouldn't detect anything)
            isi_threshold = 0.0015;
            min_isi = 0.00015;
            
            % Then run
            [fpRate, ~] = isi_violations(all_spike_times, min_time, max_time, isi_threshold, min_isi, ops);
            fpRate_list(wf_idx) = fpRate;
            
            % TODO: Reconstruct autocorrelograms
            
            %% Amplitude cutoff
%             """ Calculate approximate fraction of spikes missing from a distribution of amplitudes
%             Assumes the amplitude histogram is symmetric (not valid in the presence of drift)
%             Inspired by metric described in Hill et al. (2011) J Neurosci 31: 8699-8705
             [fraction_missing, h, b, G] = amplitude_cutoff(cur_wfs);
             % Store relevant variable
             fraction_missing_list(wf_idx) = fraction_missing;
             % Plots
             subplot(ceil(length(good_cluster_idx) / 4), 4, wf_idx);
             plot(b(1:end-1), h)
             hold on
             plot([b(G(1)) b(G(1))], [0, max(h)], '--')
%              fraction_txt = ['Fraction missing = ' num2str(fraction_missing, 3)];
%              text(b(G(1)), 0, fraction_txt, 'fontsize',10);
             title(['Cluster' num2str(cluster_phy_id)]);
            
             %% Presence ratio
%              Presence ratio is not a standard metric in the field, 
%              but it's straightforward to calculate and is an easy way to 
%              identify incomplete units. It measures the fraction of time 
%              during a session in which a unit is spiking, and ranges from 
%              0 to 0.99 (an off-by-one error in the calculation ensures that it will never reach 1.0).
             pr = presence_ratio(all_spike_times, min_time, max_time, 100, ops);
             presence_ratio_list(wf_idx) = pr;
             
             % TODO: Plot histograms

        end
        
        %% Save fraction missing plots
        screen_size = get(0, 'ScreenSize');

        origSize = get(gcf, 'Position'); % grab original on screen size
        set(gcf, 'Position', [0 0 screen_size(3) screen_size(4)] ); %set to scren size
        set(gcf,'PaperPositionMode','auto') %set paper pos for printing
        print([gwfparams.dataDir '\CSV files\' prename 'fraction_missing'], '-dpdf', '-bestfit', '-painters');
        
        %% write csv
        TT = array2table([gwfparams.good_clusters best_channels_csv shanks ...
            cluster_quality fpRate_list fraction_missing_list presence_ratio_list],...
            'VariableNames',{'Cluster' 'Best_channel' 'Shank' ...
            'Cluster_quality' 'ISI_FPRate' 'Fraction_missing', 'Presence_ratio'});
        
        % Change code into words
        TT.Cluster_quality = num2cell(TT.Cluster_quality);
        
        for i=1:length(TT.Cluster_quality)
            if TT.Cluster_quality{i} == 1
                TT.Cluster_quality(i) = {'good'};
            else
                TT.Cluster_quality(i) = {'mua'};
            end
        end
        writetable(TT, fullfile(gwfparams.dataDir, 'CSV files', [prename '_quality_metrics.csv']));
        
        tEnd = toc(t0);
        fprintf('Done in: %d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));
    end
end