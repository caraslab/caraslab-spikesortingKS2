function caraslab_reformat_OpenEphysGUI_data(input_dir, output_dir, format, data_channel_idx, adc_channel_idx)
    %epData = caras_lab_reformat_synapse_data(Tankdir,Savedir,sel);
    %   Function to reformat and save ephys data from OpenEphys GUI.
    %
    %   Input variables:
    %       Tankdir:    path to raw recording directory
    %
    %       Savedir:    path to directory where -mat and -csv files will be saved
    %
    %       format:     'oe':      OpenEphys format
    %                   'binary':   Binary format
    %  
    %   Output:
    %       Reformatted .dat file with ephys data
    %       .info file containing ssupporting information (can be read as a -mat)
    %       CSV files with TTL information
    %   Written by M Macedo-Lima July, 2021
    
    % Provide defaults if user is confident these channels can be read from
    % openephys files
    if nargin < 4
        data_channel_idx = [];
        adc_channel_idx = [];
    end
    
    if nargin < 5
        adc_channel_idx = [];
    end
    
    %Check that tank directory exists and abort if it doesn't
    if ~exist(input_dir,'dir')
        fprintf('\n Tank directory does not exist!!\n')
        return
    end


    %Check if save directory exists. If it doesn't, create it now.
    if ~exist(output_dir,'dir')
        [success,message,messageID] = mkdir(output_dir);

        %Stop if directory cannot be created, and display reason why
        if ~success
            message %#ok<*NOPRT>
            messageID
            return
        end   
    end


    %Prompt user to select folder
    datafolders_names = uigetfile_n_dir(input_dir,'Select data directory');
    blocknames = {};
    for i=1:length(datafolders_names)
        [~, blocknames{end+1}, ~] = fileparts(datafolders_names{i});
    end


    %Check that at least one block has been selected
    if isempty(blocknames)
        fprintf('\n No BLOCKS could be found!!\n')
        return
    end


    %For each block
    for i = 1:numel(blocknames)
        t0 = tic;
        cur_path.name = blocknames{i};
        cur_savedir = [output_dir filesep cur_path.name];
        data_filename = fullfile(cur_savedir, [cur_path.name '.dat']);
        adc_filename = fullfile(cur_savedir, [cur_path.name '_ADC.dat']);

        events_filename = fullfile(cur_savedir, [cur_path.name '.info']);

        cur_subdir = dir(fullfile(input_dir, cur_path.name, 'Record Node*'));  % assume 1 recording node
        fullpath = fullfile(cur_subdir.folder, cur_subdir.name);

        %Convert tank data to -mat and display elapsed time
        fprintf('\n======================================================\n')
        fprintf('Processing ephys data, %s.......\n', cur_path.name)

        % Try to read file. If missing, skip to the next folder
        mkdir(cur_savedir);
        % In case this folder is still absent
        mkdir(fullfile(cur_savedir, 'CSV files'));

        t0 = tic;
        chunk_size = 1800000; % Lower this number if running out of memory
        fid_data = fopen(data_filename,'w');
        fid_adc = fopen(adc_filename,'w');

        try
            switch format
                case 'oe'
                    handle_oe_data();
                case 'binary'
                    handle_binary_data();
            end
        catch ME
            if strcmp(ME.identifier, 'MATLAB:load:couldNotReadFile')
                fprintf('\nFile not found\n')
                continue
            else
                fprintf(2, ME.identifier)
                fprintf(2, ME.message)
                break
            end
        end

        tEnd = toc(t0);
        fprintf('\n~~~~~~\nFinished in: %d minutes and %f seconds\n~~~~~~\n', floor(tEnd/60),rem(tEnd,60));
    end
    fprintf('\n\n ##### Finished reformatting and saving data files.\n\n')

    
    %% Helpers to handle appropriate formats
    function handle_oe_data()
        % Data channels have the naming convention *CHXX.continuous
        % Read recording info to get channel order; ADC channels will also
        % be in the mix; Channels should come out in order
        session_info = get_session_info(fullpath); 
        
        % Find the index of the recording node (assume only one)
        node_idx = find(contains(session_info.processors(:,2),'Filters/Record Node'));
        
        % Something weird happens sometimes
        all_channels = session_info.processors{node_idx, 3}{1};
        if isempty(all_channels)
            all_channels = session_info.processors{node_idx, 3}{2};
        end
        
        data_channels = all_channels(contains(all_channels,'CH'));
        adc_channels = all_channels(contains(all_channels,'ADC'));
        
        % Weird bug in OpenEphys GUI sometimes names these differently
        % Tweak this mannually
        if isempty(data_channels)
            disp('No channels named CH; Make sure you provided channel numbers')
            if isempty(data_channel_idx)
                disp('Please provide array with data channel numbers, e.g. 1:64')
                return
            end
            
            data_channels = all_channels(data_channel_idx);
            adc_channels = all_channels(adc_channel_idx);
        end

        
        %% Read data channels and add to .dat
        fprintf('Concatenating data channels in chunks:\n')
        eof = 0;
        chunk_counter = 1;
        t1 = 0;
        while ~eof
            disp(['Chunk counter: ' num2str(chunk_counter) '...']);
            rawsig = [];
            t2 = t1 + chunk_size;  % Update t2
            disp(['Reading data channels...']);
            for ch_idx=1:length(data_channels)
                [cur_ch_data, ~, ~, is_eof] = load_open_ephys_data_chunked(fullfile(fullpath, data_channels{ch_idx}), t1, t2, 'samples');
                rawsig = [rawsig; cur_ch_data'];
                if is_eof
                    eof = 1;
                end

            end
            disp(['Writing to file...']);
            fwrite(fid_data, rawsig, 'int16');
            t1 = t2;  % Update t1
            chunk_counter = chunk_counter + 1;
        end
        fclose(fid_data);
        
        %% Read ADC channels and add to .dat
        if ~isempty(adc_channels)
            fprintf('Concatenating ADC channels in chunks in the following order:\n')
            eof = 0;
            chunk_counter = 1;
            t1 = 0;
            while ~eof
                disp(['Chunk number: ' num2str(chunk_counter) '...']);
                rawsig = [];
                t2 = t1 + chunk_size;  % Update t2
                disp(['Reading ADC channels...']);

                for ch_idx=1:length(adc_channels)
                    [cur_ch_data, ~, ~, is_eof] = load_open_ephys_data_chunked(fullfile(fullpath, adc_channels{ch_idx}), t1, t2, 'samples');
                    rawsig = [rawsig; cur_ch_data'];
                    if is_eof
                        eof = 1;
                    end

                end
                disp(['Writing to file...']);
                fwrite(fid_adc, rawsig, 'int16');
                t1 = t2;  % Update t1
                chunk_counter = chunk_counter + 1;
            end
            fclose(fid_adc);
        end
        
        clear rawsig cur_ch_data;
        
        %% Read DAC channels; Bundle in a .info file and output as CSV
        % unpack them in the behavior pipeline; this struct is meant to
        % look somewhat similar to the TDT epData that comes out of Synapse
        
        % For Rig2 with ePsych:
        % 0: DAC1 = sound on/off
        % 1: DAC2 = spout on/off
        % 2: DAC3 = trial start/end
        fprintf('Reading Events channels:\n')
        
        % Load only a little bit of a channel file to get the zero timestamp info
        [~, data_timestamps, ~, ~] = load_open_ephys_data_chunked(fullfile(fullpath, data_channels{1}), 0, 5, 'samples');

        [event_ids, timestamps, info] = load_open_ephys_data_faster(fullfile(fullpath, 'all_channels.events'));
        epData.event_ids = event_ids + 1;  % Convert to 1-base index
        epData.event_states = info.eventId;
        epData.timestamps = timestamps - data_timestamps(1); % Zero TTL timestamps based on the first sampled data  time
        epData.info.blockname = cur_path.name;
                
        % Grab date and timestamp from info
        block_date_timestamp = info.header.date_created;
        block_date_timestamp = datevec(block_date_timestamp, 'dd-mmm-yyyy HH:MM:SS');
        epData.info.StartTime = block_date_timestamp;  % TDT-like
        
        save(events_filename, 'epData','-v7.3');
        
        % Output each channel with events as separate csv with onset,
        % offset and duration
        unique_dacs = unique(epData.event_ids);
        for cur_event_id_idx=1:length(unique_dacs)
            cur_event_id = unique_dacs(cur_event_id_idx);
            cur_event_mask = epData.event_ids == cur_event_id;
            cur_event_states = epData.event_states(cur_event_mask);            
            cur_timestamps = epData.timestamps(cur_event_mask);
            
            cur_onsets = cur_timestamps(cur_event_states == 1);
            cur_offsets = cur_timestamps(cur_event_states == 0);
            
            % Handle DAC exceptions here
            % Skip DAC if either onset or offset are completely absent
            if isempty(cur_onsets) || isempty(cur_offsets)
                continue
            end
            
            % Remove first offset if lower than first onset 
            if cur_offsets(1) < cur_onsets(1)
                cur_offsets = cur_offsets(2:end);
            end
            
            % Remove last onset if length mismatch
            if length(cur_onsets) ~= length(cur_offsets)
                cur_onsets = cur_onsets(1:end-1);
            end
            
            % Calulate durations
            cur_durations = cur_offsets - cur_onsets;
            
            % Convert to table and output csv
            
            fileID = fopen(fullfile(cur_savedir, 'CSV files', ...
                [cur_path.name '_DAC' int2str(cur_event_id) '.csv']), 'w');

            header = {'Onset', 'Offset', 'Duration'};
            fprintf(fileID,'%s,%s,%s\n', header{:});
            nrows = length(cur_onsets);
            for idx = 1:nrows
                output_cell = {cur_onsets(idx), cur_offsets(idx), cur_durations(idx)};

                fprintf(fileID,'%f,%f,%f\n', output_cell{:});
            end
            fclose(fileID);

        end
    end


    function handle_binary_data()
        oebin_filedir = dir(fullfile(fullpath, '**', '*.oebin'));
        
        all_data = load_open_ephys_binary(fullfile(oebin_filedir.folder, oebin_filedir.name), 'continuous', 1, 'mmap');
        
        data_channel_indeces = contains({all_data.Header.channels.channel_name}, 'CH');
        adc_channel_index = contains({all_data.Header.channels.channel_name}, 'ADC');

        % Weird bug in OpenEphys GUI sometimes names these differently
        % Tweak this mannually
        if sum(data_channel_indeces) == 0
            disp(['No channels named CH; tweak mannually (Lines 277,278)'])
            allchan = [data_channel_idx adc_channel_idx];
            data_channel_indeces = ismember(allchan, data_channel_idx);
            adc_channel_index = ismember(allchan, adc_channel_idx);
        end

        
        %% Read data channels and add to .dat
        fprintf('Concatenating data channels in chunks:\n')
        eof = 0;
        chunk_counter = 1;
        t1 = 1;
        while ~eof
            disp(['Chunk counter: ' num2str(chunk_counter) '...']);
            t2 = t1 + chunk_size - 1;  % Update t2
            disp(['Reading data channels...']);
            
            if t2 < length(all_data.Timestamps)
            
                rawsig = all_data.Data.Data(1).mapped(data_channel_indeces, t1:t2);
            else
                rawsig = all_data.Data.Data(1).mapped(data_channel_indeces, t1:end);
                eof = 1;
            end
            disp(['Writing to file...']);
            fwrite(fid_data, rawsig, 'int16');
            t1 = t2 + 1;  % Update t1
            chunk_counter = chunk_counter + 1;
        end
        fclose(fid_data);
        rawsig = rawsig';  % Transpose for [n_ch n_samp] size;
        
        %% Read ADC channels and add to .dat
        if ~isempty(adc_channel_index)
            fprintf('Concatenating ADC channels in chunks in the following order:\n')
            eof = 0;
            chunk_counter = 1;
            t1 = 1;
            while ~eof
                disp(['Chunk number: ' num2str(chunk_counter) '...']);
                t2 = t1 + chunk_size - 1;  % Update t2
                disp(['Reading ADC channels...']);
                if t2 < length(all_data.Timestamps)
                    rawsig = all_data.Data.Data(1).mapped(adc_channel_index, t1:t2);
                else
                    rawsig = all_data.Data.Data(1).mapped(adc_channel_index, t1:end);
                    eof = 1;
                end
                
                disp(['Writing to file...']);
                fwrite(fid_adc, rawsig, 'int16');
                t1 = t2 + 1;  % Update t1
                chunk_counter = chunk_counter + 1;
            end
            fclose(fid_adc);
        end
        
        clear rawsig cur_ch_data;
        
        %% Read DAC channels; Bundle in a .info file and output as CSV
        % unpack them in the behavior pipeline; this struct is meant to
        % look somewhat similar to the TDT epData that comes out of Synapse
        
        % For Rig2 with ePsych:
        % 0: DAC1 = sound on/off
        % 1: DAC2 = spout on/off
        % 2: DAC3 = trial start/end
        fprintf('Reading Events channels:\n')
        dac_data = load_open_ephys_binary(fullfile(oebin_filedir.folder, oebin_filedir.name), 'events', 1);

        epData.event_states = uint16(dac_data.Data) ./ dac_data.ChannelIndex;
        epData.event_ids = dac_data.ChannelIndex;  % convert to 1-base index
        epData.timestamps = double(dac_data.Timestamps - all_data.Timestamps(1)) / dac_data.Header.sample_rate; % Zero TTL timestamps based on the first sampled data  time
        epData.info.blockname = cur_path.name;
               
        % Grab date and timestamp from info
        block_date_timestamp = info.header.date_created;
        block_date_timestamp = datevec(block_date_timestamp, 'dd-mmm-yyyy HH:MM:SS');
        epData.info.StartTime = block_date_timestamp;  % TDT-like
        
        save(events_filename, 'epData','-v7.3');
        
        % Output each channel with events as separate csv with onset,
        % offset and duration
        unique_dacs = unique(epData.event_ids);
        for cur_event_id_idx=1:length(unique_dacs)
            cur_event_id = unique_dacs(cur_event_id_idx);
            cur_event_mask = epData.event_ids == cur_event_id;
            cur_event_states = epData.event_states(cur_event_mask);            
            cur_timestamps = epData.timestamps(cur_event_mask);
            
            cur_onsets = cur_timestamps(cur_event_states == 1);
            cur_offsets = cur_timestamps(cur_event_states == 0);
            
            % Handle DAC exceptions here
            % Skip DAC if either onset or offset are completely absent
            if isempty(cur_onsets) || isempty(cur_offsets)
                continue
            end
            
            % Remove first offset if lower than first onset 
            if cur_offsets(1) < cur_onsets(1)
                cur_offsets = cur_offsets(2:end);
            end
            % Remove last onset if length mismatch
            if length(cur_onsets) ~= length(cur_offsets)
                cur_onsets = cur_onsets(1:end-1);
            end
            
            % Calulate durations
            cur_durations = cur_offsets - cur_onsets;
            
            % Convert to table and output csv
            
            fileID = fopen(fullfile(cur_savedir, 'CSV files', ...
                [cur_path.name '_DAC' int2str(cur_event_id) '.csv']), 'w');

            header = {'Onset', 'Offset', 'Duration'};
            fprintf(fileID,'%s,%s,%s\n', header{:});
            nrows = length(cur_onsets);
            for idx = 1:nrows
                output_cell = {cur_onsets(idx), cur_offsets(idx), cur_durations(idx)};

                fprintf(fileID,'%f,%f,%f\n', output_cell{:});
            end
            fclose(fileID);
        end
    end
end
