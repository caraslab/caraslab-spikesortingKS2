function caraslab_reformat_intanOEG_data(input_dir,output_dir,sel)
%epData = caras_lab_reformat_synapse_data(Tankdir,Savedir,sel);
%   Function to reformat and save ephys data from OpenEphys GUI (OEG format).
%
%   Input variables:
%       Tankdir:    path to raw recording directory
%
%       Savedir:    path to directory where -mat and -csv files will be saved
%
%       sel:        if 0 or omitted, program will cycle through all BLOCKS
%                   in the tank directory. 
%                   
%                   if 1, user will be prompted to select a BLOCK
%  
%   Uses TDTbin2mat to reformat tank data to a matlab struct. 
%   Two files are saved:    
%       (1) A -mat file containing an MxN matrix of raw voltages, where 
%               M = the number of channels
%               N = the number of samples
%
%       (2) A -info file containing supporting information, including
%               sampling rate, epocs, and timing
%       
%
%   Written by ML Caras Mar 22, 2019 
%   patched by M Macedo-Lima 9/8/20


%Default to selecting folders
if nargin < 3
    sel = 1;   
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

if ~sel
    %Get a list of all BLOCKS in the tank directory
    blocks = caraslab_lsdir(input_dir);
    blocknames = {blocks.name};
    

elseif sel  
    %Prompt user to select folder
    datafolders_names = uigetfile_n_dir(input_dir,'Select data directory');
    blocknames = {};
    for i=1:length(datafolders_names)
        [~, blocknames{end+1}, ~] = fileparts(datafolders_names{i});
    end
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
    
%     try
        
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
        if length(data_channels) == 0
            disp(['No channels named CH; tweak mannually (Lines 125,126) if needed'])
            data_channels = all_channels(1:64);  % Look into your input folder and check that the first Nchannel files are your data files
            adc_channels = all_channels(65);  % Check that your next channel is an adc file, or make this [] if you want to ignore it
        end

        
        %% Read data channels and add to .dat
        fprintf('Concatenating data channels in chunks:\n')
        %eof = {0};
        eof = 0;
        chunk_counter = 1;
        t1 = 0;
%         while ~eof
%             disp(['Chunk counter: ' num2str(chunk_counter) '...']);
%             rawsig = [];
%             t2 = t1 + chunk_size;  % Update t2
%             disp(['Reading data channels...']);
%             for ch_idx=1:length(data_channels)
%                 [cur_ch_data, ~, ~, is_eof] = load_open_ephys_data_chunked(fullfile(fullpath, data_channels{ch_idx}), t1, t2, 'samples');
%                 rawsig = [rawsig; cur_ch_data'];
%                 if is_eof
%                     eof = 1;
%                 end
%                 
%             end
%             disp(['Writing to file...']);
%             fwrite(fid_data, rawsig, 'int16');
%             t1 = t2;  % Update t1
%             chunk_counter = chunk_counter + 1;
%         end
        while ~all([eof{:}])
            disp(['Chunk counter: ' num2str(chunk_counter) '...']);
            rawsig = cell(length(data_channels),1);
            t2 = t1 + chunk_size;  % Update t2
            disp(['Reading data channels...']);
            parfor ch_idx=1:length(data_channels)
                [cur_ch_data, ~, ~, eof{ch_idx}] = load_open_ephys_data_chunked(fullfile(fullpath, data_channels{ch_idx}), t1, t2, 'samples');
                rawsig{ch_idx} = int16(cur_ch_data)';
            end
            rawsig = cell2mat(rawsig);
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
        epData.event_ids = event_ids;
        epData.event_states = info.eventId;
        epData.timestamps = timestamps - data_timestamps(1); % Zero TTL timestamps based on the first sampled data  time
        epData.info.blockname = cur_path.name;
        epData.info.StartTime = datevec(info.header.date_created,'dd-mmm-yyyy HH:MM:SS');
        
        save(events_filename, 'epData','-v7.3')
        
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
                [cur_path.name '_DAC' int2str(cur_event_id+1) '.csv']), 'w');

            header = {'Onset', 'Offset', 'Duration'};
            fprintf(fileID,'%s,%s,%s\n', header{:});
            nrows = length(cur_onsets);
            for idx = 1:nrows
                output_cell = {cur_onsets(idx), cur_offsets(idx), cur_durations(idx)};

                fprintf(fileID,'%f,%f,%f\n', output_cell{:});
            end
            fclose(fileID);

        end
        
%     catch ME
%         if strcmp(ME.identifier, 'MATLAB:load:couldNotReadFile')
%             fprintf('\nFile not found\n')
%             continue
%         else
%             fprintf(ME.identifier)
%             fprintf(ME.message)
%             break
%         end
%     end

    tEnd = toc(t0);
    fprintf('\n~~~~~~\nFinished in: %d minutes and %f seconds\n~~~~~~\n', floor(tEnd/60),rem(tEnd,60));

end


fprintf('\n\n ##### Finished reformatting and saving data files.\n\n')


end






