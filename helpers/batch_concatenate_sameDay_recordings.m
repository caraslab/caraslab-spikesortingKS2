function caraslab_concatenate_sameDay_recordings(Savedir, sel)

% Only select folders with data in them
datafolders_names = uigetfile_n_dir(Savedir,'Select data directory');

datafolders = {};
for i=1:length(datafolders_names)
    [~, datafolders{end+1}, ~] = fileparts(datafolders_names{i});
end

% TODO: Only grab folders with .dat files in them

% Sort according to dates and times in folder names
date_time = regexp(datafolders, '\d+', 'match');
recording_dates = [];
recording_times = [];
for i=1:length(date_time)
   recording_dates = [recording_dates str2num(date_time{i}{1})];
   recording_times = [recording_times str2num(date_time{i}{2})];
end

% now sort hierarchically first date, then time
temp_cell = horzcat(datafolders', num2cell([recording_dates' recording_times']) );
temp_cell = sortrows(temp_cell, [2 3]);
datafolders = temp_cell(:,1)';

unique_days = unique(recording_dates);
for day_idx=1:length(unique_days)
    cur_day_datafolders = strfind(datafolders, num2str(unique_days(day_idx)));
    cur_day_datafolders = datafolders(~cellfun('isempty', cur_day_datafolders));
    
    %Skip if cur_day only has 1 recording
    if length(cur_day_datafolders) == 1
        continue
    end
    
    output_dir = [num2str(unique_days(day_idx)) '_concat'];
    
    output_file_name = output_dir;
    
    % Print names of files to check order
    fprintf('\nConcatenating files in the following order:\n')
    for i = 1:length(cur_day_datafolders)
        fprintf('%s\n',cur_day_datafolders{i}) % print each file to make sure it's in order
        % need the sorting performed above!!!
    end

    full_output_dir = fullfile(Savedir, output_dir);
    mkdir(full_output_dir);

%     fidC        = fopen(fullfile(full_output_dir, [output_file_name '.dat']),  'w'); % Write concatenated recording
    fidC        = fopen(fullfile(full_output_dir, [output_file_name '_CLEAN.dat']),  'w'); % Write concatenated recording
    session_names = {};
    break_points = [];
    tic;
    for i = 1:numel(cur_day_datafolders)
        cur_path.name = cur_day_datafolders{i};
        cur_savedir = [Savedir filesep cur_path.name];

        %Load in configuration file (contains ops struct)
        % Catch error if -mat file is not found
        try
            load(fullfile(cur_savedir, 'config.mat'));
        catch ME
            if strcmp(ME.identifier, 'MATLAB:load:couldNotReadFile')
                fprintf('\n-mat file not found\n')
                continue
            else
                fprintf(ME.identifier)
                fprintf(ME.message)
                continue
            end
        end

        NchanTOT = ops.NchanTOT;
        NT = ops.NT;

        %Start timer

%         fprintf('\nReading raw file: %s\n', ops.fbinary)
%         fid         = fopen(ops.fbinary, 'r'); % open current raw data
        fprintf('\nReading raw file: %s\n', ops.fclean)
        fid         = fopen(ops.fclean, 'r'); % open current raw data

        session_names{end+1} = dir(ops.fbinary).name;

        if i > 1
            break_points = [break_points; get_file_size(ops.fbinary)/NchanTOT/2 + break_points(i-1)];
        else
            break_points = [break_points; get_file_size(ops.fbinary)/NchanTOT/2];
        end

        while ~feof(fid)  % read until end of file
            buff = fread(fid, [NchanTOT NT], '*int16'); % read and reshape. Assumes int16 data (which should perhaps change to an option)
            fwrite(fidC, buff, 'int16'); % write this batch to concatenated file
        end
        fclose(fid); % close the files
    end
    % Output csv breakpoints
    ret_table = cell2table(session_names', 'VariableNames', {'Session_file'});
    ret_table.Break_point = break_points;
    writetable(ret_table, fullfile(full_output_dir, [output_file_name '_breakpoints.csv']));

    fclose(fidC);

    tEnd = toc;
    fprintf('\nDone in: %d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));
end

