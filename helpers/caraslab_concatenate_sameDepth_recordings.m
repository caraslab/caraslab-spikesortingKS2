function caraslab_concatenate_sameDepth_recordings(Savedir, recording_type)
% This function concatenates recordings in provided folders into a single
% .dat file. It orders them by date, then concatenates.
% It also works when concatenating previsouly concatenated recordings. In
% this case, it will search for a breakpoints.csv file and append
% breakpoints to that
%
%   sel:        if 0 or omitted, program will cycle through all files
%               in the data directory. 
%                   
%               if 1, user will be prompted to select a file
%
% The outputs are a concatenated .dat file and a .csv file listing the
% breakpoints (in samples) at which the original recordings ended. This file is
% important to realign spout and stimulus timestamps which start at 0

%Prompt user to select folder
datafolders_names = uigetfile_n_dir(Savedir,'Select data directory');
datafolders = {};
for i=1:length(datafolders_names)
    [~, datafolders{end+1}, ~] = fileparts(datafolders_names{i});
end

% Remove SUBJ ID from folder names
subj_in_filename = 0;
for i = 1:length(datafolders)
    if contains(datafolders{i}, 'SUBJ-ID')
        subj_in_filename = 1;
        df = split(datafolders{i}, '_');
        id = df{1};
        datafolders{i} = append(df{2},'_',df{3},'_',df{4});
    end
end


% Sort according to dates then times in folder names
date_time = regexp(datafolders, '\d+', 'match');
recording_dates = [];
recording_times = [];
for i=1:length(date_time)
   recording_dates = [recording_dates str2num(date_time{i}{1})];
   
   if length(date_time{i}) > 1
      recording_times = [recording_times str2num(date_time{i}{2})];
   else
       recording_times = [recording_times 0];
   end
end

% now sort hierarchically first date, then time
temp_cell = horzcat(datafolders', num2cell([recording_dates' recording_times']) );
temp_cell = sortrows(temp_cell, [2 3]);
datafolders = temp_cell(:,1)';

% Print names of files to check order
fprintf('\nConcatenating files in the following order:\n')
for i = 1:length(datafolders)
    fprintf('%s\n',datafolders{i}) % print each file to make sure it's in order
    % need the sorting performed above!!!
end

output_folder_name = [num2str(recording_dates(1)) '-' num2str(recording_dates(end)) '_concat'];
full_output_folder_name = fullfile(Savedir, output_folder_name);
mkdir(full_output_folder_name);

fidC        = fopen(fullfile(full_output_folder_name, [output_folder_name '_CLEAN.dat']),  'w'); % Write concatenated recording
session_names = {};
break_points = [];
breakpoints_seconds = [];
tranges = [];
cumulative_tranges =[];
badchannels = [];
tic;
for i = 1:numel(datafolders)
    cur_path.name = datafolders{i};
    cur_sourcedir = fullfile(Savedir,cur_path.name);
    
    %Start timer
    
    dir_fbinary = dir(fullfile(cur_sourcedir, '*CLEAN.dat'));  % presume file is cleaned
    fbinary = fullfile(cur_sourcedir, dir_fbinary(1).name);
    
    fprintf('\nReading raw file: %s\n', fbinary)
    fid         = fopen(fbinary, 'r'); % open current raw data

    %Load in configuration file (contains ops struct)
    % Catch error if -mat file is not found
    try
        load(fullfile(cur_sourcedir, 'config.mat'));
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
    badchannels = [badchannels ops.badchannels];

    % Check if breakpoints file exists in folder; this indicates this
    % current day consists of previously concatenated recordings. If this
    % is the case, just read that file
    cur_breakpoint_file = [dir(fullfile(cur_sourcedir, '*breakpoints.csv'))];
    if ~isempty(cur_breakpoint_file)
        cur_breakpoint_table = readtable(fullfile(cur_sourcedir, cur_breakpoint_file(1).name), 'Delimiter', ',');
        session_names = [session_names cur_breakpoint_table.Session_file'];
        cur_breakpoint = cur_breakpoint_table.Break_point;
        % concatenate break_points adding the last value
        if ~isempty(break_points)
            break_points = [break_points; break_points(end) + cur_breakpoint];  
            breakpoints_seconds = [breakpoints_seconds; breakpoints_seconds(end) + cur_break_point/ops.fs];
        else
            break_points = [break_points; cur_breakpoint]; 
            breakpoints_seconds = [breakpoints_seconds; cur_break_point/ops.fs];
        end
    else
        session_names{end+1} = dir(fbinary).name;
        cur_breakpoint = get_file_size(fbinary)/NchanTOT/2;
        if i > 1
            break_points = [break_points; cur_breakpoint + break_points(i-1)];
            breakpoints_seconds = [breakpoints_seconds; cur_breakpoint/ops.fs + breakpoints_seconds(i-1)];
        else
            break_points = [break_points; cur_breakpoint];
            breakpoints_seconds = [breakpoints_seconds; cur_breakpoint/ops.fs];
        end
    end

    % create offset variables relevant to the concatenation process
    if i > 1
        previous_breakpoint = break_points(i-1);
    else
        previous_breakpoint = 0;
    end

    if isfield('ops', 'concat_cumulative_tranges')
        cumulative_tranges = [cumulative_tranges; ops.concat_cumulative_tranges];
    else 
        cumulative_tranges = [cumulative_tranges; breakpoints_seconds(end)];
    end
    %% Handle the recording file
    while ~feof(fid)  % read until end of file
        buff = fread(fid, [NchanTOT NT], '*int16'); % read and reshape. Assumes int16 data (which should perhaps change to an option)
        fwrite(fidC, buff, 'int16'); % write this batch to concatenated file
    end
    fclose(fid); % close the files
    
    %% Finish up by copying over the ePsych files
    % Make copies of the behavioral timestamps and metadata (if they exist) into the
    % CSV files folder
    % In case this folder is still absent
    mkdir(fullfile(full_output_folder_name, 'CSV files'));

    % find spout file
    spout_filename = fullfile(cur_sourcedir, 'CSV files', '*spoutTimestamps.csv');
    if ~isempty(spout_filename)
        copy_data_files_from_dir(fullfile(full_output_folder_name, 'CSV files'), spout_filename)
%         copyfile(fullfile(spout_filename.folder, spout_filename.name), fullfile(full_output_folder_name, 'CSV files', spout_filename.name));
    end

     % find trial info file
    trialInfo_filename = fullfile(cur_sourcedir, 'CSV files', '*trialInfo.csv');
    if ~isempty(trialInfo_filename)
        copy_data_files_from_dir(fullfile(full_output_folder_name, 'CSV files'), trialInfo_filename)
%         copyfile(fullfile(trialInfo_filename.folder, trialInfo_filename.name), fullfile(full_output_folder_name, 'CSV files', trialInfo_filename.name));
    end

     % find metadata file
    metadata_filename = fullfile(cur_sourcedir, 'CSV files', '*ePsychMetadata.mat');
    if ~isempty(metadata_filename)
        copy_data_files_from_dir(fullfile(full_output_folder_name, 'CSV files'), metadata_filename)
%         copyfile(fullfile(metadata_filename.folder, metadata_filename.name), fullfile(full_output_folder_name, 'CSV files', metadata_filename.name));
    end
end

%% Output csv breakpoints
% Grab subject name; this is specific for my naming convention; Should be tweaked if yours is different
split_dir = split(Savedir, filesep); 

subj_id = split(split_dir{end}, '-');
subj_id = join(subj_id(1:3), '-'); 
subj_id = subj_id{1}; 
ret_table = cell2table(session_names', 'VariableNames', {'Session_file'});
ret_table.Break_point = break_points;
ret_table.Break_point_seconds = breakpoints_seconds;
writetable(ret_table, fullfile(full_output_folder_name, 'CSV files', [subj_id '_' output_folder_name '_breakpoints.csv']));

%% Close file and save Config
fclose(fidC);

% Create new config.mat
caraslab_createconfig(Savedir,ops.chanMap, unique(badchannels), 0, recording_type, full_output_folder_name)
load(fullfile(full_output_folder_name, 'config.mat'));
ops.concat_tranges = tranges;
ops.concat_cumulative_tranges = cumulative_tranges;

save(fullfile(full_output_folder_name, 'config.mat'), 'ops');

tEnd = toc;
fprintf('\nDone in: %d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));
