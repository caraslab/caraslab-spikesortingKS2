function caraslab_concatenate_sameDepth_recordings(Savedir, sel, NchanTOT, NT)
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

% Written by M Macedo-Lima October, 2020
if ~sel
    datafolders = caraslab_lsdir(Savedir);
    datafolders = {datafolders.name};

elseif sel  
    %Prompt user to select folder
    datafolders_names = uigetfile_n_dir(Savedir,'Select data directory');
    datafolders = {};
    for i=1:length(datafolders_names)
        [~, datafolders{end+1}, ~] = fileparts(datafolders_names{i});
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

output_file_name = [num2str(recording_dates(1)) '-' num2str(recording_dates(end)) '_concat'];
full_output_dir = fullfile(Savedir, output_file_name);
mkdir(full_output_dir);

fidC        = fopen(fullfile(full_output_dir, [output_file_name '_CLEAN.dat']),  'w'); % Write concatenated recording
session_names = {};
break_points = [];
tic;
for i = 1:numel(datafolders)
    cur_path.name = datafolders{i};
    cur_dir = [Savedir filesep cur_path.name];
    
    %Start timer
    
    dir_fbinary = dir(fullfile(cur_dir, '*CLEAN.dat'));  % presume file is cleaned
    fbinary = fullfile(cur_dir, dir_fbinary(1).name);
    
    fprintf('\nReading raw file: %s\n', fbinary)
    fid         = fopen(fbinary, 'r'); % open current raw data

    
    % Check if breakpoints file exists in folder; this indicates this
    % current day consists of previously concatenated recordings. If this
    % is the case, just read that file
    cur_breakpoint_file = [dir(fullfile(cur_dir, '*breakpoints.csv'))];
    if ~isempty(cur_breakpoint_file)
        cur_breakpoint_table = readtable(fullfile(cur_dir, cur_breakpoint_file(1).name), 'Delimiter', ',');
        session_names = [session_names cur_breakpoint_table.Session_file'];
        
        % concatenate break_points adding the last value
        if ~isempty(break_points)
            break_points = [break_points; break_points(end) + cur_breakpoint_table.Break_point];  
        else
            break_points = [break_points; cur_breakpoint_table.Break_point]; 
        end
    else
        session_names{end+1} = dir(fbinary).name;
        if i > 1
            
            break_points = [break_points; get_file_size(fbinary)/NchanTOT/2 + break_points(i-1)];
        else
            break_points = [break_points; get_file_size(fbinary)/NchanTOT/2];
        end
    end
    
    while ~feof(fid)  % read until end of file
        buff = fread(fid, [NchanTOT NT], '*int16'); % read and reshape. Assumes int16 data (which should perhaps change to an option)
        fwrite(fidC, buff, 'int16'); % write this batch to concatenated file
    end
    fclose(fid); % close the files
end

% Output csv breakpoints
ret_table = cell2table(session_names(:), 'VariableNames', {'Session_file'});
ret_table.Break_point = break_points;
writetable(ret_table, fullfile(full_output_dir, [output_file_name '_breakpoints.csv']));

fclose(fidC);

tEnd = toc;
fprintf('\nDone in: %d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));
