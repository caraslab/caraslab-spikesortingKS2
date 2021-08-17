function tstart = fetch_tstart_from_behavior(csv_dir)
% This function grabs the first timestamp from a behavioral file existing inside
% the csv_dir; these files are now hardcoded to end in *spoutTimestamps.csv or  
% *trialInfo.csv
% This is done so noise can be excluded from the beginning of each file in a
% more automated manner
% If Aversive is in the name of the folder, the first spout onset is
% selected; if Aversive is not in the name of the folder, the first
% stimulus onset is selected.

try
    if contains(csv_dir, 'Aversive')
        behav_file = dir(fullfile(csv_dir, '*spoutTimestamps.csv'));

        behav_table = readtable(fullfile(behav_file.folder, behav_file.name));

        tstart = behav_table.Spout_onset(1);
    else
        behav_file = dir(fullfile(csv_dir, '*trialInfo.csv'));
        
        % Do something funy here to get around having behavioral files
        % inside the concat folders. Throw an exception if behav_file has
        % more than one item; It is a 'silent' exception caught below
        if length(behav_file) > 1
            ME = MException('MATLAB:narginchk:notEnoughInputs', ...
                'Processing a concatenated file folder...');
            throw(ME)
        end
        behav_table = readtable(fullfile(behav_file.folder, behav_file.name));

        tstart = behav_table.Trial_onset(1);

        % Subtract 5 seconds from tstart
        tstart = tstart - 5;
    end
catch ME
    % If there is no behavior file just output 0.
    if strcmp(ME.identifier, 'MATLAB:narginchk:notEnoughInputs')
        fprintf('No behavior file found for: %s\n Set trange manually before sorting\n', csv_dir)
        tstart = 0;
    else
        throw(ME)
    end
end