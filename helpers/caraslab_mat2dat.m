function caraslab_mat2dat(Savedir)
%
% This function converts -mat files to 16 bit integer -dat files, the
% required input format for kilosort. 
% This function also runs a quick RMS-based bad channel detector. The
% output is not used for anything but it tells the user of unknown
% potential bad channels
%
%Written by ML Caras Mar 27 2019
% Patched by M Macedo-Lima October, 2020


%Prompt user to select folder
datafolders_names = uigetfile_n_dir(Savedir,'Select data directory');
datafolders = {};
for i=1:length(datafolders_names)
    [~, datafolders{end+1}, ~] = fileparts(datafolders_names{i});
end


%For each data folder...
for i = 1:numel(datafolders)
    clear ops temp dat matfile temp_raw
    
    cur_path.name = datafolders{i};
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
    
    %Find the -mat file to convert
    mat_file = ops.rawdata;
    
    %Load -mat data file
    fprintf('Loading -mat file: %s.......\n', mat_file)
    tic;
    try
        temp = load(mat_file);
%         temp = matfile(mat_file);
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
    
    tEnd = toc;
    fprintf('Loaded in: %d minutes and %f seconds\n', floor(tEnd/60),rem(tEnd,60));
    
    % RMS detection of bad channels
    % Has to be done before amplitude range normalization
    ops.igood = caraslab_rms_badChannels(temp.rawsig);
    
    fprintf('Writing raw binary file: %s.......\n', ops.fbinary)
    t0 = tic;
    
    fid = fopen(ops.fbinary,'w');
    if fid == -1
        fprintf('Cannot create binary file!')
        return
    end
    
    % The TDT outputs a very small amplitude. After visual inspection,
    % multiplying it by this 30k gain factor works very well.
    fwrite(fid, temp.rawsig(:)*30000, 'int16'); 
    fclose(fid);
    tEnd = toc(t0);
    fprintf('Finished in: %d minutes and %f seconds\n', floor(tEnd/60),rem(tEnd,60));
    
    %Save configuration file
    save(fullfile(cur_savedir, 'config.mat'),'ops')
end