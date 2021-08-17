function caraslab_mat2datChunked(Savedir)
%
% This function converts -mat files to 16 bit integer -dat files, the
% required input format for kilosort. 
% This function also runs a quick RMS-based bad channel detector.
%Written by ML Caras Mar 27 2019
% Patched by M Macedo-Lima August, 2021

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
    try
%         temp = load(mat_file);
        temp = matfile(mat_file);
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
   
    
    % Process file chunk-by-chunk using ops.NT as chunksize
    sampsToRead = size(temp, 'rawsig', 2); % total number of samples to read
    NT = ops.NT; 
    Nbatch = ceil(sampsToRead / NT); % number of data batches
    [~, ~, ~, ~, NchanTOTdefault] = loadChanMap(ops.chanMap); % function to load channel map file
    NchanTOT = getOr(ops, 'NchanTOT', NchanTOTdefault); % if NchanTOT was left empty, then overwrite with the default

    
    % Store bad channels in each batch; if a channel is good in less than
    % half of the batches, exclude it
    goodchannels_by_batch = zeros(NchanTOT, Nbatch);
    fprintf('Writing raw binary file: %s.......\n', ops.fbinary)
    fprintf('Processing %d chunks.......\n', Nbatch)
    fid = fopen(ops.fbinary,'w');
    if fid == -1
        fprintf('Cannot create binary file!')
        return
    end
    t0 = tic;
    lineLength_toUpdate = 0;
    for ibatch = 1:Nbatch
        offset = NT * (ibatch-1);
        
        if offset+NT > sampsToRead
            buff = temp.rawsig(:,offset+1:end);
        else
            buff = temp.rawsig(:,offset+1:offset+NT);
        end
        

        % RMS detection of bad channels
        goodchannels_by_batch(:, ibatch) = caraslab_rms_badChannels(buff, 0);

        % The TDT outputs a very small amplitude. After visual inspection,
        % multiplying it by this 30k gain factor works very well.
        fwrite(fid, buff(:)*30000, 'int16'); 
        
        % Update last line in console as a poor-man's progress bar
        fprintf(repmat('\b',1,lineLength_toUpdate));
        lineLength_toUpdate = fprintf('Completed %d out %d of chunks.......\n', ibatch, Nbatch);
    end
    fclose(fid);
    tEnd = toc(t0);
    
    % Count the number of batches in which each channel is good nad
    % eliminate the ones that are good in less than half the batches
    goodchannels_by_batch = sum(goodchannels_by_batch, 2)/Nbatch;
    ops.igood = goodchannels_by_batch > 0.5;
    
    [badchannels, ~] = find(ops.igood == 0);
    fprintf('RMS noise of channel(s) %d is too high or too low! Removed from analysis!\n',badchannels)
    
    
    fprintf('Finished in: %d minutes and %f seconds\n', floor(tEnd/60),rem(tEnd,60));
    %Save configuration file
    save(fullfile(cur_savedir, 'config.mat'),'ops')
end