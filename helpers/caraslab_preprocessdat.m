function caraslab_preprocessdat(Savedir, inspect_artifact_removal_only)
% This function takes .dat files and employs in this order:
% 1. Comb filter (if ops.comb==1)
% 2. Median-CAR filter (if ops.CAR==1)
% 3. Kilosort-inspired GPU-based chunkwise filter
% 4. Saves a filename_CLEAN.dat file
%
%Input variables:
%
%       Savedir: path to folder containing data directories. Each directory
%                should contain a binary (-dat) data file and
%                a kilosort configuration (config.mat) file. 
%
%       sel:    if 0 or omitted, program will cycle through all folders
%               in the data directory.    
%
%               if 1, program will prompt user to select folder

%Written by M Macedo-Lima 10/05/20

%Prompt user to select folders
datafolders_names = uigetfile_n_dir(Savedir,'Select data directory');
datafolders = {};
for i=1:length(datafolders_names)
    [~, datafolders{end+1}, ~] = fileparts(datafolders_names{i});
end


%For each data folder...
for i = 1:numel(datafolders)
    clear ops
    
    cur_path.name = datafolders{i};
    cur_savedir = [Savedir filesep cur_path.name];

    %Load in configuration file (contains ops struct)
    % Catch error if -mat file is not found
    try
        load(fullfile(cur_savedir, 'config.mat'));

    catch ME
        if strcmp(ME.identifier, 'MATLAB:load:couldNotReadFile')
            fprintf('\nConfig file not found\n')
            continue
        else
            fprintf(ME.identifier)
            fprintf(ME.message)
            continue
        end
    end
    
    
    % Load config paramaters
    [chanMap, ~, ~, ~, NchanTOTdefault] = loadChanMap(ops.chanMap); % function to load channel map file
    ops.NchanTOT = getOr(ops, 'NchanTOT', NchanTOTdefault); % if NchanTOT was left empty, then overwrite with the default
    NchanTOT = ops.NchanTOT; % total number of channels in the raw binary file, including dead, auxiliary etc
    % MML edit
    if isfield(ops, 'igood')
        igood = logical(ops.igood);
    else
        igood = true(size(chanMap));
    end

    % MML edit: further remove defective channels
    if isfield(ops, 'badchannels')
        igood(ops.badchannels) = 0;
    end

    
%     NT       = ops.NT ; % number of timepoints per batch
    NT = 3*ops.fs;  % Can be significantly higher than kilosort's here, but lower it if running out of RAM
    
    bytes       = get_file_size(ops.fbinary); % size in bytes of raw binary
    nTimepoints = floor(bytes/NchanTOT/2); % number of total timepoints
    ops.tstart  = ceil(ops.trange(1) * ops.fs); % starting timepoint for processing data segment
    ops.tend    = min(nTimepoints, ceil(ops.trange(2) * ops.fs)); % ending timepoint
    ops.sampsToRead = ops.tend-ops.tstart; % total number of samples to read
    ops.twind = ops.tstart * NchanTOT*2; % skip this many bytes at the start

    Nbatch      = ceil(ops.sampsToRead /(NT-ops.ntbuff)); % number of data batches
    ops.Nbatch = Nbatch;
    NTbuff      = NT + 4*ops.ntbuff; % we need buffers on both sides for filtering

    % set up the parameters of the filter
    if isfield(ops,'fslow')&&ops.fslow<ops.fs/2
        [b1, a1] = butter(3, [ops.fshigh/ops.fs,ops.fslow/ops.fs]*2, 'bandpass'); % butterworth filter with only 3 nodes (otherwise it's unstable for float32)
    else
        [b1, a1] = butter(3, ops.fshigh/ops.fs*2, 'high'); % the default is to only do high-pass filtering at 150Hz
    end

    if getOr(ops, 'comb', 0)  % MML edit; comb filter
        N  = 407;    % Order
        BW = 2;    % Bandwidth
        Fs = ops.fs;  % Sampling Frequency
        h = fdesign.comb('Notch', 'N,BW', N, BW, Fs);
        comb_filter = design(h, 'butter');
        comb_b1= comb_filter.Numerator;
        comb_a1= comb_filter.Denominator;
    end

    %Start timer
    t0 = tic;
    fprintf('Reading raw file and applying filters for file: %s\n', ops.fbinary)
    fidC        = fopen(ops.fclean,  'w'); % MML edit; write processed data for phy
    fid         = fopen(ops.fbinary, 'r'); % open for reading raw data
    

    for ibatch = 1:Nbatch
        clear buff dataRAW datr noiseBuff
        % we'll create a binary file of batches of NT samples, which overlap consecutively on ops.ntbuff samples
        % in addition to that, we'll read another ops.ntbuff samples from before and after, to have as buffers for filtering
        % MML edit: this part is weird in the original kilosort2 code so I
        % rewrote it in a way it makes sense to me. I described the issue
        % in the link below, and Marius (Kilosort developer) didn't 
        % necessary think it was a bug, but my confusion remains. 
        % https://github.com/MouseLand/Kilosort/issues/223
        offset = max(0, ops.twind + 2*NchanTOT*(NT*(ibatch-1) - 2*ops.ntbuff)); % number of samples to start reading at.
        if offset==0
            ioffset = 0; % The very first batch has no pre-buffer, and has to be treated separately
        else
            ioffset = 2*ops.ntbuff;
        end
        fseek(fid, offset, 'bof'); % fseek to batch start in raw file

        buff = fread(fid, [NchanTOT NTbuff], 'int16'); % read and reshape. Assumes int16 data (which should perhaps change to an option)

        if isempty(buff)
            break; % this shouldn't really happen, unless we counted data batches wrong
        end

        nsampcurr = size(buff,2); % how many time samples the current batch has

        if nsampcurr<NTbuff
            buff(:, nsampcurr+1:NTbuff) = repmat(buff(:,nsampcurr), 1, NTbuff-nsampcurr); % pad with zeros, if this is the last batch
        end

        % Finally start filtering...
        % Can't use GPU acceleration for comb filter yet...
        if getOr(ops, 'comb', 0)  % MML edit; comb filter
            buff = buff';  % MML edit: transpose sooner
            buff = filter(comb_b1, comb_a1, buff);
            dataRAW = gpuArray(buff); % move int16 data to GPU
        else
            dataRAW = gpuArray(buff); % move int16 data to GPU
            dataRAW = dataRAW';
        end

        dataRAW = single(dataRAW); % convert to float32 so GPU operations are fast
        % subtract the mean from each channel
        dataRAW = dataRAW - mean(dataRAW, 1); % subtract mean of each channel

        % CAR, common average referencing by median
        if getOr(ops, 'CAR', 1)
            % MML edit:take median of good channels only
            dataRAW = dataRAW - median(dataRAW(:, chanMap(igood)), 2); % subtract median across channels
        end

        datr = filter(b1, a1, dataRAW); % causal forward filter

        datr = flipud(datr); % reverse time
        datr = filter(b1, a1, datr); % causal forward filter again
        datr = flipud(datr); % reverse time back

        datr    = datr(ioffset + (1:NT),:); % remove timepoints used as buffers
        
        if getOr(ops, 'rm_artifacts', 1)
%             inspect_results = 1;
%             fprintf('Removing artifacts.......\n')
%             t0 = tic;
            warning off;
            [datr, ~, exit_flag] = caraslab_artifact_reject(datr, ops, inspect_artifact_removal_only);
            warning on;
            if exit_flag
                close all
                fclose(fid);
                fclose(fidC);
                return
            end
%             tEnd = toc(t0);
%             fprintf('Finished in: %d minutes and %f seconds\n', floor(tEnd/60),rem(tEnd,60));
        end
        
        datr = datr';
        
        % Convert noisy parts in the beginning of file (before ops.tstart)
        % to gaussian noise using the
        % median and std of the first batch
        % This prevents kilosort from getting biased by the higher noise
        % band in those often long parts
        if ibatch == 1
            if ops.tstart > 0
                % Gather on CPU side right away because GPU might run out
                % of memory
                
                % take median of medians of good channels
                first_median = gather(double(median(median(datr(chanMap(igood), :), 2))));
                % take median std of good channels
                first_std = gather(double(median(std(datr(chanMap(igood), :), [], 2))));
                fseek(fid, 0, 'bof'); % fseek to batch start in raw file
                noiseBuff = fread(fid, [NchanTOT ops.tstart-1], 'int16'); % read and reshape. Assumes int16 data (which should perhaps change to an option)
                
                noiseBuff =  first_std*randn(size(noiseBuff)) + first_median;
                fwrite(fidC, noiseBuff, 'int16'); % write this batch to clean file
            end
        end
        
        % Convert noisy channels to gaussian noise
        
        % take median of medians of good channels
        cur_median = gather(double(median(median(datr(chanMap(igood), :), 2))));
        % take median std of good channels
        cur_std = gather(double(median(std(datr(chanMap(igood), :), [], 2))));
        
        datr(chanMap(~igood), :) = gpuArray(cur_std*randn(sum(~igood), size(datr, 2)) + cur_median);

        datr  = gather(int16(datr)); % convert to int16, and gather on the CPU side
        
        fwrite(fidC, datr, 'int16'); % write this batch to clean file
    end
    fclose(fid); % close the files
    fclose(fidC);

    tEnd = toc(t0);
    fprintf('Done in: %d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));
    
end