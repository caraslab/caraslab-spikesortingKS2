
[b1, a1] = butter(3, [150/24410,3000/24410]*2, 'bandpass'); % butterworth filter with only 3 nodes (otherwise it's unstable for float32)

% MML edit
% Redo the filters without adding offsets for visualization purposes on
% Phy

ops.fbinary = '/media/matheus/132bfc10-ead6-48da-986e-007a5a3d1d87/Matt/Sorted/Sanes lab chanMap test/Aversive-AM-200204-135330/Aversive-AM-200204-135330.dat';
ops.fclean = '/media/matheus/132bfc10-ead6-48da-986e-007a5a3d1d87/Matt/Sorted/Sanes lab chanMap test/Aversive-AM-200204-135330/Aversive-AM-200204-135330_CLEAN.dat';
load('/media/matheus/132bfc10-ead6-48da-986e-007a5a3d1d87/Matt/Sorted/Sanes lab chanMap test/Aversive-AM-200204-135330/rez.mat')
ops = rez.ops;

NchanTOT = ops.NchanTOT;

% ops.ntbuff = 64;
% ops.NT = 64*1024+ ops.ntbuff;
NT = ops.NT;
NTbuff      = NT + 4*ops.ntbuff; % we need buffers on both sides for filtering

% ops.fs = 24410;

ops.trange = [0 Inf]; % time range to sort
bytes       = get_file_size(ops.fbinary); % size in bytes of raw binary
nTimepoints = floor(bytes/NchanTOT/2); % number of total timepoints
ops.tstart  = ceil(ops.trange(1) * ops.fs); % starting timepoint for processing data segment
ops.tend    = min(nTimepoints, ceil(ops.trange(2) * ops.fs)); % ending timepoint
ops.sampsToRead = ops.tend-ops.tstart; % total number of samples to read
ops.twind = ops.tstart * NchanTOT*2; % skip this many bytes at the start

Nbatch      = ceil(ops.sampsToRead /(NT-ops.ntbuff)); % number of data batches


fidC        = fopen(ops.fclean,  'w'); % MML edit; write processed data for phy
fid         = fopen(ops.fbinary, 'r'); % open for reading raw data
for ibatch = 1:Nbatch
    % we'll create a binary file of batches of NT samples, which overlap consecutively on ops.ntbuff samples
    % in addition to that, we'll read another ops.ntbuff samples from before and after, to have as buffers for filtering
%     offset = max(0, ops.twind + 2*NchanTOT*((NT - ops.ntbuff) * (ibatch-1) - 2*ops.ntbuff)); % number of samples to start reading at.
    offset = max(0, ops.twind + 2*NchanTOT*(NT*(ibatch-1) - 2*ops.ntbuff)); % number of samples to start reading at.
    if offset==0
        ioffset = 0; % The very first batch has no pre-buffer, and has to be treated separately
    else
        ioffset = 2*ops.ntbuff;
    end
    fseek(fid, offset, 'bof'); % fseek to batch start in raw file

    buff = fread(fid, [NchanTOT NTbuff], '*int16'); % read and reshape. Assumes int16 data (which should perhaps change to an option)
    
    % DEBUG
%     a_pre = buff(1, 65600-100:65600);
%     b_pre = buff(1, 65601:65601 + 100);
    
    if isempty(buff)
        break; % this shouldn't really happen, unless we counted data batches wrong
    end
    

%     if offset==0
%         buff = [zeros(NchanTOT, ioffset) buff]; % pad with starting zeros on first batch
%     end
%     
    nsampcurr = size(buff,2); % how many time samples the current batch has
    
    if nsampcurr<NTbuff
        buff(:, nsampcurr+1:NTbuff) = repmat(buff(:,nsampcurr), 1, NTbuff-nsampcurr); % pad with zeros, if this is the last batch
    end
% 
    dataRAW = gpuArray(buff); % move int16 data to GPU
    dataRAW = single(dataRAW); % convert to float32 so GPU operations are fast
    dataRAW = dataRAW';  % transpose for filtering
    datr = filter(b1, a1, dataRAW); % causal forward filter

    datr = flipud(datr); % reverse time
    datr = filter(b1, a1, datr); % causal forward filter again
    datr = flipud(datr); % reverse time back
    
%     datr    = gpufilter(buff, ops, ops.chanMap(ops.igood), 1); % apply filters and median subtraction
    % FOR DEBUG
%     datr = dataRAW;
 
    datr    = datr(ioffset + (1:NT),:); % remove timepoints used as buffers
%     if offset ==0
%         datr    = datr((1:NT),:); % remove timepoints used as buffers
%     else
%         datr    = datr( (NT - offset/NchanTOT/2) + (1:NT),:); % remove timepoints used as buffers
%     end
    % DEBUG
%     a = datr(end-100:end, 1);    
% %     b = datr(1:101, 1);
%     figure
%     hold on
%     plot([a_pre b_pre])
%     plot([a; b])
%     
    datr = datr';
    
    datr  = gather(int16(datr)); % convert to int16, and gather on the CPU side
    fwrite(fidC, datr(:), 'int16'); % write this batch to clean file
end
fclose(fid); % close the files
fclose(fidC);

% 
% 
% % set up the parameters of the filter
% if isfield(ops,'fslow')&&ops.fslow<ops.fs/2
%     [b1, a1] = butter(3, [ops.fshigh/ops.fs,ops.fslow/ops.fs]*2, 'bandpass'); % butterworth filter with only 3 nodes (otherwise it's unstable for float32)
% else
%     [b1, a1] = butter(3, ops.fshigh/ops.fs*2, 'high'); % the default is to only do high-pass filtering at 150Hz
% end
% 
% if getOr(ops, 'comb', 0)  % MML edit; comb filter
%     N  = 407;    % Order
%     BW = 2;    % Bandwidth
%     Fs = ops.fs;  % Sampling Frequency
%     h = fdesign.comb('Notch', 'N,BW', N, BW, Fs);
%     comb_filter = design(h, 'butter');
%     comb_b1= comb_filter.Numerator;
%     comb_a1= comb_filter.Denominator;
% end
% 
% % Can't use GPU acceleration for comb filter yet...
% buff = buff';  % MML edit: transpose sooner
% if getOr(ops, 'comb', 0)  % MML edit; comb filter
%     buff = filter(comb_b1, comb_a1, buff);
%     dataRAW = gpuArray(buff); % move int16 data to GPU
% else
%     dataRAW = gpuArray(buff); % move int16 data to GPU
% 
% end
% % MML edit: reorganized this for comb filter
% %     dataRAW = dataRAW';
% %     dataRAW = dataRAW(:, chanMap); % subsample only good channels 
% dataRAW = single(dataRAW); % convert to float32 so GPU operations are fast
% % MML edit: even bad channels are filtered but excluded outside of this
% % function; except for CAR. See below
% % dataRAW = dataRAW(:, chanMap); % subsample only good channels  
% 
% % subtract the mean from each channel
% dataRAW = dataRAW - mean(dataRAW, 1); % subtract mean of each channel
% 
% % CAR, common average referencing by median
% if getOr(ops, 'CAR', 1)
% %     dataRAW = dataRAW - median(dataRAW, 2); % subtract median across channels
%     % MML edit:take median of good channels only
%     dataRAW = dataRAW - median(dataRAW(:, chanMap), 2); % subtract median across channels
% end
% 
% % next four lines should be equivalent to filtfilt (which cannot be used because it requires float64)
% datr = filter(b1, a1, dataRAW); % causal forward filter
% 
% datr = flipud(datr); % reverse time
% datr = filter(b1, a1, datr); % causal forward filter again
% datr = flipud(datr); % reverse time back