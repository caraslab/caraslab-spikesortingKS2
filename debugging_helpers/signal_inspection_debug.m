load('config.mat')
fid         = fopen(ops.fbinary, 'r'); % open for reading raw data

offset = round(510*ops.fs);
% igood = ops.igood;
[chanMap, xc, yc, kcoords, NchanTOTdefault] = loadChanMap(ops.chanMap); % function to load channel map file

fseek(fid, offset, 'bof'); % fseek to batch start in raw file
% close all
NTbuff = 400*ops.fs;
buff = fread(fid, [ops.NchanTOT NTbuff], '*int16'); % read and reshape. Assumes int16 data (which should perhaps change to an option)
fclose(fid)
dataRAW = gpuArray(buff); % move int16 data to GPU
dataRAW = dataRAW';

dataRAW = single(dataRAW); % convert to float32 so GPU operations are fast
% subtract the mean from each channel
dataRAW = dataRAW - mean(dataRAW, 1); % subtract mean of each channel

if getOr(ops, 'CAR', 1)
    % MML edit:take median of good channels only
    dataRAW = dataRAW - median(dataRAW(:, chanMap(igood)), 2); % subtract median across channels
end

datr = filter(b1, a1, dataRAW); % causal forward filter

datr = flipud(datr); % reverse time
datr = filter(b1, a1, datr); % causal forward filter again
datr = flipud(datr); % reverse time back
figure
for ch=1:NchanTOT
    plot(datr(:,ch) + ch*50000)
    hold on 
end