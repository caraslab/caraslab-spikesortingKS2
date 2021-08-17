igood = rez.ops.igood;
[chanMap, ~, ~, ~, NchanTOTdefault] = loadChanMap(ops.chanMap); % function to load channel map file
[b1, a1] = butter(3, ops.fshigh/ops.fs*2, 'high'); % the default is to only do high-pass filtering at 150Hz

trange = [ops.trange(1)+60, ops.trange(1)+120];

buff = rawsig(:,ceil(trange(1)*ops.fs:trange(2)*ops.fs));
dataRAW = gpuArray(buff); % move int16 data to GPU
dataRAW = dataRAW';
dataRAW = single(dataRAW); % convert to float32 so GPU operations are fast
% subtract the mean from each channel
dataRAW = dataRAW - mean(dataRAW, 1); % subtract mean of each channel
dataRAW = dataRAW - median(dataRAW(:, chanMap(igood)), 2); % subtract median across channels

datr = filter(b1, a1, dataRAW); % causal forward filter
datr = flipud(datr); % reverse time
datr = filter(b1, a1, datr); % causal forward filter again
datr = flipud(datr); % reverse time back

%Find all the peaks (spikes and artifacts) greater than 20 std above noise (might need to be tweaked)
% stdbkg = median((abs(sig)/0.6745));
abs_sig_median = median(abs(datr), 1);
abs_sig_std = std(abs(datr), [], 1);
thresh = abs_sig_median + 20*abs_sig_std;


ch = 31;
plot(datr(:,ch))
hold on
yline(gather(abs_sig_median(ch)))
yline(gather(thresh(ch)))
xticks(linspace(0, length(datr(:,ch))*ops.fs, length(datr(:,ch))))