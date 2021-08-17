input_file = '/mnt/132bfc10-ead6-48da-986e-007a5a3d1d87/Matt/Sorted/SUBJ-ID-26-200614-103221/200720_concat/200720_concat_CLEAN.dat';

best_channel = 20;

n_channel = 64;

sr = 24410;
plot_time = 300;
nt = plot_time * sr;

offset_time = 0;

fo = fopen(input_file);

% Read response from first playback
% for cur_idx, cur_time in enumerate(key_df[key_df.iloc[:, 1] == wav_file].iloc[:, 0]):

offset_bytes = ceil(offset_time * sr * 2 * n_channel);

fseek(fo, offset_bytes, 'bof');
cur_buff = fread(fo, [n_channel nt], '*int16');
fclose(fo);

Y = cur_buff(best_channel+1, :);

Y = Y';
fprintf('Bandpass filtering cleaned data...\n')

% Y = Y - mean(Y, 1);

% bandpass and comb
% hp = 150;   %High pass (Hz)
% lp = 3000;  %Low pass (Hz)
% [b1, a1] = butter(3, [hp/sr,lp/sr]*2, 'bandpass');
% Y = filter(b1, a1, Y);
% Y = flipud(Y);
% Y = filter(b1, a1, Y);
% Y = flipud(Y);
% 

N  = 407;    % Order
BW = 2;    % Bandwidth
Fs = sr;  % Sampling Frequency
h = fdesign.comb('Notch', 'N,BW', N, BW, Fs);
comb_filter = design(h, 'butter');
comb_b1= comb_filter.Numerator;
comb_a1= comb_filter.Denominator;
comb_filter = design(h, 'butter');

Y_filt = filter(comb_b1, comb_a1, Y);
Y= Y';
Y_filt = Y_filt';

% 
% Y_sine = chunkwiseDeline(Y', sr, [60, 180, 300], 10);
% Y_sine = Y_sine';

fft_Y = fft(Y);
fft_Y_filt = fft(Y_filt);
% fft_Y_sine = fft(Y_sine);

n = length(Y);          % number of samples
f = (0:n-1)*(sr/n);     % frequency range

power_Y = abs(fft_Y).^2/n;    % power of the DFT
power_Y_filt = abs(fft_Y_filt).^2/n;    % power of the DFT
% power_Y_sine = abs(fft_Y_sine).^2/n;    % power of the DFT

plot(f,power_Y)
hold on
plot(f,power_Y_filt)
% hold on
% plot(f,power_Y_sine)
xlabel('Frequency')
ylabel('Power')

figure
plot(Y)
hold on
plot(Y_filt)
