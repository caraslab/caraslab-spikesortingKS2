% Quick n' dirty script to output voltage recordings by channel and color
% them according to timestamps detected through unit sorting

clusters_to_plot = [51];
cluster_colors = {'black'};
channels_to_plot = [1];

nChanTOT = 4;  % channels in probe

cur_time = 356;

sr = 30000;

% Before and after time chunk
pre_chunk_s = 1;
post_chunk_s = 20;

% Before and after waveform timestamp
pre_wf_s = 0.001;
post_wf_s = 0.003;

% IO paths to find your file
subj_id = 'SUBJ-ID-84';
recording_name = '210813_concat_L2770';
input_dir = '/mnt/CL_4TB_2/Matt/Acute_opto_recordings/SUBJ-ID-84/210813_concat_L2770';
input_file = [input_dir filesep recording_name '_CLEAN300hz.dat'];

% Key file to shade around stimulus, e.g. optogenetic stimulation
key_file = '/mnt/CL_4TB_2/Matt/Acute_opto_recordings/SUBJ-ID-84/2021-08-13_13-26-07_L2770/CSV files/2021-08-13_13-26-07_L2770_DAC1.csv';
key_file = readtable(key_file);

shading_color = [245, 114, 66]/255;  % RGB amber
shading_color = [107, 159, 198]/255; % RGB blue

do_shading = 1;  % Set to 0 if no stimulus shading is wanting

% Shouldn't need to tweak anything from here down
gwfparams.cluster_quality = tdfread(fullfile(input_dir, 'cluster_info.tsv'));

gwfparams.chanMap = readNPY(fullfile(input_dir, 'channel_map.npy')); % this is important in esp if you got rid of files. 

nt = ceil((post_chunk_s + pre_chunk_s) * sr);

offset_bytes = ceil((cur_time - pre_chunk_s) * sr * 2 * nChanTOT);

fo = fopen(input_file);
fseek(fo, offset_bytes, 'bof');
cur_buff = fread(fo, [nChanTOT nt], '*int16');
fclose(fo);

figure
scaling_factor = 500;
for ch_idx=1:length(channels_to_plot)
    ch_to_plot = channels_to_plot(ch_idx);
    snip_to_plot = cur_buff(ch_to_plot, :);

    h = plot(1:length(snip_to_plot), scaling_factor*ch_idx+snip_to_plot, 'black');
    
    hold on
    
    % Color wfs
    for cluster_idx=1:length(clusters_to_plot)
        cluster = clusters_to_plot(cluster_idx);
        timestamp_file = [input_dir filesep subj_id '_' recording_name '_cluster' int2str(cluster) '.txt'];
        cluster_timestamps = fscanf(fopen(timestamp_file, 'r'),'%f');
        relevant_timestamps = cluster_timestamps(...
                                cluster_timestamps > cur_time - pre_chunk_s & ...
                                cluster_timestamps < cur_time + post_chunk_s);
        cur_color = cluster_colors{cluster_idx};
        for timestamp_idx=1:length(relevant_timestamps)
            timestamp = relevant_timestamps(timestamp_idx);

            x_time_idx_start = floor((timestamp - pre_wf_s - (cur_time - pre_chunk_s))*sr)+1;
            
            x_time_idx_end = x_time_idx_start + ceil((pre_wf_s + post_wf_s)*sr);
            if x_time_idx_end-1 > length(snip_to_plot)
                x_time_idx_end = length(snip_to_plot);
            end
            
            plot(x_time_idx_start:x_time_idx_end-1, scaling_factor*ch_idx+snip_to_plot(x_time_idx_start:x_time_idx_end-1), 'Color', rgb(cur_color));
        end
    end
    
    if do_shading
        % Apply shading during stims
        relevant_timestamps = key_file((key_file.Onset > cur_time - pre_chunk_s) & ...
                                        (key_file.Onset < cur_time + post_chunk_s),:);


        fill_YY = [min(scaling_factor*ch_idx+snip_to_plot), max(scaling_factor*ch_idx+snip_to_plot)];

        YY = repelem(fill_YY, 1, 2);  
        for event_idx=1:height(relevant_timestamps)
            fill_XX = ([relevant_timestamps(event_idx,:).Onset relevant_timestamps(event_idx,:).Offset] - cur_time + pre_chunk_s)*sr;
            XX = [fill_XX, fliplr(fill_XX)];
            f = fill(XX, YY, shading_color);
            % Choose a number between 0 (invisible) and 1 (opaque) for facealpha.  
            set(f,'facealpha',.3,'edgecolor','none')
        end
    end
    
    % Save as .wav file too for fun
%     Fs = sr;                                       % Sampling Frequency (Hz)
%     Fn = Fs/2;                                              % Nyquist Frequency (Hz)
%     Wp = 8000/Fn;                                           % Passband Frequency (Normalised)
%     Ws = 8500/Fn;                                           % Stopband Frequency (Normalised)
%     Rp =   1;                                               % Passband Ripple (dB)
%     Rs = 150;                                               % Stopband Ripple (dB)
%     [n,Ws] = cheb2ord(Wp,Ws,Rp,Rs);                         % Filter Order
%     [z,p,k] = cheby2(n,Rs,Ws,'low');                        % Filter Design
%     [soslp,glp] = zp2sos(z,p,k);                            % Convert To Second-Order-Section For Stability
%     filtered_sound = filtfilt(soslp, glp, double(snip_to_plot));
    
    audio_data = rescale(resample(double(snip_to_plot), 44100, sr), -1, 1);
    audiowrite([input_dir '/' recording_name '_sampleTrace_ch' int2str(ch_to_plot) '.wav'], audio_data, 44100)
                                
end
% Transform xticklabels from samples to time
ax = ancestor(h, 'axes');
ax.XAxis.Exponent = 0;
ax.XAxis.TickLabelFormat = '%.0f';
xtl = cell2mat(cellfun(@str2num,xticklabels,'UniformOutput',false));
xticklabels([round(xtl / sr, 2)]);
print([input_dir '/' recording_name '_sampleTrace'], '-dpdf', '-bestfit', '-painters');



