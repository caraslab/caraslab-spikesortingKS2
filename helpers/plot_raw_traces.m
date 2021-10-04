
cluster_n = 1355;

recording_name = 'SUBJ-ID-231_210709_concat';
input_dir = '/mnt/CL_4TB_2/Matt/OFC_PL_recording/Sorting/SUBJ-ID-231/210709_concat';

input_file = dir(fullfile(input_dir, '*_CLEAN300hz.dat'));

input_file = fullfile(input_file.folder, input_file.name);

output_file = [input_dir '/' recording_name '_cluster' int2str(cluster_n) '_trace.txt'];

gwfparams.cluster_quality = tdfread(fullfile(input_dir, 'cluster_info.tsv'));

gwfparams.chanMap = readNPY(fullfile(input_dir, 'channel_map.npy')); % this is important in esp if you got rid of files. 

best_channel = gwfparams.cluster_quality.ch(gwfparams.cluster_quality.cluster_id == cluster_n);  
best_channel_idx = find(gwfparams.chanMap == best_channel);

pre_chunk_s = 2;
post_chunk_s = 2;

n_channel = 64;

sr = 30000;

nt = ceil((post_chunk_s + pre_chunk_s) * sr);



cur_time = 959.9;

offset_bytes = ceil((cur_time - pre_chunk_s) * sr * 2 * n_channel);

fo = fopen(input_file);
fseek(fo, offset_bytes, 'bof');
cur_buff = fread(fo, [n_channel nt], '*int16');
fclose(fo);

snip_to_plot = double(cur_buff(best_channel_idx, :));
% downsample 100x
N=30;
snip_to_plot = gpuArray(interp1(1:length(snip_to_plot), snip_to_plot, linspace(1,length(snip_to_plot), length(snip_to_plot)/N + 1)));

x_time = gpuArray(linspace(-pre_chunk_s, post_chunk_s, length(snip_to_plot)));
figure
plot(x_time, snip_to_plot, 'black')
print([input_dir '/' recording_name '_cluster' int2str(cluster_n) '_trace'], '-dpdf', '-bestfit', '-painters');

dlmwrite(output_file, snip_to_plot);
