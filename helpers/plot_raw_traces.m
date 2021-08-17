
cluster_n = 3288;

recording_name = 'SUBJ-ID-26_200716_concat';
input_dir = '/mnt/132bfc10-ead6-48da-986e-007a5a3d1d87/Matt/Sorted/SUBJ-ID-26-200614-103221/200716_concat';

input_file = [input_dir '/' '200716_concat_CLEAN300hz.dat'];


output_file = [input_dir '/' recording_name '_cluster' int2str(cluster_n) '_trace.txt'];

gwfparams.cluster_quality = tdfread(fullfile(input_dir, 'cluster_info.tsv'));

gwfparams.chanMap = readNPY(fullfile(input_dir, 'channel_map.npy')); % this is important in esp if you got rid of files. 

best_channel = gwfparams.cluster_quality.ch(gwfparams.cluster_quality.cluster_id == cluster_n);  
best_channel_idx = find(gwfparams.chanMap == best_channel);

pre_chunk_s = 0;
post_chunk_s = 5;

n_channel = 64;

sr = 24414.0625;

nt = ceil((post_chunk_s + pre_chunk_s) * sr);



cur_time = 1390;

offset_bytes = ceil((cur_time - pre_chunk_s) * sr * 2 * n_channel);

fo = fopen(input_file);
fseek(fo, offset_bytes, 'bof');
cur_buff = fread(fo, [n_channel nt], '*int16');
fclose(fo);

snip_to_plot = gpuArray(cur_buff(best_channel_idx, :));
x_time = gpuArray(linspace(-pre_chunk_s, post_chunk_s, length(snip_to_plot)));
figure
plot(x_time, snip_to_plot, 'black')
print([input_dir '/' recording_name '_cluster' int2str(cluster_n) '_trace'], '-dpdf', '-bestfit', '-painters');

dlmwrite(output_file, snip_to_plot);
