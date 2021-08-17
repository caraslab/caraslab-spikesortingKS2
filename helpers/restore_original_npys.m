function restore_original_npys(Savedir, rez)
    spike_times = rez.spike_times;
    spike_clusters = rez.spike_clusters;
    spike_templates = rez.spike_templates;
    amplitudes = rez.amplitudes;
    pc_features = rez.pc_features;
        
    writeNPY(spike_times, fullfile(Savedir, 'spike_times.npy'));
    writeNPY(spike_clusters, fullfile(Savedir, 'spike_clusters.npy'));
    writeNPY(spike_templates, fullfile(Savedir, 'spike_templates.npy'));
    writeNPY(amplitudes, fullfile(Savedir, 'amplitudes.npy'));
    writeNPY(pc_features, fullfile(Savedir, 'pc_features.npy'));
end