function [spike_times, spike_clusters, spike_templates, amplitudes, pc_features] = ...
        remove_spikes(spike_times, spike_clusters, spike_templates, amplitudes, pc_features, spikes_to_remove)
% 
%         """
%         Removes spikes from Kilosort outputs
% 
%         Inputs:
%         ------
%         spike_times : numpy.ndarray (num_spikes x 0)
%             Spike times in samples 
%         spike_clusters : numpy.ndarray (num_spikes x 0)
%             Cluster IDs for each spike time
%         spike_templates : numpy.ndarray (num_spikes x 0)
%             Template IDs for each spike time
%         amplitudes : numpy.ndarray (num_spikes x 0)
%             Amplitude value for each spike time
%         pc_features : numpy.ndarray (num_spikes x num_pcs x num_channels)
%             Pre-computed PCs for blocks of channels around each spike
%         spikes_to_remove : numpy.ndarray
%             Indices of spikes to remove
% 
%         Outputs:
%         --------
%         spike_times : numpy.ndarray (num_spikes - spikes_to_remove x 0)
%         spike_clusters : numpy.ndarray (num_spikes - spikes_to_remove x 0)
%         spike_templates : numpy.ndarray (num_spikes - spikes_to_remove x 0)
%         amplitudes : numpy.ndarray (num_spikes - spikes_to_remove x 0)
%         pc_features : numpy.ndarray (num_spikes - spikes_to_remove x num_pcs x num_channels)
% 
%         """
        spike_times(spikes_to_remove) = [];
        spike_clusters(spikes_to_remove) = [];
        spike_templates(spikes_to_remove) = [];
        amplitudes(spikes_to_remove) = [];
        pc_features(spikes_to_remove, :, :) = [];
     
    end