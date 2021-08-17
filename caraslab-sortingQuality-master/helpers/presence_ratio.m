function [pr] = presence_ratio(spike_train, min_time, max_time, num_bins, ops)
%     From Allen Brain Institute repository
%     Translated from Python by M Macedo-Lima,  December 2020

%     """Calculate fraction of time the unit is present within an epoch.
% 
%     Inputs:
%     -------
%     spike_train : array of spike times
%     min_time : minimum time for potential spikes
%     max_time : maximum time for potential spikes
% 
%     Outputs:
%     --------
%     presence_ratio : fraction of time bins in which this unit is spiking
% 
%     """
    
    % MML edit: for concatenated recordings, the presence ratio is
    % confounded by periods of noise that were not removed before sorting
    % So we remove those periods here
    if nargin > 4
        ops = ops;
    else
        ops = [];
    end
    
    if ~isempty(ops)
        % Read the concat_tranges to subtract the max_time from the offsets
        % Also, use the same offset to "correct" spike times and the
        % min_time

        % Subtract the offsets from spike times
        if isfield(ops, 'concat_tranges')
            total_offset = sum(ops.concat_tranges, 1);
            max_time = max_time - total_offset(1);
            for dummy_idx=size(ops.concat_cumulative_tranges, 1):-1:1
                cur_offset_time = ops.concat_cumulative_tranges(dummy_idx, 1);
                cur_offset_duration = ops.concat_tranges(dummy_idx, 1);
                offset_spikes = spike_train > cur_offset_time;  % only spike times higher than current offset need to be shifted
                spike_train(offset_spikes) = spike_train(offset_spikes) - cur_offset_duration;
            end
        else
            max_time = max_time - ops.trange(1);
            spike_train = spike_train - ops.trange(1);
        end
    end
    
    [h, ~] = histcounts(spike_train, linspace(min_time, max_time, num_bins));
    
    pr = sum(h > 0) / num_bins;
end