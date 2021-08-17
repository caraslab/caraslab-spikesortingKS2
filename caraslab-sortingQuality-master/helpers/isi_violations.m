
function [fpRate, num_violations] = isi_violations(spike_train, min_time, max_time, isi_threshold, min_isi, ops)
%     From Allen Brain Institute repository
%     Translated from Python by M Macedo-Lima,  December 2020

%     """Calculate ISI violations for a spike train.
% 
%     Based on metric described in Hill et al. (2011) J Neurosci 31: 8699-8705
% 
%     modified by Dan Denman from cortex-lab/sortingQuality GitHub by Nick Steinmetz
% 
%     Inputs:
%     -------
%     spike_train : array of spike times
%     min_time : minimum time for potential spikes
%     max_time : maximum time for potential spikes
%     isi_threshold : threshold for isi violation
%     min_isi : threshold for duplicate spikes
% 
%     Outputs:
%     --------
%     fpRate : rate of contaminating spikes as a fraction of overall rate
%         A perfect unit has a fpRate = 0
%         A unit with some contamination has a fpRate < 0.5
%         A unit with lots of contamination has a fpRate > 1.0
%     num_violations : total number of violations
% 
%     """

    % MML edit: for concatenated recordings, the presence ratio is
    % confounded by periods of noise that were not removed before sorting
    % (consider removing these in the future). Because of this we need to
    % revisit the file before concatenation to fibure out where the
    % non-noisy recording starts. So we need to run this for every file
    % then compile the results
    if nargin > 5
        ops = ops;
    else
        ops = [];
    end
    
    if ~isempty(ops)

        if isfield(ops, 'concat_tranges')
            % Read the concat_tranges to subtract the max_time from the offsets
            trange_or_tranges = ops.concat_tranges;
        else
            % for non-concatenated files, just read the offset
            trange_or_tranges = ops.trange;
        end
             
        
        for row=1:size(trange_or_tranges, 1)
            max_time = max_time - trange_or_tranges(row, 1);
        end
    end
    
    duplicate_spikes = spike_train(diff(spike_train) <= min_isi);
    
    spike_train = spike_train(~ismember(spike_train, duplicate_spikes));
    
    isis = diff(spike_train);

    num_spikes = length(spike_train);
    num_violations = sum(isis < isi_threshold);
    violation_time = 2*num_spikes*(isi_threshold - min_isi);
    total_rate = num_spikes / (max_time - min_time);
    violation_rate = num_violations/violation_time;
    fpRate = violation_rate/total_rate;
    
end
