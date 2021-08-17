function [fraction_missing, h, b, G] = amplitude_cutoff(wfs)
%   From the AllenInstitude github. 
%   Translated from Python by M Macedo-Lima, 11/30/2020

%     """ Calculate approximate fraction of spikes missing from a distribution of amplitudes
%     Assumes the amplitude histogram is symmetric (not valid in the presence of drift)
%     Inspired by metric described in Hill et al. (2011) J Neurosci 31: 8699-8705
    
    abs_wfs = abs(wfs);
    amplitudes = nan(size(wfs, 1), 1);
    for dummy_idx=1:size(wfs, 1)
        cur_amp = wfs(dummy_idx, abs_wfs(dummy_idx,:) == max(abs_wfs(dummy_idx,:), [], 2));
        % Rarely a waveform is made of NaNs for some unknown reason; Just
        % skip it
        if isempty(cur_amp)
            continue
        end
        cur_amp = cur_amp(1) - mean(wfs(dummy_idx, 1:10));
        amplitudes(dummy_idx) = cur_amp;
    end
    
    % MML edit:
    % Occasionally a waveform will be inverted; check the mean of
    % amplitudes and if it's >0, it probably means the waveform peak is
    % above 0, so invert the amplitude values;
    if mean(amplitudes) > 0
        amplitudes = -amplitudes;
    end
    
    [h, b] = histcounts(amplitudes, 500, 'Normalization', 'pdf');

    pdf = conv(h,  gausswin(3));
    support = b(1:end-1);

    peak_index = find(pdf == max(pdf));

    % Find where minimum values are after the peak
    % MML edit: truncate at zero
    [~, zero_bin] = min(abs(b));
    abs_pdf = abs(pdf(peak_index:end) - pdf(1));
    G = min(find(abs_pdf  == min(abs_pdf), 1) + peak_index, zero_bin);

    bin_size = mean(diff(support));

    fraction_missing = sum(pdf(G:end))*bin_size;
    
    % If more than 50% of spikes are missing, an accurate estimate isn't possible
    fraction_missing = min([fraction_missing, 0.5]);
end