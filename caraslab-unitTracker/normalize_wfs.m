function [bestChan_wfs, bestChan_wf_mean, shank, allShankChan_wf_mean, shank_channels] = normalize_wfs(wf, gwfparams, wf_idx)
    % This function levels wfs to be around 0V and rescales so that
    % their peak is at -1 or 1

    cluster_phy_id = wf.unitIDs(wf_idx);
    % Get channel with highest amplitude for this unit

    % This usually works but sometimes phy fails at detecting the best
    % channel, so let's do it manually
%         best_channel = gwfparams.cluster_quality.ch(...
%             gwfparams.cluster_quality.cluster_id == cluster_phy_id);

    temp_wfs = wf.waveFormsMean(wf_idx, :,:);
    temp_wfs = abs(temp_wfs);

    [max_amplitude, ~] = max(temp_wfs, [], 3);
    clear temp_wfs

    % Grab best channel index
    [~, best_channel] = max(max_amplitude);
    best_channel_0in = best_channel - 1;
    best_channel_idx = find(gwfparams.chanMap == best_channel_0in);
    % Get shank
    shank = gwfparams.cluster_quality.sh(...
        gwfparams.cluster_quality.cluster_id == cluster_phy_id);  
    
    % What are the channels in current shank?
    shank_channels_0ind = gwfparams.chanMap(gwfparams.channelShanks == shank);
%     [~ , shank_channels_idx] = intersect(gwfparams.chanMap, shank_channels_0ind);
    shank_channels = shank_channels_0ind + 1;
%     shank_channels = sort(shank_channels);
%     
    wf_npoints = length(wf.waveFormsMean(1, 1,:));
    
    allShankChan_wf_mean = nan(wf_npoints, length(shank_channels_0ind));
    

    % Get best channel and all wf shapes too
    % Get waveforms at current channel
    cur_chan_wfs = wf.waveForms(wf_idx, :, best_channel_idx,:);
    cur_chan_wfs = squeeze(cur_chan_wfs);
    % Find and remove nan wfs
    cur_chan_wfs(any(isnan(cur_chan_wfs), 2),:) = [];	

    % Normalize amplitude so that peak is at -1 or 1
    cur_chan_wf_mean = wf.waveFormsMean(wf_idx, best_channel_idx,:);
    cur_chan_wf_mean = squeeze(cur_chan_wf_mean);

    % Recenter and rescale wfs
    cur_chan_wf_mean = cur_chan_wf_mean - mean(cur_chan_wf_mean(1:5));
    abs_wf = abs(cur_chan_wf_mean);
    peak_idx = find(abs_wf == max(abs_wf));
    peak_value = abs(cur_chan_wf_mean(peak_idx));

    if length(peak_value) > 1
        peak_value = peak_value(1);
    end

    bestChan_wf_mean = cur_chan_wf_mean ./ peak_value;

    bestChan_wfs = cur_chan_wfs ./ peak_value;
    
    % Get means for each channel and normalize by best channel peak
    for ch_dummy_idx=1:length(shank_channels_0ind)
        cur_ch_0ind = shank_channels_0ind(ch_dummy_idx);
        cur_ch_idx = find(gwfparams.chanMap == cur_ch_0ind);
        
        if isempty(cur_ch_idx)
            continue;
        end
        
        % Normalize amplitude so that peak is at -1 or 1
        cur_chan_wf_mean = wf.waveFormsMean(wf_idx, cur_ch_idx,:);
        cur_chan_wf_mean = squeeze(cur_chan_wf_mean);

        % Recenter and rescale wfs
        cur_chan_wf_mean = cur_chan_wf_mean - mean(cur_chan_wf_mean(1:5));
        abs_wf = abs(cur_chan_wf_mean);
        peak_idx = find(abs_wf == max(abs_wf));
        
        peak_value = abs(cur_chan_wf_mean(peak_idx));
        if length(peak_value) > 1
            peak_value = peak_value(1);
        end
        
        cur_chan_wf_mean = cur_chan_wf_mean ./ peak_value;

        allShankChan_wf_mean(:, ch_dummy_idx) = cur_chan_wf_mean;
    end


end