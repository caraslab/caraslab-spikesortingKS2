function [signal, thresholds, exit_flag] = caraslab_artifact_reject(signal, ops, inspect_results)
warning('off','all')
%[finalsignal,rejectthresh] = caraslab_artifact_reject(rawsignal,fs)
%
% This function performs automatic artifact removal on raw in vivo
% extracellular physiology data. Artifact removal proceeds as follows:
% First, all peaks (both spikes and artifacts) with amplitudes greater than
% 5 median-based stdevs of the baseline noise are identified. Next, the
% median-based stdev of the amplitude of those peaks is calculated, and a
% rejection threshold is set to 2X that value. Any peaks with amplitudes
% exceeding this rejection threshold are considered to be artifacts. A
% window is set around each artifact peak, and the signal is set to zero
% within that window. This initial pass is somewhat conservative in 
% identifying artifacts. Therefore, a second pass is run, where the peak
% identifcation is performed on the initial cleaned signal, and remaining
% artifcats are identified and zeroed. 
%
% Input variables:
%   rawsignal:  MxN matrix of raw voltage values, where
%                   M = number of samples, and 
%                   N = number of channels
%
%   fs:         Sampling rate of data (in Hz)
%
% Output variables:
%   finalsignal: MxN matrix of cleaned voltage values,where 
%                   M = number of samples, and 
%                   N = number of channels. 
%                Note that cleaned signal is high-pass filtered at 100 Hz.
%
%   rejectthresh: artifact rejection threshold (V)
%
%
% Written by ML Caras Apr 1 2019

%Determine the number of channels
numchans = ops.NchanTOT;


%Define the number of samples before and after each artifact peak to be
%normalized. This value was empirically determined.
win = 0.01*ops.fs;


% Save original signal here in case it needs to  be recovered after
% inspection
backup_signal = signal;

% Threshold in multiples of std
std_threshold = ops.std_threshold;

exit_flag = 0;


%---------------------------------------------------------------------
%First cleaning pass
% fprintf('First pass artifact rejection...\n')
[signal, thresholds] = rm_artifacts(signal,numchans, win, std_threshold);

%---------------------------------------------------------------------
%Second cleaning pass
%  fprintf('Second pass artifact rejection...\n')
[signal, thresholds] = rm_artifacts(signal,numchans, win, std_threshold);

if inspect_results
    [new_threshold, exit_flag] = inspect_plot(backup_signal, signal, thresholds, ops);
    
    fprintf('New threshold chosen: %f. Change the config file and rerun caraslab_preprocessdat with inspect_results = 0 to apply it\n', new_threshold)
end


% TODO: turn this into a parfor function?
% TODO2: after median and std are taken, it can be chunked and GPU-based...
function [ret_sig, thresh] = rm_artifacts(sig,numchans,win, std_threshold)
    ret_sig = sig;
    %Find all the peaks (spikes and artifacts) greater than ops.std_threshold std above noise (might need to be tweaked)
    % stdbkg = median((abs(sig)/0.6745));
    abs_sig = abs(sig);
    
    % Find local medians and stds first then take the grand median of all
    % values; this way, brief moments of noise will not affect the true
    % signal median and std
    n_blocks = 100;
    block_size = size(sig, 1) / n_blocks;
    median_by_block = nan(n_blocks, numchans);
    std_by_block = nan(n_blocks, numchans);
    abs_median_by_block = nan(n_blocks, numchans);
    abs_std_by_block = nan(n_blocks, numchans);
    block_counter = 1;
    for block_start = 1:block_size:size(sig, 1)
        % If it's the last block
        if block_start+block_size>size(sig,1)
            block_size = size(sig,1) - block_start;
        end
        
        median_by_block(block_counter, :) = median(sig(block_start:block_start+block_size, :), 1);
        std_by_block(block_counter, :) = std(sig(block_start:block_start+block_size, :), [], 1);  % better to use mad() maybe?
        abs_median_by_block(block_counter, :) = median(abs_sig(block_start:block_start+block_size, :), 1);
        abs_std_by_block(block_counter, :) = std(abs_sig(block_start:block_start+block_size, :), [], 1);
        block_counter = block_counter + 1;
    end
    
    % Grand medians across blocks
    sig_median = median(median_by_block, 1);
    sig_std = median(std_by_block, 1);
    abs_sig_median = median(abs_median_by_block, 1);
    abs_sig_std = median(abs_std_by_block, 1);
    
    thresh = abs_sig_median + std_threshold*abs_sig_std;

    for ch = 1:numchans
        cleansig_ch = sig(:,ch);
%         cleansig_ch(abs(cleansig_ch) > thresh(ch)) = sig_std(ch).*randn() + sig_median(ch);
        peak_indeces = find(abs(cleansig_ch) > thresh(ch) == 1);
%         [p,l]=findpeaks(abs(cleansig_ch),ops.fs,'MinPeakHeight', thresh(ch),'MinPeakDistance',0.05);
        for i = 1:win/2:numel(peak_indeces)
            samp = peak_indeces(i);
            try
%                 cleansig_ch(samp-win+1:samp+win) = zeros(length(cleansig_ch(samp-win+1:samp+win)), 1);
                cleansig_ch(samp-win+1:samp+win) = sig_std(ch).*randn(length(cleansig_ch(samp-win+1:samp+win)), 1) + sig_median(ch);
            catch ME
                if strcmp(ME.identifier, 'parallel:gpu:array:InvalidValue')
                    if samp-win < 0
%                         cleansig_ch(1:samp+win) = zeros(length(cleansig_ch(1:samp+win)), 1);
                        cleansig_ch(1:samp+win) = sig_std(ch).*randn(length(cleansig_ch(1:samp+win)),1) + sig_median(ch);

                    elseif samp+win > length(cleansig_ch)
%                         cleansig_ch(samp-win+1:end) = zeros(length(cleansig_ch(samp-win+1:end)), 1);
                        cleansig_ch(samp-win+1:end) = sig_std(ch).*randn(length(cleansig_ch(samp-win+1:end)),1) + sig_median(ch);
                    end
                end
            end
        end
    ret_sig(:,ch) = cleansig_ch;
    end
    

end


function [new_threshold, exit_flag] = inspect_plot(sig, clean_sig, thresholds, ops)
    % Decimate for faster visualization?
    
    cur_threshold = ops.std_threshold;
    exit_flag = 0;
    time_vec = linspace(0, size(sig, 1)/ops.fs, size(sig, 1));
    
    f = gcf;
    set(f, 'Position', get(0, 'Screensize')/2);
    
    hold on;
    channel_to_plot = 1;
    cax1 = plot(time_vec, sig(:, channel_to_plot), 'color', 'black');
    cax2 = plot(time_vec, clean_sig(:, channel_to_plot), 'color', 'blue');
    cax3 = yline(thresholds(channel_to_plot), ':r', 'LineWidth', 3);
    cax4 = yline(-thresholds(channel_to_plot), ':r', 'LineWidth', 3);
    cax1.Color(4) = 0.2; cax2.Color(4) = 0.2;

    % Threshold slider
    b = uicontrol('Parent',f,'Style','slider','Position',[81,54,419,23],...
              'value',cur_threshold, 'min',0, 'max',200);
    bgcolor = f.Color;
    bl1 = uicontrol('Parent',f,'Style','text','Position',[50,54,23,23],...
                    'String','0','BackgroundColor',bgcolor);
    bl2 = uicontrol('Parent',f,'Style','text','Position',[500,54,23,23],...
                    'String','200','BackgroundColor',bgcolor);
    bl3 = uicontrol('Parent',f,'Style','text','Position',[240,25,100,23],...
                    'String','Threshold','BackgroundColor',bgcolor);
                
    addlistener(b,'ContinuousValueChange',@(hObject, event) changeThresholdAndPlot(hObject, event));
    
    function changeThresholdAndPlot(hObject,event)
        nthreshold = get(hObject,'Value');
        nops = ops;
        nops.NchanTOT = 1;
        nops.std_threshold = nthreshold;
        [new_y, new_thresh] = caraslab_artifact_reject(sig(:, channel_to_plot), nops, 0);
        set(cax2,'ydata', gather(new_y));
        set(cax3,'Value', new_thresh);
        set(cax4,'Value', -new_thresh);
        drawnow;
    end
    
    % Continue button
    continue_button = uicontrol('Parent',f,'String','Continue','Callback',@OnContinueButtonCallback, 'Position',[0,100,100,30]);
    function OnContinueButtonCallback(~,~)
        uiresume(f);
    end
    
    % Exit button
    exit_button = uicontrol('Parent',f,'String','Exit','Callback',@OnExitButtonCallback, 'Position',[0,50,100,30]);
    function OnExitButtonCallback(~,~)
        uiresume(f);
        exit_flag = 1;
    end
    
    % Channel selector
    ddLabel = uicontrol('Parent',f,'Style','text','Position',[0,230,100,50],...
                    'String','Channel:','BackgroundColor',bgcolor);
    dd = uicontrol('Parent', f, 'Style','popupmenu', 'String', strsplit(num2str(1:size(sig, 2))),...
                     'Value',1, 'Position',[30,200,50,50]);

    addlistener(dd,'Action',@(hObject, event) changeChannelAndPlot(hObject, event));
    function changeChannelAndPlot(hObject,event, hplot1, hplot2)
        channel_to_plot = get(hObject,'Value');
        set(cax1,'ydata', gather(sig(:, channel_to_plot)));
        set(cax2,'ydata', gather(clean_sig(:, channel_to_plot)));
    end
    
    uiwait(f)
    
    new_threshold = get(b,'Value');
    
    close(f)
end
end
