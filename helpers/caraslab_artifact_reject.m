function [signal] = caraslab_artifact_reject(signal, std_threshold)
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
numchans = size(signal,2);

%Define the number of samples before and after each artifact peak to be
%zeroed. This value was empirically determined.
win = 244; 

%Create initial highpass filter parameters
% hp = 100;  %(Hz)
% [b1, a1] = butter(3,(hp/fs)*2, 'high');
% 
% %Highpass filter the raw signal
% fprintf('Highpass filtering signal...\n')
% sigflt = filter(b1, a1, rawsignal);
% sigflt = flipud(sigflt);
% sigflt = filter(b1, a1, sigflt);
% sigflt = flipud(sigflt);


%---------------------------------------------------------------------
%First cleaning pass
% fprintf('First pass artifact rejection...\n')
signal = rm_artifacts(signal,numchans,win, std_threshold);

%---------------------------------------------------------------------
%Second cleaning pass
%  fprintf('First pass artifact rejection...\n')
signal = rm_artifacts(signal,numchans,win, std_threshold);

% TODO: turn this into a parfor function?
% TODO2: after median and std are taken, it can be chunked and GPU-based...
function [sig] = rm_artifacts(sig,numchans,win, std_threshold)
    %Find all the peaks (spikes and artifacts) greater than 50 std above noise (might need to be tweaked)
    % stdbkg = median((abs(sig)/0.6745));
    abs_sig_median = median(abs(sig), 1);
    abs_sig_std = std(abs(sig), [], 1);
    thresh = abs_sig_median + std_threshold*abs_sig_std;
    
    % To replace the artifacts with median signal noise
    sig_median = median(sig, 1);
    sig_std = std(sig, [], 1);

    for ch = 1:numchans
        cleansig_ch = sig(:,ch);
%         cleansig_ch(abs(cleansig_ch) > thresh(ch)) = sig_std(ch).*randn() + sig_median(ch);
        peak_indeces = find(abs(cleansig_ch) > thresh(ch) == 1);
        for i = 1:numel(peak_indeces)
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
        end
        sig(:,ch) = cleansig_ch;

        
        
        
        
%         cleansig_ch = sig(:,ch);
%     %     fprintf('Cleaning channel %d\n',ch)
%         [pks,locs] = findpeaks(abs(cleansig_ch),'MinPeakHeight',thresh(ch));
% 
%         %Determine threshold for artifact rejection
%     %     rejectthresh = 2*median(pks/0.6745);
% %         rejectthresh = 2*median(pks);
%         %Find violations
% %         idx = find(pks>rejectthresh);
% 
%         %Define a window of samples around each violation and set signal to 0
% %         for i = 1:numel(idx)
%         for i = 1:numel(locs)
% 
% %             samp = locs(idx(i));
%             samp = locs(i);
%             try
% %                 cleansig_ch(samp-win:samp+win) = 0;
%                 cleansig_ch(samp-win+1:samp+win) = sig_std(ch).*randn(length(cleansig_ch(samp-win+1:samp+win)), 1) + sig_median(ch);
%             catch ME
%                 if strcmp(ME.identifier, 'parallel:gpu:array:InvalidValue')
%                     if samp-win < 0
% %                         cleansig_ch(1:samp+win) = 0;
%                         cleansig_ch(1:samp+win) = sig_std(ch).*randn(length(cleansig_ch(1:samp+win)),1) + sig_median(ch);
%                     elseif samp+win > length(cleansig_ch)
% %                         cleansig_ch(samp-win:end) = 0;
%                         cleansig_ch(samp-win+1:end) = sig_std(ch).*randn(length(cleansig_ch(samp-win+1:end)),1) + sig_median(ch);
%                     end
%                 end
%             end
%         end
%         sig(:,ch) = cleansig_ch;
%     if nargout>1
%         varargout{1} = rejectthresh;
%     end
end

warning('on','all')
end
