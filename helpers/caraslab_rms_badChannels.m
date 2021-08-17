function  [igood] = caraslab_rms_badChannels(sig)
%refsig = caraslab_CAR(sig,nchans,badchans)
%
% Common average referencing is implemented using matrix multiplication
% in the format D x C = R, where:
%
% C is a single time-slice of data in the shape [n x 1]. In other
% words, it is the value from all n channels sampled at a single point
% in time.
%
% D is a n x n matrix. Each row in the matrix defines the weights of the
% individual channels.
%
% R is the referenced output in the shape [n x 1].
%
% Input variables:
%   sig:    MxN matrix of voltage signals, where 
%               M = number of samples, and  
%               N = number of channels.
%           Data should already have had artifacts removed using 
%           caraslab_artifact_reject.m, and have been bandpass filtered. 
%
%   nchans: total number of recording channels (both active and dead)
%
%   badchans: identity of channels that aren't connected because 
%             they're used for wireless transmission
%
%
%Written by ML Caras Apr 1 2019
% Patched by M Macedo-Lima 9/8/20

%--------------------------------------------------------------------------


    %Identify bad channels based on RMS of signal. Limits are taken directly
    %from Ludwig et al. 2009.
    fprintf('Detecting bad channels...\n')
    tic
    rmsvals = rms(sig, 2);
    minrms = 0.3*mean(rmsvals);
    maxrms = 2*mean(rmsvals);
    toohigh = find(rmsvals > maxrms); %RMS noise floor too high
    toolow = find(rmsvals < minrms); %RMS noise floor too low

    badchans = [];
    if ~isempty(toohigh)
        badchans = [badchans; toohigh];
        fprintf('RMS noise of channel(s) %d is too high! Removed from analysis!\n',toohigh)
    end

    if ~isempty(toolow)
        badchans = [badchans; toolow];
        fprintf('RMS noise of channel(s) %d is too low! Removed from analysis!\n',toolow)
    end

    %Remove redundant channels
    badchans = unique(badchans); 
    igood = ~ismember(1:size(sig,1), badchans);

    fprintf('Done in %3.0fs!\n', toc);

end

