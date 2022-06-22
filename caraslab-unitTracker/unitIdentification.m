function [survival, score, corrscore, wavescore, autoscore, basescore, correlations, auto] = ...
    unitIdentification(shank, unit, spiketimes, wmean, gwfparamses, shank_channelses, channel_map, varargin)

    %% MML notes: 
    % This function is an
    % implementation of their algorithm to Neuronexus probe data, i.e.
    % potential-same units are located on the same shank
    % Modifications to this algorithm:
    %   1. Waveforms on the same probe shank are treated as potential-same.
    %   2. Recording duration adapts to noise removal information at the
    %       beginning of recordings from the TDT recording system. Important in
    %       computeDailyBaserates.m
    %   3. Probe geometry is a key factor in computeWaveScore. Waveforms are
    %       compared between and across all channels (e.g. 16x16 comparisons in a 16ch shank)
    %       In these, I implemented a weighted Euclidean distance to confer more
    %       weight to comparisons between closer channels and maximum weight to
    %       same-channel comparisons

    %% Original docstring by Fraser, 2011
    % THE SPIKE TIMES MUST BE SORTED IN ASCENDING ORDER!!!  Failure to do this
    % will cause unpredictable results.
    %
    % You will need to compile relativeHist.c
    %
    % Each argument is a #days x 1 cell array where each element represents one
    % day of recording.
    % 
    % channel{day}: #neurons vector of channel numbers
    % unit{day}: #neurons vector of sort id numbers, would be between 1 and
    %   4 if your system has the capacity to sort 4 neurons on each channel.
    % spiketimes{day}: #neurons cell array of #spikes vector of spike
    %   times.  spiketimes{day}{i} would be the vector of spiketimes for
    %   channel{day}{i} unit number unit{day}{i}.
    % wmean{day}: #neurons cell array of #samples vector representing a mean
    %   waveform.  wmean{day}{i} would be the mean waveform for channel{day}{i}
    %   unit number unit{day}{i}.
    %
    % [survival] indicates where the same neuron survives from day to
    %   day.  survival{d}(i,j) == 1 if the i'th unit on day d represents the same
    %   neuron as the j'th unit on day d+1.
    %
    % [score] is the posterior probability of being the same neuron according
    %   to the Gaussian mixture model-based classifier.
    %
    % [corrscore] is the value for pairwise correlation score for identity
    %   computed in the last iteration, which you will need if you want to make
    %   plots of the parameters used for classification.
    %
    % [wavescore, autoscore, basescore, correlations] can be fed directly to
    %   iterateSurvival to save time if for some reason you need to run the
    %   algorithm again with exactly the same parameters.
    %
    % Options:
    % 'plot': call gaussianMix at the end, which makes a fancy plot of some of
    %   the similarity scores
    % 'dickey': Compute autocorrelation similarity using the method of Dickey, 
    %   A. S., A. Suminski, et al. (2009). J Neurophysiol 102(2).
    % 
    % Warnings: 
    %
    % It is important that the spiketimes vectors represent contiguous
    % recordings.  If there are gaps in your data you need to excise them so
    % that spiketimes represents a pseudo-continuous recording.  This is
    % necessary because the mean firing rate of each neuron is calculated as
    % length(spiketimes{day}{i})/range(spiketimes{day}{i}).  
    %
    % If you give this algorithm units that are poorly discriminated, it will
    % make a best effort but the results will be somewhat ill-defined due to
    % the fact that some sorted units may be combinations of multiple neurons
    % with the combination changing from day to day.  Having said that, the
    % algorithm is fairly robust to this assumption so you should not feel that
    % everything needs to be perfectly isolated.  A good criterion is that if
    % you feel confident a sorted unit is >90% from one neuron then it should
    % be included.
    %
    % This function can be very memory-intensive if you feed it a lot of data
    % at once.  For every gigabyte of memory you have you should be able to do
    % at least 10 sessions with 100 neurons/session.
    %
    % IT IS VERY IMPORTANT THAT YOU MINIMIZE THE PRESENCE OF MOTION ARTIFACT
    % AND ELECTRICAL CROSS-TALK IN YOUR DATA.  If the same unit is being
    % recorded on multiple channels, or if there is heavy motion artifact
    % causing simultaneous spikes on many channels, this will screw up the
    % cross-correlation based classifier.  Occasional motion artifact is
    % acceptable as long as it represents a negligable portion of the overall
    % spikes.  If you suspect that motion artifact is causing problems in your
    % data you should turn on the 'plot' option and look at the distribution of
    % pairwise scores.  If the blue distribution is significantly away from zero
    % for the pairwise scores then there may be something funny going on.
    %
    % George Fraser
    % fraser.george.w@gmail.com
    % 1/11/2011
    %%
    if(ismember('dickey',varargin))
        dickey = computeDailyDickey(spiketimes);
        autoscore = computeDickeyScore(dickey);
    else
        auto = computeDailyAutocorrelations(spiketimes);
        autoscore = computeAutoScore(auto);
    end
    base = computeDailyBaserates(spiketimes, gwfparamses);
    basescore = computeBaseScore(base);
    correlations = computeDailyCorrelations(spiketimes);
    wavescore = computeWaveScore(wmean, channel_map, shank_channelses);

    [survival, corrscore, wavescore, autoscore, basescore, score] = iterateSurvival(shank, unit, correlations, wavescore, autoscore, basescore);

    if(ismember('plot',varargin))
        gaussianMix(shank, unit, corrscore, wavescore, autoscore, basescore, survival);
    end

end