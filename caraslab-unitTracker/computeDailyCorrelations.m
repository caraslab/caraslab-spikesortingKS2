function correlations = computeDailyCorrelations(spiketimes)
fprintf('Computing daily pairwise cross-correlograms\n');
edges = 0:.100:1;

correlations = cell(size(spiketimes));
for iid=1:length(spiketimes)
    fprintf('Day %d \t',iid);
    
    % Each element of this cell array will be a matrix representing the
    % cross-correlograms between this neuron and every other neuron in the
    % population.  
    C = cell(1,length(spiketimes{iid}));
    for ii1=1:length(spiketimes{iid})
        % Rows ~ other neurons, cols ~ xcorr times
        C{ii1} = nan(length(spiketimes{iid}),length(edges)-1);
    end
    for ii1=1:length(spiketimes{iid})
        fprintf('.');
        for ii2=ii1+1:length(spiketimes{iid})
            h = relativeHist(spiketimes{iid}{ii1}-.5,spiketimes{iid}{ii2},edges);
            % Pre-normalize so that we can quickly compute the pearson 
            % correlation between the shapes of two cross-correlogrmas
            n = length(h);
            h = h-sum(h)/n;
            h = h/sqrt(sum(h.^2)/(n-1));
            h = reshape(h,1,numel(h));
            C{ii2}(ii1,:) = h;
            C{ii1}(ii2,:) = h;
        end
    end

    correlations{iid} = C;
    fprintf('\n');
end