function [survival, score] = computeSurvival(shank, unit, corrscore, wavescore, autoscore, basescore, survival)

    sameChan = cell(size(corrscore));
    sameUnit = cell(size(corrscore));
    for iid=1:length(corrscore)
        sameChan{iid} = bsxfun(@eq,shank{iid}(:),shank{iid+1}(:)');
        sameUnit{iid} = sameChan{iid} & bsxfun(@eq,unit{iid}(:),unit{iid+1}(:)');
    end

    data = [unroll(corrscore) unroll(wavescore) unroll(autoscore) unroll(basescore)];
    good = ~isnan(sum(data));
    fprintf('Using %d classifiers, %d others are nan\n',sum(good),sum(~good));
    % Remove missing classifiers
    data = data(:,good);

    % MML edit: subtract min of each column 
    % Tried this strategy while troubleshooting; not to be used
    % datamin = min(data, [], 1); 
    % datamax = max(data, [], 1); 
    % 
    % denom = datamax - datamin;
    % denom(denom == 0) = 1; % Change zeros to 1 to avoid zero-division
    % 
    % % Scale to the max of each column 
    % data = (data - datamin) ./ denom;
    % corrscore = cellfun(@(x)(x-datamin(1)) ./ denom(1),corrscore,'UniformOutput',false); 
    % wavescore = cellfun(@(x)(x-datamin(2)) ./ denom(2),wavescore,'UniformOutput',false); 
    % autoscore = cellfun(@(x)(x-datamin(3)) ./ denom(3),autoscore,'UniformOutput',false); 
    % basescore = cellfun(@(x)(x-datamin(4)) ./ denom(4),basescore,'UniformOutput',false); 

    % Approximates expectation-maximization with partly specified labels
    % The E-step is making a hard assignment which is not exactly right but in
    % the context of this type of data where most labels are specified it makes
    % no difference
    unrollSameChan = unroll(sameChan);
    C = unroll(survival);
    for i=1:10
        [C,~,P] = classify(data, data, C,'quadratic');
        % Only same-channel is possible
        C = C & unrollSameChan;
    end

    % Identify a threshold for a 1% FP rate
    % MML comment: Do they mean 5% FP rate? See quantile() call below
    negative = data(~unroll(sameChan),:);
    [~,~,Pneg] = classify(negative,data,C,'quadratic');

    % MML edit: If classification fails just return
    if size(Pneg,2) == 1
    %     survival = reroll(corrscore, C);
        score = reroll(corrscore, C);
        return
    end

    negative = Pneg(:,2); 
    threshold = quantile(negative,.95);
    C = P(:,2)>threshold;
    C = C & unrollSameChan;

    survival = reroll(corrscore, C);
    score = reroll(corrscore, P(:,2));

    % Occasionally survival will indicate something impossible, like the same 
    % unit becomes two different units or vice versa.  We need to fix all those 
    % cases.
    for iid=1:length(survival)
        for iic=shank{iid}(:)'
            left = shank{iid}==iic;
            right = shank{iid+1}==iic;
            survival{iid}(left,right) = takeBest(score{iid}(left,right), threshold);
        end
    end
end

function survival = takeBest(similarity, thresh)
    % Posterior must be > thresh
    [pre,post] = find(similarity > thresh);
    pre = unique(pre); post = unique(post);

    similarity(similarity < thresh) = nan;
    sim = similarity(pre,post);

    survival = zeros(size(similarity));
    if isempty(pre) && isempty(post)
        return
    end

    survival(pre,post) = eye(numel(pre),numel(post));
    survival = survival > 0;
    if(numel(pre)>numel(post))
        P = perms(pre);
        score = nan(size(P,1),1);
        for i=1:size(P,1)
            score(i) = sum(sum(sim(survival(P(i,:),post)), 'omitnan'));
        end
        P = P(find(score==max(score),1),:);
        survival(pre,post) = survival(P,post);
    else
        P = perms(post);
        score = nan(size(P,1),1);
        for i=1:size(P,1)
            score(i) = sum(sum(sim(survival(pre,P(i,:))), 'omitnan'));
        end
        P = P(find(score==max(score),1),:);
        survival(pre,post) = survival(pre,P);
    end
end
