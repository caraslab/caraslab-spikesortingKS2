function score = computeCorrScore(correlations, survival)
score = cell(numel(correlations)-1,1);
for iid=1:length(score)
    score{iid} = zeros(length(correlations{iid}),length(correlations{iid+1}));
end

for iid=1:length(score)
    corr1 = correlations{iid};
    corr2 = correlations{iid+1};
    % Rearrange the pairwise cross correlations so they correspond to
    % putative same neurons across days
    [r,c] = find(survival{iid});
    if(isempty(r))
        continue;
    end
    corr1 = cellfun(@(X) X(r,:), corr1, 'UniformOutput', false);
    corr2 = cellfun(@(X) X(c,:), corr2, 'UniformOutput', false);
    
    for unit1=1:length(correlations{iid})
        for unit2=1:length(correlations{iid+1})
            left = corr1{unit1};
            right = corr2{unit2};
            n = size(left,2);
            r = atanh(sum(left.*right,2) ./ (n-1));
            score{iid}(unit1,unit2) = nanmean(r);
%             if(survival{iid}(unit1,unit2))
%                 for i=1:min(25,size(corr1{unit1},1))
%                     subplot(5,5,i);
%                     plot([corr1{unit1}(i,:); corr2{unit2}(i,:)]');
%                 end
%                 r2 = nan(size(r));
%                 for i=1:length(r)
%                     r2(i) = corr(corr1{unit1}(i,:)',corr2{unit2}(i,:)');
%                 end
%                 keyboard;
%             end
        end
    end
    if(any(~isfinite(score{iid}(:)))) % This means that only one unit survived
        score{iid}(~isfinite(score{iid})) = 0;
    end
end