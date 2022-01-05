function score = computeAutoScore(auto)
score = cell(numel(auto)-1,1);

for iid=1:length(score)
    if isempty(auto{iid})
        S = [];
    else
        S = corr(auto{iid},auto{iid+1});
    end
    % An empty autocorrelogram gets a correlation of 0
    S(isnan(S)) = 0;
    score{iid} = atanh(S);
end