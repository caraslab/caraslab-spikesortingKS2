function score = computeBaseScore(base)
score = cell(numel(base)-1,1);

for iid=1:length(score)
    score{iid} = bsxfun(@minus,base{iid+1}(:)',base{iid}(:));
end