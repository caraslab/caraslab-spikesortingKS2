function [survival, corrscore, wavescore, autoscore, basescore, score] = iterateSurvival(shank, unit, correlations, wavescore, autoscore, basescore)
% [survival, score] = iterateSurvival(channel, unit, correlations, wavescore, autoscore, basescore)
% Each of these arguments is a #days x 1 cell array of:
% channel: #neurons vector of channel numbers
% unit: #neurons vector of sort id numbers
% spiketimes: #neurons cell array of #spikes vector of spike times
% 
% wavescore is a #days-1 cell array where each entry is a #cells-today x 
%   #cells-tomorrow peak of cross-correlogram between wave shapes


survival = cell(numel(correlations)-1,1);
for iid=1:length(survival)
    survival{iid} = bsxfun(@eq,shank{iid}(:),shank{iid+1}(:)') & bsxfun(@eq,unit{iid}(:),unit{iid+1}(:)');
end

corrscore = computeCorrScore(correlations, survival);

for iterations=1:50
    [new_survival, score] = computeSurvival(shank, unit, corrscore, wavescore, autoscore, basescore, survival);
    
    change = unroll(survival)-unroll(new_survival);
    if(sum(abs(change))==0)
        break;
    else
        fprintf('%d identities changed\n',sum(abs(change(:))));
        survival = new_survival;
        corrscore = computeCorrScore(correlations, survival);
    end
end