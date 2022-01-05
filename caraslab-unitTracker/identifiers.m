function id = identifiers(survival)
% id = identifiers(survival)
%
% survival : survival matrix output from unitIdentification
% id : length #days cell array where id{i} is a length #neurons vector of 
%     integers which uniquely identify the same neurons across multiple
%     days.  So if id{i}(j) == id{m}(n), then day i unit j is the same
%     neuron as day m unit n.
id = cell(numel(survival)+1,1);
todayIds = 1:size(survival{1},1);
maxId = length(todayIds);
id{1} = todayIds';
for iid=2:length(id)
    yesterdayIds = todayIds;
    
    newSlots = ~sum(survival{iid-1}>0);
    newIds = maxId+(1:sum(newSlots));
    maxId = maxId+sum(newSlots);
    [r,c] = find(survival{iid-1});
    todayIds = nan(1,size(survival{iid-1},2));
    todayIds(c) = yesterdayIds(r);
    todayIds(newSlots) = newIds;
    
    id{iid} = todayIds';
end