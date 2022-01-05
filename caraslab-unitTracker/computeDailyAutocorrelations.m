function auto = computeDailyAutocorrelations(spiketimes)
auto = cell(size(spiketimes));

edges = 0:.005:.100;

for iid=1:length(spiketimes)
    auto{iid} = nan(length(edges)-1,length(spiketimes{iid}));
    for iic=1:length(spiketimes{iid})
        auto{iid}(:,iic) = relativeHist(spiketimes{iid}{iic},spiketimes{iid}{iic},edges);
    end
end