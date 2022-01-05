function base = computeDailyBaserates(spiketimes, gwfparamses)
base = cell(size(spiketimes));

for day=1:length(spiketimes)
    base{day} = nan(1,length(spiketimes{day}));
    for iic=1:length(spiketimes{day})
        max_time    = range(spiketimes{day}{iic});
        ops = gwfparamses{day}{iic}.ops;
        % Subtract periods of noise removal
        if isfield(ops, 'concat_tranges')
            % Read the concat_tranges to subtract the max_time from the offsets
            trange_or_tranges = ops.concat_tranges;
        else
            % for non-concatenated files, just read the offset
            trange_or_tranges = ops.trange;
        end

        for row=1:size(trange_or_tranges, 1)
            max_time = max_time - trange_or_tranges(row, 1);
        end

        base{day}(iic) = log(length(spiketimes{day}{iic})./ max_time);
        
        % old code
        % base{day}(iic) = log(length(spiketimes{day}{iic})./range(spiketimes{day}{iic}));
    end
end