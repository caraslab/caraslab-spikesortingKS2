function string = unroll(daily)
daily = cellfun(@(day) day(:), daily, 'UniformOutput', false);
string = cell2mat(daily(:));
end