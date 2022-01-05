function daily = reroll(daily, string)
caret = 0;
for iid=1:length(daily)
    daily{iid}(:) = string(caret+(1:numel(daily{iid})));
    caret = caret + numel(daily{iid});
end
end