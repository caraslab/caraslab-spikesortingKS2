function dickeyscore = computeDickeyScore(dickey)

% These are the normalization parameters required in Dickey's algorithm.
% They were calculated using a long series of recordings from a monkey
% named Frank.  Frank enjoys grapes, tang, and grooming his cagemate Tupac.
varNorm = [0.299147685915073;0.0635765890720025;0.0739105580696342;0.710622858394210;0.0122298060112537;0.125068458304676;4.21665762254236e-05;0.00236002865176362;];

dickeyscore = cell(numel(dickey)-1,1);
for iid=1:length(dickeyscore)
    dickeyscore{iid} = zeros(size(dickey{iid},2),size(dickey{iid+1},2));
    for unit1=1:size(dickey{iid},2)
        for unit2=1:size(dickey{iid+1},2)
            normed = (dickey{iid+1}(:,unit2)-dickey{iid}(:,unit1)).^2;
            normed = normed ./ varNorm;
            dickeyscore{iid}(unit1,unit2) = log(sqrt(sum(normed)));
        end
    end
    dickeyscore{iid}(~isfinite(dickeyscore{iid})) = min(dickeyscore{iid}(isfinite(dickeyscore{iid})));
end
