function [out, cmbs] = SpaceCorps(in)
% Usage out = SpaceCorps(in)
% Find spatial relations between fish
% Distance histograms, range overlap

    ctrs{1} = -150:5:250;
    ctrs{2} = -150:5:250;

numfish = length(in.fish); % How many fish in this recording

for j = numfish:-1:1 % For each fish
    tmp(j).valididx = find(~isnan(in.fish(j).freq(:,2)));
    tmp(j).xy(1,:) = in.fish(j).x(tmp(j).valididx);
    tmp(j).xy(2,:) = in.fish(j).y(tmp(j).valididx);  
    out(j).f1h = hist3(tmp(j).xy', ctrs);
end

if numfish > 1 % We have more than one fish
    
    combos = combnk(1:numfish, 2); % All pairwise combinations of fish

    for p = length(combos):-1:1 % For each pair of fish

        cmbs(p).fishnums = combos(p,:); % Save the output combo
        tt = intersect(tmp(combos(p,1)).valididx, tmp(combos(p,2)).valididx);
        cmbs(p).dist = pdist2(tmp(combos(p,1)).xy(tt,:), tmp(combos(p,2)).xy(tt,:));

    end

end
