function out = SpaceCorps(in)
% Usage out = SpaceCorps(in)
% Find spatial relations between fish
% Distance histograms, range overlap

    ctrs{1} = -150:5:250;
    ctrs{2} = -150:5:250;


combos = combnk(1:numfish, 2); % All pairwise combinations of fish

for p = length(combos):-1:1 % For each pair of fish

    out(p).fishnums = combos(p,:); % Save the output combo

    % Need to get only data points for which we have a valid frequency reading
    ttf1 = find(~isnan(in.fish(combos(p,1)).freq(:,2))); % First fish of combo
    ttf2 = find(~isnan(in.fish(combos(p,2)).freq(:,2))); % Second fish of combo

    F1(1,:) = in.fish(combos(p,1)).x(ttf1);
    F1(2,:) = in.fish(combos(p,1)).y(ttf1);

    F2(1,:) = in.fish(combos(p,2)).x(ttf2);
    F2(2,:) = in.fish(combos(p,2)).y(ttf2);

    out(p).f1h = hist3(F1', ctrs);
    out(p).f2h = hist3(F2', ctrs);
    
end

