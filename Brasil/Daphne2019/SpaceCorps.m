function out = SpaceCorps(in)
% Usage out = SpaceCorps(in)
% Find spatial relations between fish
% Distance histograms, range overlap

    ctrs{1} = -150:5:250;
    ctrs{2} = -150:5:250;

numfish = length(in.fish); % How many fish in this recording

for j = numfish:-1:1 % For each pair of fish
    F1(1,:) = in.fish(j).x(~isnan(in.fish(j).freq(:,2)));
    out(j).f1h = hist3(F1', ctrs);
    clear F1
end

% if numfish > 1 % We have more than one fish
%     
% combos = combnk(1:numfish, 2); % All pairwise combinations of fish
% 
% for p = length(combos):-1:1 % For each pair of fish
% 
%     out(p).fishnums = combos(p,:); % Save the output combo
% 
% end
% 
% end
