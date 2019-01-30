function [out, cmbs] = SpaceCorps(in)
% Usage out = SpaceCorps(in)
% Find spatial relations between fish
% Distance histograms, range overlap

% Spatial centers for histogram
    ctrs{1} = -150:5:250;
    ctrs{2} = -150:5:250;

% Distance centers for histogram
    dctrs = 1:5:300;
    
numfish = length(in.fish); % How many fish in this recording

for j = numfish:-1:1 % For each fish
    % Get valid data (check frequency data)
    out(j).valididx = find(~isnan(in.fish(j).freq(:,2)));
    tmp(j).xy(1,:) = in.fish(j).x(out(j).valididx);
    tmp(j).xy(2,:) = in.fish(j).y(out(j).valididx);  

   % Spatial histogram
    out(j).fhist = hist3(tmp(j).xy', ctrs);
        
end

if numfish > 1 % We have more than one fish
    
    combos = combnk(1:numfish, 2); % All pairwise combinations of fish

    for p = length(combos):-1:1 % For each pair of fish

        cmbs(p).fishnums = combos(p,:); % Save the output combo
        
        % Get only entries for which we have data for both fish
        [~, idx1, idx2] = intersect(out(combos(p,1)).valididx, out(combos(p,2)).valididx);
       
        % Calculate the interfish-distance
        for jj = 1:length(idx1) % THERE HAS GOT TO BE A MORE EFFICIENT WAY OF DOING THIS...
            cmbs(p).dist(jj) = pdist2(tmp(combos(p,1)).xy(:,idx1(jj))', tmp(combos(p,2)).xy(:,idx2(jj))');
        end
        
        % Distance histogram for each pair of fish
        
            cmbs(p).dhist = hist(cmbs(p).dist, dctrs);
        
        % Compare spatial histograms
        
        

    end

end
