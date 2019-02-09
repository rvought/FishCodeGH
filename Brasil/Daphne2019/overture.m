function out = overture(in, fishlist)
% Calculates the overlap in movements for all pairs of fish

if nargin == 2
   numfish = length(fishlist);
else
   numfish = length(in);
   fishlist = 1:numfish;
end
    combos = combnk(1:numfish, 2); % All pairwise combinations of fish

    for p = length(fishlist):-1:1 % For each pair of fish



diffhist = in(1).realhist - in(2).realhist;

overlapfishone = 1 - (sum(diffhist(diffhist > 0)) / sum(sum(in(1).realhist)));
overlapfishtwo = 1 - (sum(diffhist(diffhist < 0)) / sum(sum(in(2).realhist)));



    end
    
    

