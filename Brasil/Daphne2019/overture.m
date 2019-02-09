function out = overture(in, fishlist)
% Calculates the overlap in movements for all pairs of fish

if nargin == 2
   numfish = length(fishlist);
else
   numfish = length(in);
   fishlist = 1:numfish;
end

    combos = combnk(fishlist, 2); % All pairwise combinations of fish

    for p = numfish:-1:1 % For each fish

        for n = length(combos) % For each pair of fish

            if ~isempty(combos(:,n) == p)
                
            out(p).combo(n).diffhist = in(fishlist(combos(1,n))).realhist - in(fishlist(combos(1,n))).realhist;
            
            if fishlist(combos(1,n) == p
                overlapfishone = 1 - (sum(diffhist(diffhist > 0)) / sum(sum(in(1).realhist)));
            else
            overlapfishtwo = 1 - (sum(diffhist(diffhist < 0)) / sum(sum(in(2).realhist)));
            end
            
            end
        end
    end
    
    

