function out = overture(in, fishlist)
% Calculates the overlap in movements for all pairs of fish

if nargin == 2
   numfish = length(fishlist);
else
   numfish = length(in);
   fishlist = 1:numfish;
end

out(1).overlap = []; out(1).overfish = [];

    combos = combnk(fishlist, 2)
    % All pairwise combinations of fish

    for p = numfish:-1:1 % For each fish

        for n = length(combos) % For each pair of fish

            if ~isempty(find([combos(:,n)] == p, 1))
                
            out(p).combo(n).diffhist = in(fishlist(combos(1,n))).realhist - in(fishlist(combos(2,n))).realhist;
            
            if fishlist(combos(1,n)) == p
                out(p).overlap(end+1) = 1 - (sum(out(p).combo(n).diffhist(out(p).combo(n).diffhist > 0)) / sum(sum(in(fishlist(combos(1,n))).realhist)));
                out(p).overfish(end+1) = fishlist(combos(2,n));
            else
                out(p).overlap(end+1) = 1 - abs((sum(out(p).combo(n).diffhist(out(p).combo(n).diffhist < 0)) / sum(sum(in(fishlist(combos(2,n))).realhist))));
                out(p).overfish(end+1) = fishlist(combos(1,n));
            end
            
            end
        end
    end
    
    

