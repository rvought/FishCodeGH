function out = overture(in, fishlist)
% Calculates the overlap in movements for all pairs of fish

if nargin == 2
   numfish = length(fishlist);
else
   numfish = length(in);
   fishlist = 1:numfish;
end

out(1).overlap = []; out(1).overfish = [];

combos = combnk(fishlist, 2); % All pairwise combinations of fish
    
    for p = numfish:-1:1 % For each fish

        for n = 1:length(combos) % For each pair of fish
            
            if ~isempty(find(combos(n,:) == p, 1)) % Use only cases with focal fish in it (p)

                % Need a check to ensure input is not empty and that
                % matrices have same dimensions

                out(p).combo(n).diffhist = in(fishlist(combos(n,1))).realhist - in(fishlist(combos(n,2))).realhist;
            
            if fishlist(combos(n,1)) == p
                out(p).overlap(end+1) = 1 - (sum(out(p).combo(n).diffhist(out(p).combo(n).diffhist > 0)) / sum(sum(in(fishlist(combos(n,1))).realhist)));
                out(p).overfish(end+1) = fishlist(combos(n,2));
            else
                out(p).overlap(end+1) = 1 - abs((sum(out(p).combo(n).diffhist(out(p).combo(n).diffhist < 0)) / sum(sum(in(fishlist(combos(n,2))).realhist))));
                out(p).overfish(end+1) = fishlist(combos(n,1));
            end
            
            end
        end
    end
    
    

