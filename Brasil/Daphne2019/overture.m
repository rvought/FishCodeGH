function [out, plt] = overture(in, fishlist)
% Calculates the overlap in movements for all pairs of fish

if nargin == 2
   numfish = length(fishlist);
else
   numfish = length(in);
   fishlist = 1:numfish;
end

out(1).overlap = []; out(1).overfishnums = []; plt.alloverlaps = [];

combos = combnk(fishlist, 2); % All pairwise combinations of fish
    
    for p = numfish:-1:1 % For each fish
        
        % Do the self comparisons - Self-Jiggled and Self-Randomized
        
        out(p).selfjighist = in(fishlist(p)).realhist - in(fishlist(p)).jighist;
        out(p).selfjig = 1 - (sum(out(p).selfjighist(out(p).selfjighist > 0)) / sum(sum(in(fishlist(p)).realhist)));
        
        out(p).selfrndhist = in(fishlist(p)).realhist - in(fishlist(p)).randhist;
        out(p).selfrnd = 1 - (sum(out(p).selfrndhist(out(p).selfjighist > 0)) / sum(sum(in(fishlist(p)).realhist)));
              
        
        for n = 1:length(combos) % For each pair of fish
            
            if ~isempty(find(combos(n,:) == p, 1)) % Use only cases with focal fish in it (p)

                % Need a check to ensure input is not empty and that
                % matrices have same dimensions

                out(p).combo(n).diffhist = in(fishlist(combos(n,1))).realhist - in(fishlist(combos(n,2))).realhist;
                out(p).combo(n).diffhist = in(fishlist(combos(n,1))).jighist - in(fishlist(combos(n,2))).jighist;
                out(p).combo(n).diffhist = in(fishlist(combos(n,1))).realhist - in(fishlist(combos(n,2))).realhist;
            
            if fishlist(combos(n,1)) == p
                out(p).overfishnums(end+1) = fishlist(combos(n,2));
                out(p).realoverlap(end+1) = 1 - (sum(out(p).combo(n).diffhist(out(p).combo(n).diffhist > 0)) / sum(sum(in(fishlist(combos(n,1))).realhist)));
                out(p).jigoverlap(end+1) = 1 - (sum(out(p).combo(n).diffhist(out(p).combo(n).diffhist > 0)) / sum(sum(in(fishlist(combos(n,1))).realhist)));
                out(p).rndoverlap(end+1) = 1 - (sum(out(p).combo(n).diffhist(out(p).combo(n).diffhist > 0)) / sum(sum(in(fishlist(combos(n,1))).realhist)));
            else
                out(p).overlap(end+1) = 1 - abs((sum(out(p).combo(n).diffhist(out(p).combo(n).diffhist < 0)) / sum(sum(in(fishlist(combos(n,2))).realhist))));
                out(p).overfishnums(end+1) = fishlist(combos(n,1));
            end
            
            end
        end
        
        plt.alloverlaps = [plt.alloverlaps, out(p).overlap];        
        plt.alloverlaps = [plt.alloverlaps, out(p).overlap];        
        
    end % For each fish
    
    
%% Plot

lms = 0:0.1:1;

figure(5); clf; histogram(plt.alloverlaps, 'BinEdges', lms);

figure(6); clf;
subplot(211); histogram([out.selfjig], 'BinEdges', lms)
subplot(212); histogram([out.selfrnd], 'BinEdges', lms)
