function out = quickgridplot(in, timrange, freqrange)
% Usage: out = quickgridplot(in, timrange, freqrange)
% This plots frequency and position of fish in a grid recording
% Give this only one recording - e.g. cave(5)

if nargin < 2
    for k=length(in.fish):-1:1
        maxtim(k) = in.fish(k).fish(end,1);
    end
    timrange = [0 max(maxtim)];
end
if nargin < 3
    for k=length(in.fish):-1:1
        maxfreq(k) = max(in.fish(k).fish(:,2));
        minfreq(k) = min(in.fish(k).fish(:,2));
    end
    freqrange = [min(minfreq)-50 max(maxfreq)+50];
end

% Frequency plot
figure(1); clf; hold on; ylim(freqrange);

for j=1:length(in.fish)
    
    tt = find(in.fish(j).freq(:,1) > timrange(1) & in.fish(j).freq(:,1) < timrange(2));
    plot(in.fish(j).freq(tt,1), in.fish(j).freq(tt,2), '.', 'MarkerSize', 8);
    
    
end

out = 1;
