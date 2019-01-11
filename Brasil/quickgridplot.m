function out = quickgridplot(in, freqrange, timrange)
% Usage: out = quickgridplot(in, freqrange, timrange)
% This plots frequency and position of fish in a grid recording



% Frequency plot
figure(1); clf; hold on; ylim(freqrange);

for j=1:length(in.fish)
    
    tt = find(in.fish(j).freq(:,1) > timrange(1) & in.fish(j).freq(:,1) < timrange(2));

    plot(in.fish(j).freq(tt,1), in.fish(j).freq(tt,2), '.', 'MarkerSize', 2);
    
    
end

out = 1;
