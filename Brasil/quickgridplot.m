function out = quickgridplot(in, timrange, freqrange)
% Usage: out = quickgridplot(in, timrange, freqrange)
% This plots frequency and position of fish in a grid recording
% Give this only one recording - e.g. cave(5)
% timrange (e.g. [0 120]) and freqrange (e.g. [200 500]) are optional.

% Set default time and frequency ranges if not provided by user
if nargin < 2
    for k=length(in.fish):-1:1
        maxtim(k) = in.fish(k).freq(end,1);
    end
    timrange = [0 max(maxtim)];
end
if nargin < 3
    for k=length(in.fish):-1:1
        maxfreq(k) = max(in.fish(k).freq(:,2));
        minfreq(k) = min(in.fish(k).freq(:,2));
    end
    freqrange = [min(minfreq)-50 max(maxfreq)+50];
end

% Plot
figure(1); clf; hold on; ylim(freqrange);
figure(2); clf; hold on;

for j=1:length(in.fish)
    
    tt = find(in.fish(j).freq(:,1) > timrange(1) & in.fish(j).freq(:,1) < timrange(2));
    
    % Frequency plot
    figure(1);
    plot(in.fish(j).freq(tt,1), in.fish(j).freq(tt,2), '.', 'MarkerSize', 8);
    
    % Position plot
    figure(2);
    plot(in.fish(j).x(tt), in.fish(j).y(tt), '.', 'MarkerSize', 8);
    
        
end

% grid plot

for j=1:length(in.fish)
    
    tt = find(in.fish(j).freq(:,1) > timrange(1) & in.fish(j).freq(:,1) < timrange(2));
        
end

out = 1;
