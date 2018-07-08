function scrambledeggs = shuffler(eggs)
% out = shuffler(spiketimes)
% where spikes are the spike times

scrambledeggs = zeros(1,length(eggs)); % initialize to speed up the script

% Get inter-spike-intervals

inter = diff(eggs); % Interspike intervals in seconds

% scramble the interspike intervals

trnie = inter(randperm(length(inter))); % trnie is scrambled inter

% Reconstruct a new spike train with these new intervals

firstspike = abs(randn(1) * 0.005); % Random time for the first spike in the range of about 0.001 to 0.025 seconds

scrambledeggs(1) = firstspike; % For convenience, take the first spike at 10 msec

for i=1:length(trnie); 
    scrambledeggs(i+1) = scrambledeggs(i) + trnie(i); 
end





    
    
    