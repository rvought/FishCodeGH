function out = MajorTom(in, tim, feesh)
% This generates stats on frequency and movement for each individual fish

if nargin < 3
    feesh = 1:length(in.fish);
end
if nargin < 2
    tim(2) = in.fish(1).freq(end,1);
    tim(1) = 0;
end

%% Mean frequency and variability

for j = length(feesh):-1:1
    
   idx = find(in.fish(j).freq(:,1) > tim(1) & in.fish(j).freq(:,1) < tim(2));
   out(j).meanfreq = mean(in.fish(j).freq(idx(~isnan(in.fish(j).freq(idx,2))),2)); 
   out(j).varfreq = var(in.fish(j).freq(idx(~isnan(in.fish(j).freq(idx,2))),2)); 
   out(j).numfreq = length(in.fish(j).freq(~isnan(in.fish(j).freq(idx,2)),2)); 
    
end

%% Movement distance and velocity