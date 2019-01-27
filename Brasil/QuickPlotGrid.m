function QuickPlotGrid(in, freqs, tims)
% Usage: QuickPlotGrid(in, freqs, tims)
% where freqs and tims are optional. 
% Freqs in Hz, 

if nargin == 1 
    freqs = [200 600]; % Default frequency range for EODs
end    

freqs = sort(freqs);

   figure(1); clf; hold on;
   figure(2); clf; hold on;
   
for j = 1:length(in.fish)

        tt = find(~isnan(in.fish(j).freq(:,2))); % all valid data
        
        if nargin == 3 % use only data in the range specified by the user
            tims = sort(tims);
            uu = find(in.fish(j).freq(:,1) > tims(1) & in.fish(j).freq(:,1) < tims(2));
            tt = intersect(tt, uu);
        end
   
   
   figure(1); plot(in.fish(j).freq(tt,1), in.fish(j).freq(tt,2), '.', 'MarkerSize', 8); ylim(freqs);
   figure(2); plot(in.fish(j).x(tt), in.fish(j).y(tt), '.', 'MarkerSize', 8); 
        axis([-200, 150, -100, 250]);
    
end