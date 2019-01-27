function QuickPlotGrid(in)

   figure(1); clf; hold on;
   figure(2); clf; hold on;
   
for j = 1:length(in.fish)
   
   tt = find(~isnan(in.fish(j).freq(:,2)));
   figure(1); plot(in.fish(j).freq(:,1), in.fish(j).freq(:,2), '*'); ylim([200 600]);
   figure(2); plot(in.fish(j).x(tt), in.fish(j).y(tt), '*');
    
end