function QuickPlotGrid(in, win)

   figure(1); clf; hold on;
   figure(2); clf; hold on;
   
for j = 1:length(in.fish)
    
   figure(1); plot(in.t, in.fish(j).freq, '*'); 
   figure(2); plot(in.fish(j).x, in.fish.y, '*-');
    
end