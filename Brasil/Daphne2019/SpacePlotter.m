function SpacePlotter(in, casu, num)

figure(5); clf; 

subplot(221); surf(in(num).realhist); view(90,-90);
subplot(222); surf(in(num).randhist); view(90,-90)

if casu == 1 % cave
    xminedge = -150; xscale = 400; 
    yminedge = -150; yscale = 275;
end
if casu == 2 % surface
    xminedge = -150; xscale = 300;  
    yminedge = -100; yscale = 275;
end



subplot(223); plot(in(num).realXY(1,:), in(num).realXY(2,:), '.', 'MarkerSize', 8); 
    xlim([xminedge-20, xminedge+xscale+20]); ylim([yminedge-20, yminedge+yscale+20]);
subplot(224); plot(in(num).rndXY(1,:), in(num).rndXY(2,:), '.', 'MarkerSize', 8); 
    xlim([xminedge-20, xminedge+xscale+20]); ylim([yminedge-20, yminedge+yscale+20]);

    
    