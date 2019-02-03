function SpacePlotter(in, casu, num)

figure(5); clf; 

subplot(221); surf(in(num).realhist); view(90,-90);
subplot(222); surf(in(num).randhist); view(90,-90)


subplot(223); plot(in(num).realXY(1,:), in(num).realXY(2,:), '.', 'MarkerSize', 8); 
    xlim([in(1).xbounds(1)-20, in(1).xbounds(2)+20]); ylim([in(1).ybounds(1)-20, in(1).ybounds(2)+20]);
subplot(224); plot(in(num).rndXY(1,:), in(num).rndXY(2,:), '.', 'MarkerSize', 8); 
    xlim([in(1).xbounds(1)-20, in(1).xbounds(2)+20]); ylim([in(1).ybounds(1)-20, in(1).ybounds(2)+20]);

    
    