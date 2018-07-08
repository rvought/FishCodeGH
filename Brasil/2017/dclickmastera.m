function out = dclickmastera(data, clks)
% Useage out = dclickmaster(data, clks)

numcolors = 20;
oclrs = hsv(numcolors); % 20 colors for plotting up to 20 different fishes
clrs = zeros(numcolors, 3);
    shuff = [11 2 7 13 18 1 16 12 10 17 3 8 14 19 20 4 16 5 15 9];
    for i=1:20; clrs(i,:) = oclrs(shuff(i),:); end; 

numclicks = length(clks);

% Plot the frequency data
figure(1); clf; hold on;

for i=1:length(data)
    plot(data(i).tim, data(i).freq, 'Color', clrs(i,:), 'LineWidth', 2);
end

% Plot Nicole's clicks

figure(1); hold on;

for i=1:length(clks)
    text(mean([clks(1).clicktimes(i*2) clks(1).clicktimes((i*2)-1)]), clks(i).meanfreq, num2str(i));
end

figure(2); clf; hold on;

    for k = 1:length(data)
        plot(data(k).x, data(k).y, 'Color', clrs(k,:), 'LineWidth', 1);
    end

    ax = axis;

out = [];

% Now things get rough - do the XY

EventNum = 1;

while EventNum ~= 99

    EventNum

    currfish = input('Enter fish number: ');

    startim = clks(1).clicktimes((EventNum*2)-1);
    endtim = clks(1).clicktimes(EventNum*2);
    
    repeater = 0;
        
% Loop up to the current click time

pad = 10;

for cc = startim-pad:0.2:endtim

    figure(4); clf;
 
        subplot(121); hold on; axis(ax);
     
     fulltt = find(data(currfish).tim > startim-pad & data(currfish).tim < endtim+pad);
     plot(data(currfish).tx(fulltt), data(currfish).ty(fulltt), 'y-', 'LineWidth', 10);
        
for k = 1:length(data)
     tt = find(data(k).tim > cc & data(k).tim < cc+pad);
     if k == currfish; wid = 5; 
     else
         wid = 1; 
     end
     if isempty(tt) == 0
     plot(data(k).tx(tt), data(k).ty(tt), 'Color', clrs(k,:), 'LineWidth', wid);
     plot(data(k).tx(tt(end)), data(k).ty(tt(end)), 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', clrs(k,:), 'MarkerSize', 10);
     end
end
       subplot(122); hold on; %xlim([clks(EventNum).tims(EventNum)-30 clks.tims(EventNum)+30]);

     plot(data(currfish).tim(fulltt), data(currfish).freq(fulltt), 'y-', 'LineWidth', 10);
       
for k = 1:length(data)
     tt = find(data(k).tim > cc & data(k).tim < cc+10);
     if k == currfish; wid = 5; else; wid = 1; end;
     if isempty(tt) == 0
     plot(data(k).tim(tt), data(k).freq(tt), 'Color', clrs(k,:), 'LineWidth', wid);
     plot(data(k).tim(tt(end)), data(k).freq(tt(end)), 'o', 'MarkerEdgeColor', clrs(k,:), 'MarkerFaceColor', clrs(k,:));
     end;
end
       
    pause(0.045);

end

    repeater = input('Enter next EventNum if OK, 00 to repeat, 99 to end. ');
    EventNum = repeater;
    

end % End loop for all events

