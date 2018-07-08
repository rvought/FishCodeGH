function out = clickmaster(data, clks)
% Useage out = clickmaster(data, clks)
% This is for playback of signals

numcolors = 50;
    oclrs = hsv(numcolors); % 50 colors for plotting up to 50 different fishes
    clrs = zeros(numcolors, 3);
    shuff = randperm(numcolors);
    %shuff = [11 2 7 13 18 1 16 12 10 17 3 8 14 19 20 4 16 5 15 9];
    for i=1:numcolors; clrs(i,:) = oclrs(shuff(i),:); end; 

numclicks = length(clks.freqs);

% Plot the frequency data
figure(1); clf; hold on;
for i=1:length(data)
    plot(data(i).tim, data(i).freq, '*', 'Color', clrs(i,:));
    text(data(i).tim(1)+0.1, mean(data(i).freq)+2, num2str(i), 'Color', clrs(i,:), 'FontSize', 18);
end

% Plot Nicole's clicks

figure(1); hold on;
    for i=1:length(clks.tims)
        text(clks.tims(i), clks.freqs(i), num2str(i), 'FontSize', 16);
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

    EventNum = input('Which event? (Enter 99 to quit) ');

    currfish = input('Enter fish number to highlight: ');

    currtim = clks.tims(EventNum);
    repeater = 0;
    while repeater < 20
    
        startim = currtim - 20;
        startdur = 10;

% Loop up to the current click time
for cc = 5:0.2:40

    endtim = startim + cc;
    figure(4); clf;
 
        subplot(121); hold on; axis(ax);
     
     fulltt = find(data(currfish).tim > currtim-20 & data(currfish).tim < currtim+20);
     plot(data(currfish).tx(fulltt), data(currfish).ty(fulltt), 'y-', 'LineWidth', 10);
        
for k = 1:length(data)
     tt = find(data(k).tim > endtim-startdur & data(k).tim < endtim);
     if k == currfish; wid = 5; else; wid = 2; end;
     if isempty(tt) == 0
     plot(data(k).tx(tt), data(k).ty(tt), 'Color', clrs(k,:), 'LineWidth', wid);
     plot(data(k).tx(tt(end)), data(k).ty(tt(end)), 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', clrs(k,:));
     end;
end
       subplot(122); hold on; xlim([clks.tims(EventNum)-30 clks.tims(EventNum)+30]);

     plot(data(currfish).tim(fulltt), data(currfish).freq(fulltt), 'y-', 'LineWidth', 10);
       
for k = 1:length(data)
     tt = find(data(k).tim > startim & data(k).tim < endtim);
     if k == currfish; wid = 5; else; wid = 2; end;
     if isempty(tt) == 0
     plot(data(k).tim(tt), data(k).freq(tt), 'Color', clrs(k,:), 'LineWidth', wid);
     plot(data(k).tim(tt(end)), data(k).freq(tt(end)), 'o', 'MarkerEdgeColor', clrs(k,:), 'MarkerFaceColor', clrs(k,:));
     end;
end
       
    pause(0.05);
end

    repeater = input('Enter 99 if OK, 00 to repeat. ');
    
    end

end % End loop for all events

