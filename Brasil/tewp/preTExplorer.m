function preTExplorer(dFdist, orig)


for j=1:length(dFdist)

    figure(1); clf;
    
    subplot(311); % plot original frequency data
    hold on;
    plot(orig(dFdist(j).fishnums(1)).tim, orig(dFdist(j).fishnums(1)).EOD, '.b', 'MarkerSize', 2);
    plot(orig(dFdist(j).fishnums(1)).tim, orig(dFdist(j).fishnums(2)).EOD, '.m', 'MarkerSize', 2);
    xlim([0 orig(dFdist(j).fishnums(1)).tim(end)]);
    ylim([250 500]);
    text(100, 450, num2str(dFdist(j).fishnums));
    text(100, 400, num2str(j));

    subplot(312); % plot Distance and dF data
    xlim([0 orig(dFdist(j).fishnums(1)).tim(end)]);
    hold on;
    yyaxis right;
    plot(orig(dFdist(j).fishnums(1)).tim, dFdist(j).distance);
    ylabel('distance')
    yyaxis left;
    plot(orig(dFdist(j).fishnums(1)).tim, dFdist(j).dF);
    ylabel('dF')
    
    subplot(325);
    plot(orig(dFdist(j).fishnums(1)).xy(:,1), orig(dFdist(j).fishnums(1)).xy(:,2), '.b', 'MarkerSize', 4);
    xlim([-250, 250]);
    ylim([-250, 250]);
    subplot(326);
    plot(orig(dFdist(j).fishnums(2)).xy(:,1), orig(dFdist(j).fishnums(2)).xy(:,2), '.m', 'MarkerSize', 4);
    xlim([-250, 250]);
    ylim([-250, 250]);
    
    a = input('foobar 9 to end: ');  
    if a == 9; break; end
    
end


    