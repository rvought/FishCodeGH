% Calculate CAVE dFs TIME CONSUMING - ONLY DO THIS ONCE prior to running this script!
% [caveDF, CalldFs] = dFanalysis(cave);
% [srfDF, SalldFs] = dFanalysis(srf);


% All dFs histogram
numbins = 20; 
maxfreq = 200;   
ctrs = 0:maxfreq/numbins:maxfreq;
plotspots = ctrs(2:end) - maxfreq/(2*numbins);

Asurfhist = histcounts(SalldFs, ctrs);
    totalsurf = sum(Asurfhist);
    Asurfhist = Asurfhist / totalsurf;
Acavehist = histcounts(CalldFs, ctrs);
    totalCave = sum(Acavehist);
    Acavehist = Acavehist / totalCave;

figure(1); clf; hold on;
    bar(plotspots, Asurfhist, 'b');
    bar(plotspots, -Acavehist, 'r');
    ylim([-0.2, 0.2]);

    
% Low dFs histogram
numbins = 10; 
maxfreq = 10;   
ctrs = 0:maxfreq/numbins:maxfreq;
plotspots = ctrs(2:end) - maxfreq/(2*numbins);

surfhist = histcounts(SalldFs, ctrs);
    surfhist = surfhist / totalsurf;
cavehist = histcounts(CalldFs, ctrs);
    cavehist = cavehist / totalCave;

figure(2); clf; hold on;
    bar(plotspots, surfhist, 'b');
    bar(plotspots, -cavehist, 'r');
    ylim([-0.03, 0.03]);

% Best low-dF example    
% Surface entry 5, fish 10 and 11 and pair 172 in the output of dFanalysis

figure(3); clf;
subplot(211); hold on;

    % Get valid data
    idxA = find(~isnan(srf(5).fish(10).freq(:,2)));
    idxB = find(~isnan(srf(5).fish(11).freq(:,2)));
    % Within the specified time range
    startim = 0; endtim = 300;
    idxA = idxA(find(srf(5).fish(10).freq(idxA,1)) > startim & srf(5).fish(10).freq(idxA,1) < endtim);
    idxB = idxB(find(srf(5).fish(11).freq(idxB,1)) > startim & srf(5).fish(11).freq(idxB,1) < endtim);
    
    % Plot
    yyaxis left;    
        plot(srf(5).fish(10).freq(idxA,1), srf(5).fish(10).freq(idxA,2), '.g', 'MarkerSize', 10);
        plot(srf(5).fish(11).freq(idxB,1), srf(5).fish(11).freq(idxB,2), '.c', 'MarkerSize', 10);
        ylim([360, 400]);
    yyaxis right;
        st = find(srfDF(5).pair(172).sharedtims > startim & srfDF(5).pair(172).sharedtims < endtim);
        plot(srfDF(5).pair(172).sharedtims(st), srfDF(5).pair(172).descartes(st), '.k', 'MarkerSize', 10)
        ylim([0, 300]);
        

ax(1) = subplot(223); hold on;
    plot(srf(5).fish(11).x(idxB), srf(5).fish(11).y(idxB), '.', 'MarkerSize', 10, 'Color', '[0.5, 0.5, 0.5]');
    plot(srf(5).fish(10).x(idxA), srf(5).fish(10).y(idxA), '.g', 'MarkerSize', 10);
    
ax(2) = subplot(224); hold on;
    plot(srf(5).fish(10).x(idxA), srf(5).fish(10).y(idxA), '.', 'MarkerSize', 10, 'Color', '[0.5, 0.5, 0.5]');
    plot(srf(5).fish(11).x(idxB), srf(5).fish(11).y(idxB), '.c', 'MarkerSize', 10);
    
linkaxes(ax, 'xy'); ylim([-100 125]); xlim([-125 100]);

clear Asurfhist surfhist Acavehist cavehist plotspots numbins maxfreq ctrs idxA idxB st startim endtim;
    