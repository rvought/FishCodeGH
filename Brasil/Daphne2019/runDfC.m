% This script produces figure 4 for the manuscript - pairwise dF interactions
% between fish in cave and surface habitats. It does not look at the
% correlation between distance and dF. Relies on dFanalysis.m

% Calculate CAVE dFs TIME CONSUMING - ONLY DO THIS ONCE prior to running this script!
% load SurfaceDataRev2018a.mat; load CaveDataRev2018a.mat;
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


fprintf('Surface mean dF %2.4f var dF %2.4f and n= %i. \n', mean(SalldFs), std(SalldFs), length(SalldFs));
fprintf('Cave mean dF %2.4f var dF %2.4f and n= %i. \n', mean(CalldFs), std(CalldFs), length(CalldFs));
[aa,bb,cc,dd] = ttest2(SalldFs, CalldFs);
fprintf('T-test surface versus cave dfs: P=%2.4f, tstat=%4.2f and df=%i. \n', bb, dd.tstat, dd.df);

%% Get low dF stats
caveLessthanTen = 0; caveLessthanFive = 0;
surfLessthanTen = 0; surfLessthanFive = 0;

for j=1:length(caveDF)
    if ~isempty(caveDF(j).pair)
    for k=1:length(caveDF(j).pair)
        if caveDF(j).pair(k).dFmean < 10
            caveLessthanTen = caveLessthanTen + 1;
        if caveDF(j).pair(k).dFmean < 5
            caveLessthanFive = caveLessthanFive + 1;
        end
        end
    end
    end
end
for j=1:length(srfDF)
    if ~isempty(srfDF(j).pair)
    for k=1:length(srfDF(j).pair)
        if srfDF(j).pair(k).dFmean < 10
            surfLessthanTen = surfLessthanTen + 1;
        if srfDF(j).pair(k).dFmean < 5
            surfLessthanFive = surfLessthanFive + 1;
        end
        end
    end
    end
end

%% Plot
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
    