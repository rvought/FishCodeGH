%% STARTER

clear caveDistance cavedF cavedFDistCC cavehist cavemindists clrs ctrs dFdist distmper
clear figA figB fishvar j jj kk numcolors orig out p pVal plotspots r 
clear surfDistance surfacemindists surfdFDistCC surfhist tt xlims


WindowDuration = 300;
OverlapDuration = 150;

%% CAVE

cavedFDistCC = [];
caveDistance = [];
cavedF = [];
cavemindists = [];

figure(1); clf;

for jj = [3 4 5 6 7 8 10 11 12 14]
    
    % Calculate all dFs and distances between pairs of fish
    [dFdist, orig] = preTE(cave(jj));

    % Extract minimum distances 
    distmper{dFdist(end).fishnums(2)} = []; 
    for kk = 1:length(dFdist)
        if ~isempty(dFdist(kk).distance)
            distmper{dFdist(kk).fishnums(1)}(end+1) = nanmin(dFdist(kk).distance);
            distmper{dFdist(kk).fishnums(2)}(end+1) = nanmin(dFdist(kk).distance);
        end
    end
    for kk = 1:length(distmper)
        cavemindists(end+1) = min(distmper{kk});        
    end
        clear distmp;
    
    % Calculate correlations between distance and dF for each pairwise
    % interaction
    out = getCorr(dFdist, orig, WindowDuration, OverlapDuration, [1 jj]);
    
    for j=1:length(out) 
        cavedFDistCC = [cavedFDistCC out(j).dfdistCC];
        caveDistance = [caveDistance out(j).meanDist];
        cavedF = [cavedF out(j).meandF];
    end

end

    figure(1);
    subplot(221); plot(caveDistance, cavedFDistCC, 'b.', 'MarkerSize', 4); 
    subplot(222); plot(cavedF, cavedFDistCC, 'm.', 'MarkerSize', 4); 
    
    subplot(223); plot(caveDistance, abs(cavedFDistCC), 'k.', 'MarkerSize', 4); 
    subplot(224); plot(cavedF, abs(cavedFDistCC), 'r.', 'MarkerSize', 4); 



% PEARSON CORRELATION COEFICIENT, Effect of mean distance on correlation between dF and distance    
% Is there a correlation between mean distance and the correlation between
% dF and distance?
[r, pVal] = corrcoef(caveDistance, cavedFDistCC); 
r = r(2);
pVal = pVal(2);
fprintf('Cave distance CC %2.4f, %2.4f \n', r, pVal)
ab = length(find(cavedFDistCC > 0));
bel = length(find(cavedFDistCC < 0));
fprintf('Cave above %i, Cave below %i, mean %2.4f \n', ab, bel, mean(cavedFDistCC))


% PEARSON CORRELATION COEFICIENT, Effect of mean dF on STRENGTH of correlation between dF and distance    
% Is there a correlation between mean distance and the correlation between
% dF and distance?
[r, pVal] = corrcoef(cavedF, cavedFDistCC); 
r = r(2);
pVal = pVal(2);
fprintf('Cave dF CC %2.4f, %2.8f \n', r, pVal)


% PEARSON CORRELATION COEFICIENT, Effect of mean dF on SIGN of correlation between dF and distance    
% Is there a correlation between mean distance and the correlation between
% dF and distance?
[r, pVal] = corrcoef(cavedF, abs(cavedFDistCC)); 
r = r(2);
pVal = pVal(2);
fprintf('Cave dF absCC %2.4f, %2.4f \n', r, pVal)


%% SURFACE

surfdFDistCC = [];
surfDistance = [];
surfdF = [];
surfacemindists = [];

figure(2); clf;

for jj = [1 2 3 4 5]
    
    % Calculate all dFs and distances between pairs of fish
    [dFdist, orig] = preTE(srf(jj));

    % Extract minimum distances 
    distmper{dFdist(end).fishnums(2)} = []; 
    for kk = 1:length(dFdist)
        if ~isempty(dFdist(kk).distance)
            distmper{dFdist(kk).fishnums(1)}(end+1) = nanmin(dFdist(kk).distance);
            distmper{dFdist(kk).fishnums(2)}(end+1) = nanmin(dFdist(kk).distance);
        end
    end
    for kk = 1:length(distmper)
        if ~isempty(distmper{kk})
        surfacemindists(end+1) = min(distmper{kk});        
        end
    end
        clear distmp;

    % Calculate correlations between distance and dF
    out = getCorr(dFdist, orig, WindowDuration, OverlapDuration, [2 jj]);    
    
    for j=1:length(out) 
        surfdFDistCC = [surfdFDistCC out(j).dfdistCC];
        surfDistance = [surfDistance out(j).meanDist];
        surfdF = [surfdF out(j).meandF];
    end

end

    figure(2);
    subplot(221); plot(surfDistance, surfdFDistCC, 'b.', 'MarkerSize', 4); 
    subplot(222); plot(surfdF, surfdFDistCC, 'm.', 'MarkerSize', 4); 
    
    subplot(223); plot(surfDistance, abs(surfdFDistCC), 'k.', 'MarkerSize', 4); 
    subplot(224); plot(surfdF, abs(surfdFDistCC), 'r.', 'MarkerSize', 4); 

% PEARSON CORRELATION COEFICIENT, Effect of mean distance on correlation between dF and distance    
% Is there a correlation between mean distance and the correlation between
% dF and distance?
[r, pVal] = corrcoef(surfDistance, surfdFDistCC); 
r = r(2);
pVal = pVal(2);
fprintf('Surface distance CC %2.4f, %2.4f \n', r, pVal)
ab = length(find(surfdFDistCC > 0));
bel = length(find(surfdFDistCC < 0));
fprintf('Surfae above %i, Surface below %i, mean %2.4f \n', ab, bel, mean(surfdFDistCC))



% PEARSON CORRELATION COEFICIENT, Effect of mean dF on STRENGTH of correlation between dF and distance    
% Is there a correlation between mean distance and the correlation between
% dF and distance?
[r, pVal] = corrcoef(surfdF, surfdFDistCC); 
r = r(2);
pVal = pVal(2);
fprintf('Surface dF CC %2.4f, %2.4f \n', r, pVal)


% PEARSON CORRELATION COEFICIENT, Effect of mean dF on SIGN of correlation between dF and distance    
% Is there a correlation between mean distance and the correlation between
% dF and distance?
% [r, pVal] = corrcoef(surfdF, surfdFDistCC); 
% r = r(2)
% pVal = pVal(2)


%% Plot histograms

% Distance
figure(4); clf;
xlims = [0 20 40 60 80 100 120 140 160 180 200 220 240];
subplot(221); histogram(caveDistance, xlims, 'Normalization', 'probability'); axis([0, 250, 0, 0.3]); 
subplot(223); histogram(surfDistance, xlims, 'Normalization', 'probability'); axis([0, 250, 0, 0.3]);

xlims = [0 2 4 6 8 10 12 14 16 18 20 22 24];
subplot(222); histogram(cavemindists, xlims, 'Normalization', 'probability'); axis([0, 25, 0, 0.3]); 
subplot(224); histogram(surfacemindists, xlims, 'Normalization', 'probability'); axis([0, 25, 0, 0.3]); 


%% Plot examples

[dFdist, ~] = preTE(cave(3));
figure(3); subplot(212); 
% yyaxis left; plot(dFdist(19).tim, dFdist(19).distance, '.', 'MarkerSize', 8); ylim([50 250]);
% yyaxis right; plot(dFdist(19).tim, dFdist(19).dF, '.', 'MarkerSize', 8); ylim([12 17]);
yyaxis left; plot(dFdist(19).tim(dFdist(19).tim ~= 0), dFdist(19).distance(dFdist(19).tim ~= 0), '-'); ylim([50 250]);
yyaxis right; plot(dFdist(19).tim(dFdist(19).tim ~= 0), dFdist(19).dF(dFdist(19).tim ~= 0), '-'); ylim([12 17]);
xlim([575 875]);

[dFdist, ~] = preTE(cave(4));
figure(3); subplot(211); 
% yyaxis left; plot(dFdist(18).tim, dFdist(18).distance, '.', 'MarkerSize', 8); ylim([0 200]);
% yyaxis right; plot(dFdist(18).tim, dFdist(18).dF, '.', 'MarkerSize', 8); ylim([84 88]);
yyaxis left; plot(dFdist(18).tim(dFdist(18).tim ~=0), dFdist(18).distance(dFdist(18).tim ~=0), '-'); ylim([0 200]);
yyaxis right; plot(dFdist(18).tim(dFdist(18).tim ~=0), dFdist(18).dF(dFdist(18).tim ~=0), '-'); ylim([84 88]);
xlim([300 600]);

%% Plot dF histograms

figure(5); clf; hold on;
    figA = figure(5);
    figA.Renderer = 'Painters';
    ctrs = 0:5:250;
    surfhist = histcounts(surfdF, ctrs);
        surfhist = surfhist / sum(surfhist);    
    cavehist = histcounts(cavedF, ctrs);
        cavehist = cavehist / sum(cavehist);
        
    plotspots = ctrs(2:end) - (5/2);

    bar(plotspots, surfhist, 'b');
    bar(plotspots, -cavehist, 'r');
    ylim([-0.1, 0.1]);
    
figure(6); clf;
    figB = figure(6);
    figB.Renderer = 'Painters';
subplot(211); histogram(cavedF, 0:1:20)
subplot(212); histogram(surfdF, 0:1:20)

[dFdist, orig] = preTE(srf(5));
tt{1} = srf(5).fish(10).freq(:,1) > 170 & srf(5).fish(10).freq(:,1) < 470;
tt{2} = srf(5).fish(11).freq(:,1) > 170 & srf(5).fish(11).freq(:,1) < 470;

figure(7); clf;
subplot(211); hold on;
plot(srf(5).fish(10).freq(tt{1},1)-170, srf(5).fish(10).freq(tt{1},2), '.');
plot(srf(5).fish(11).freq(tt{2},1)-170, srf(5).fish(11).freq(tt{2},2), '.');
ylim([380 395]);
subplot(223); hold on;
plot(srf(5).fish(10).x(tt{1}), srf(5).fish(10).y(tt{1}), '.');
plot(srf(5).fish(11).x(tt{2}), srf(5).fish(11).y(tt{2}), '.', 'Color', [0.7,0.7,0.7]);

subplot(224); hold on;
plot(srf(5).fish(10).x(tt{1}), srf(5).fish(10).y(tt{1}), '.', 'Color', [0.7,0.7,0.7]);
plot(srf(5).fish(11).x(tt{2}), srf(5).fish(11).y(tt{2}), '.');

