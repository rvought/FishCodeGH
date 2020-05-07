% This script produces figure 5 for the manuscript -  interactions
% between dF and distance fish in cave and surface habitats. Relies on dFanalysis.m

% Calculate CAVE dFs TIME CONSUMING - ONLY DO THIS ONCE prior to running this script!
% load SurfaceDataRev2018a.mat; load CaveDataRev2018a.mat;
% [caveDF, CalldFs] = dFanalysis(cave);
% [srfDF, SalldFs] = dFanalysis(srf);

% data = caveDF;
% data = srfDF;


realCorrs = []; shuffCorrs = []; shiftCorrs = [];
idxPs = []; idxKs = []; timStarts = []; meandF = [];
Fs = 4.8828; % This is the sample rate that emerged from the grid analysis.

CorrWindow = 300; % Time in seconds for correlation analysis
StepSz = 100; % Time in seconds for the step

%for kk = 3 % Placeholder for our data
for kk = 1:length(data)
    
if ~isempty(data(kk).pair)
    
for j = 1:length(data(kk).pair) % For each pair of fish
    
    tim = 1/Fs:1/Fs:length(data(kk).pair(j).descartes)/Fs;
    
    if ~isempty(tim)
    
    stepnum = floor( (tim(end)-StepSz) / StepSz );
    
    % Get the real Correlation coefficients and time-scrambled Correlation`
    % coefficients.
    
    tmp = 1:stepnum;
    fakies = tmp(randperm(length(tmp))); 
    clear tmp;
    
    for z = 1:stepnum
        
        % Our time window
        tt = find(tim > StepSz * (z-1) & tim < (StepSz * (z-1)) + CorrWindow);

        % Our shuffled time window
        tf = find(tim > StepSz * (fakies(z)-1) & tim < (StepSz * (fakies(z)-1)) + CorrWindow);
        
        % Our shifted time window
        if z ~= stepnum
        ts = find(tim > StepSz * z & tim < (StepSz * z) + CorrWindow);
        else
        ts = find(tim > 0 & tim < CorrWindow); % The first window
        end
        
        makethemthesamelength = min([length(tt), length(tf), length(ts)]);
        tt = tt(1:makethemthesamelength);
        tf = tf(1:makethemthesamelength);
        ts = ts(1:makethemthesamelength);

        
        tmp = corrcoef(data(kk).pair(j).descartes(tt), data(kk).pair(j).dF(tt));
            realCorrs(end+1) = tmp(2);
        tmp = corrcoef(data(kk).pair(j).descartes(tt), data(kk).pair(j).dF(tf));
            shuffCorrs(end+1) = tmp(2);
        tmp = corrcoef(data(kk).pair(j).descartes(tt), data(kk).pair(j).dF(ts));
            shiftCorrs(end+1) = tmp(2);
        clear tmp;
        
        idxPs(end+1) = j;
        idxKs(end+1) = kk;
        timStarts(end+1) = StepSz * (z-1);
        
        meandF(end+1) = mean(data(kk).pair(j).dF(tt));
        
    end    
    
    end % to make sure...

end % For each pair

end % Has pair data

end % For each entry 

%% Plot the histgram data

%figure(1); clf; 
figure(3); clf;
%cenbins = -1:0.33333:1;
%pltbins = -0.83335:0.3333:0.83335;
cenbins = -1:0.2:1;
pltbins = -0.9:0.2:0.9;
% ax(1) = subplot(311); histogram(realCorrs, cenbins);
% ax(2) = subplot(312); histogram(shuffCorrs, cenbins);
% ax(3) = subplot(313); histogram(shiftCorrs, cenbins);
% linkaxes(ax, 'xy'); xlim([-1,1]);

    histogram(shuffCorrs, cenbins, 'FaceColor', 'red'); 
%    histogram(shiftCorrs, cenbins, 'FaceColor', 'red'); 
        hold on; 
    histogram(realCorrs, cenbins, 'FaceColor', 'blue');
    
% TESTING significance.  We have 3908 samples for CAVE.

% MORE STRONG POSITIVE CORRELATIONS than expected?
posthresh = 0.6;
gtgt(1,:) = [length(find(shuffCorrs > posthresh)), length(find(shuffCorrs < posthresh))];
%gtgt(1,:) = [length(find(shiftCorrs > posthresh)), length(find(shiftCorrs < posthresh))];
gtgt(2,:) = [length(find(realCorrs > posthresh)), length(find(realCorrs < posthresh))];
[~,p,stats] = fishertest(gtgt);
fprintf('Fisher Test more POSITIVE correlations than expected? p=%2.8f, oddsratio=%2.4f \n', p, stats.OddsRatio); 

% MORE STRONG NEGATIVE CORRELATIONS than expected?
negthresh = -posthresh;
ltlt(1,:) = [length(find(shuffCorrs < negthresh)), length(find(shuffCorrs > negthresh))];
%ltlt(1,:) = [length(find(shiftCorrs < negthresh)), length(find(shiftCorrs > negthresh))];
ltlt(2,:) = [length(find(realCorrs < negthresh)), length(find(realCorrs > negthresh))];
[~,p,stats] = fishertest(ltlt);
fprintf('Fisher Test more NEGATIVE correlations than expected? p=%2.8f, oddsratio=%2.4f \n', p, stats.OddsRatio); 


%figure(2); clf;
figure(4); clf;

totalnum = length(realCorrs);

    shuffHist = histcounts(shuffCorrs, cenbins); 
    shiftHist = histcounts(shiftCorrs, cenbins); 
        hold on; 
    realHist = histcounts(realCorrs, cenbins);
    plot(pltbins, shuffHist/totalnum, '.r-', 'MarkerSize', 16); 
    plot(pltbins, shiftHist/totalnum, '.m-', 'MarkerSize', 16); 
        hold on;
    plot(pltbins, realHist/totalnum, '.b-', 'MarkerSize', 16);
    ylim([0 0.31]);
    


[h,p,ci,stats] = vartest2(realCorrs, shuffCorrs);
fprintf('vartest2 realCorrs versus shuffCorss p=%1.6f\n', p);

fprintf('Mean dF=%3.4f and std=%3.4f and N=%i\n', mean(meandF), std(meandF), length(meandF));

clear stepnum tt tf ts fakies z tim makethemthesamelength numbins j kk cenbins


%% What is the phase lag between movement and dF for highly correlated epochs
% posShift = []; negShift = [];
% posCorrThresh = 0.85;
% negCorrThresh = -0.85;
% posIDX = find(realCorrs > posCorrThresh);
% negIDX = find(realCorrs < negCorrThresh);
% 
% figure(4); clf; hold on;
% 
% for z=1:length(posIDX)    
%     tim = 1/Fs:1/Fs:length(data(idxKs(z)).pair(idxPs(z)).descartes);
%     tt = find(tim > timStarts(z) & tim < timStarts(z) + CorrWindow);
%     
%     descartesdata = data(idxKs(z)).pair(idxPs(z)).descartes(tt) - mean(data(idxKs(z)).pair(idxPs(z)).descartes(tt));
%     dFdata = data(idxKs(z)).pair(idxPs(z)).dF(tt) - mean(data(idxKs(z)).pair(idxPs(z)).dF(tt));
%     xc = xcorr(descartesdata, dFdata);
%     % plot(abs(xc)/max(abs(xc)));
%     wid = length(xc);
%     [~, maxidx] = max(abs(xc)); 
%     if maxidx > 0.4*length(xc) && maxidx < 0.6*length(xc)
%         posShift(end+1) = maxidx - length(data(idxKs(z)).pair(idxPs(z)).descartes(tt));
%     end
% end
% for z=1:length(negIDX)    
%     tim = 1/Fs:1/Fs:length(data(idxKs(z)).pair(idxPs(z)).descartes);
%     tt = find(tim > timStarts(z) & tim < timStarts(z) + CorrWindow);
%     
%     descartesdata = data(idxKs(z)).pair(idxPs(z)).descartes(tt) - mean(data(idxKs(z)).pair(idxPs(z)).descartes(tt));
%     dFdata = data(idxKs(z)).pair(idxPs(z)).dF(tt) - mean(data(idxKs(z)).pair(idxPs(z)).dF(tt));
%     xc = xcorr(descartesdata, dFdata);
%     plot(abs(xc)/max(abs(xc)));
%     wid = length(xc);
%     [~, maxidx] = max(abs(xc)); 
%     if maxidx > 0.4*length(xc) && maxidx < 0.6*length(xc)
%         negShift(end+1) = maxidx - length(data(idxKs(z)).pair(idxPs(z)).descartes(tt));
%     end
% end
% 
% 
% % % % Fill in missing data.  This is dangerous - need reality check somewhere!!
% % % 
% % %         dist(dist == 0) = NaN;
% % %         % dist = fillmissing(dist, 'pchip');
% % %         dist = fillmissing(dist, 'linear','EndValues','nearest');
% % %         
% % %         dF(dF == 0) = NaN;
% % %         %dF = fillmissing(dF, 'pchip');
% % %         dF = fillmissing(dF, 'linear','EndValues','nearest');
% 
% % aaaa = fillmissing(eod1(timtim > curstarteod & timtim < curstarteod+wndo), 'linear','EndValues','nearest');
% 
% %         RReod = corrcoef(aaaa, bbbb);
% %        curreodCorr(loopreod) = RReod(2);
% 
% [b,a] = butter(5, 0.1/(2*Fs), 'low');
% 
% % Plot excellent examples - this is correlation coefficient of 0.XXXX
% kval = 12; pval = 12;
% tim = 1/Fs:1/Fs:length(caveDF(kval).pair(pval).descartes)/Fs;
% figure(3); clf; 
% subplot(211); hold on;
% yyaxis left; plot(tim, caveDF(kval).pair(pval).descartes, '.b'); 
% plot(tim, filtfilt(b,a, caveDF(kval).pair(pval).descartes), '-b');
% yyaxis right; plot(tim, caveDF(kval).pair(pval).dF, '.m');
% plot(tim, filtfilt(b,a, caveDF(kval).pair(pval).dF), '-m');
% xlim([600, 1100]);
% 
% % Plot extreme examples - this is correlation coefficient of -XXXX
% kval = 14; pval = 18;
% tim = 1/Fs:1/Fs:length(caveDF(kval).pair(pval).descartes)/Fs;
% % figure(3); clf; hold on;
% subplot(212); hold on;
% yyaxis left; plot(tim, caveDF(kval).pair(pval).descartes, '.b'); 
% plot(tim, filtfilt(b,a, caveDF(kval).pair(pval).descartes), '-b');
% yyaxis right; plot(tim, caveDF(kval).pair(pval).dF, '.m');
% plot(tim, filtfilt(b,a, caveDF(kval).pair(pval).dF), '-m');
% xlim([150, 650]);
% 
% 

