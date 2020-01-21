% This script produces figure 5 for the manuscript -  interactions
% between dF and distance fish in cave and surface habitats. Relies on dFanalysis.m

% Calculate CAVE dFs TIME CONSUMING - ONLY DO THIS ONCE prior to running this script!
% load SurfaceDataRev2018a.mat; load CaveDataRev2018a.mat;
% [caveDF, CalldFs] = dFanalysis(cave);
% [srfDF, SalldFs] = dFanalysis(srf);

realCorrs = []; shuffCorrs = []; shiftCorrs = [];
idxPs = []; idxKs = [];
Fs = 4.8828;

CorrWindow = 200; % Time in seconds for correlation analysis
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
    end    
    
end 

end % For each pair

end % Has pair data

end % For each entry 

% Plot the data

figure(1); clf; 
cenbins = -1:0.2:1;
ax(1) = subplot(311); histogram(realCorrs, cenbins);
ax(2) = subplot(312); histogram(shuffCorrs, cenbins);
ax(3) = subplot(313); histogram(shiftCorrs, cenbins);
linkaxes(ax, 'xy'); xlim([-1,1]);


clear stepnum CorrWindow StepSz tt tf ts fakies z tim makethemthesamelength numbins j kk cenbins

% % % Fill in missing data.  This is dangerous - need reality check somewhere!!
% % 
% %         dist(dist == 0) = NaN;
% %         % dist = fillmissing(dist, 'pchip');
% %         dist = fillmissing(dist, 'linear','EndValues','nearest');
% %         
% %         dF(dF == 0) = NaN;
% %         %dF = fillmissing(dF, 'pchip');
% %         dF = fillmissing(dF, 'linear','EndValues','nearest');

% aaaa = fillmissing(eod1(timtim > curstarteod & timtim < curstarteod+wndo), 'linear','EndValues','nearest');

%         RReod = corrcoef(aaaa, bbbb);
%        curreodCorr(loopreod) = RReod(2);

[b,a] = butter(5, 0.1/(2*Fs), 'low');

% Plot extreme examples - this is correlation coefficient of 0.9048
kval = 6; pval = 1;
tim = 1/Fs:1/Fs:length(caveDF(kval).pair(pval).descartes)/Fs;
figure(2); clf; plot(tim, caveDF(kval).pair(pval).descartes, '.'); 
yyaxis right; plot(tim, caveDF(kval).pair(pval).dF, '.');
plot(tim, filtfilt(b,a, caveDF(kval).pair(pval).dF), '-k');

% Plot extreme examples - this is correlation coefficient of -0.9256
kval = 5; pval = 21;
tim = 1/Fs:1/Fs:length(caveDF(3).pair(pval).descartes)/Fs;
figure(3); clf; plot(tim, caveDF(kval).pair(pval).descartes, '.'); 
yyaxis right; plot(tim, caveDF(kval).pair(pval).dF, '.');
plot(tim, filtfilt(b,a, caveDF(kval).pair(pval).dF), '-k');



