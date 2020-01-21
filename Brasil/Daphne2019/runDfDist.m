% This script produces figure 5 for the manuscript -  interactions
% between dF and distance fish in cave and surface habitats. Relies on dFanalysis.m

% Calculate CAVE dFs TIME CONSUMING - ONLY DO THIS ONCE prior to running this script!
% load SurfaceDataRev2018a.mat; load CaveDataRev2018a.mat;
% [caveDF, CalldFs] = dFanalysis(cave);
% [srfDF, SalldFs] = dFanalysis(srf);
realCorrs = []; fakeCorrs = [];
Fs = 4.8828;
CorrWindow = 300; % Time in seconds for correlation analysis
StepSz = 150; % Time in seconds for the step

kk = 3; % Placeholder for our data

for j = 1:length(data(kk).pair) % For each pair of fish
    
    tim = 1/Fs:1/Fs:length(data(kk).pair(j).descartes)/Fs;
    
    stepnum = floor( (tim(end)-StepSz) / StepSz );
    
    % Get the real Correlation coefficients
    for z = 1:stepnum
        tt = find(tim > StepSz * (z-1) & tim < (StepSz * (z-1)) + CorrWindow);
        realCorrs(end+1) = corrcoef(data(kk).pair(j).descartes(tt), data(kk).pair(j).dF(tt));
    end
    
    % Get time scrambled (fake) Correlations (bootstrap)
    
    
end



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





