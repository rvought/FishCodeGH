% This script produces figure 5 for the manuscript -  interactions
% between dF and distance fish in cave and surface habitats. Relies on dFanalysis.m

% Calculate CAVE dFs TIME CONSUMING - ONLY DO THIS ONCE prior to running this script!
% load SurfaceDataRev2018a.mat; load CaveDataRev2018a.mat;
% [caveDF, CalldFs] = dFanalysis(cave);
% [srfDF, SalldFs] = dFanalysis(srf);

realCorrs = []; shuffCorrs = []; shiftCorrs = [];

Fs = 4.8828;

CorrWindow = 300; % Time in seconds for correlation analysis
StepSz = 150; % Time in seconds for the step

for kk = 3:4 % Placeholder for our data

for j = 1:length(data(kk).pair) % For each pair of fish
    
    tim = 1/Fs:1/Fs:length(data(kk).pair(j).descartes)/Fs;
    
    stepnum = floor( (tim(end)-StepSz) / StepSz );
    
    % Get the real Correlation coefficients and time-scrambled Correlation
    % coefficients.
    
    tmp = 1:stepnum;
    fakies = tmp(randperm(length(tmp))); 
    clear tmp;
    
    for z = 1:stepnum
        
        % Our time window
        tt = find(tim > StepSz * (z-1) & tim < (StepSz * (z-1)) + CorrWindow);

        % Our shuffled time window
        tf = find(tim > StepSz * (fakies(z)-1) & tim < (StepSz * (fakies(z)-1)) + CorrWindow);
        makethemthesamelength = min([length(tt), length(tf)]);
        tt = tt(1:makethemthesamelength);
        tf = tf(1:makethemthesamelength);
        
        % Our shifted time window
        if z ~= stepnum
        tt = find(tim > StepSz * z & tim < (StepSz * z) + CorrWindow);
        else
        tt = find(tim > 0 & tim < CorrWindow); % The first window
        end
        
        tmp = corrcoef(data(kk).pair(j).descartes(tt), data(kk).pair(j).dF(tt));
            realCorrs(end+1) = tmp(2);
        tmp = corrcoef(data(kk).pair(j).descartes(tt), data(kk).pair(j).dF(tf));
            shuffCorrs(end+1) = tmp(2);
        tmp = corrcoef(data(kk).pair(j).descartes(tt), data(kk).pair(j).dF(tf));
            shiftCorrs(end+1) = tmp(2);
            
        clear tmp;
    end
    
    
    
end

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





