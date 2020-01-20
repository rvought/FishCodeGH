% This script produces figure 5 for the manuscript -  interactions
% between dF and distance fish in cave and surface habitats. Relies on dFanalysis.m

% Calculate CAVE dFs TIME CONSUMING - ONLY DO THIS ONCE prior to running this script!
% load SurfaceDataRev2018a.mat; load CaveDataRev2018a.mat;
% [caveDF, CalldFs] = dFanalysis(cave);
% [srfDF, SalldFs] = dFanalysis(srf);

kk = 1;

for j = 1:length(data(kk).pair) % For each pair of fish
    
    
    
    
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





