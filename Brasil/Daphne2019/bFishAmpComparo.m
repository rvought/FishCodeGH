function stts = bFishAmpComparo(cfish, sfish)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

OutlierLevel = 0.005;

CaveAmps = []; SurfaceAmps = []; 

for j=1:length(sfish) 
    for k=1:length(sfish(j).fish) 
        SurfaceAmps(end+1) = sfish(j).fish(k).dipStrength; 
    end
end

for j=1:length(cfish)
    for k=1:length(cfish(j).fish)
        CaveAmps(end+1) = cfish(j).fish(k).dipStrength; 
    end
end

stts.meanSurfaceAmp = mean(SurfaceAmps(SurfaceAmps < OutlierLevel));
stts.meanCaveAmp = mean(CaveAmps(CaveAmps < OutlierLevel));
stts.stdSurfaceAmp = std(SurfaceAmps(SurfaceAmps < OutlierLevel));
stts.stdCaveAmp = std(CaveAmps(CaveAmps < OutlierLevel));

[stts.H,stts.P,stts.CI,stts.STATS] = ttest2(CaveAmps(CaveAmps < OutlierLevel), SurfaceAmps(SurfaceAmps <OutlierLevel));

    fprintf('Amplitudes different between cave and surface pVal = %1.4f \n', stts.P);
    fprintf('Cave mean & std %1.4f %1.4f \n', stts.meanCaveAmp*100, stts.stdCaveAmp*100);
    fprintf('Surface mean & std %1.4f %1.4f \n', stts.meanSurfaceAmp*100, stts.stdSurfaceAmp*100);


figure(1); clf;
    ctrs = 0:0.0001:0.004;
    ax(1) = subplot(211); histogram(CaveAmps(CaveAmps < OutlierLevel),ctrs); 
        text(0.002,20, 'CaveFish EOD dipole Strengths');
    ax(2) = subplot(212); histogram(SurfaceAmps(SurfaceAmps < OutlierLevel), ctrs);
        text(0.002,20, 'Surface Fish EOD dipole Strengths');
    linkaxes(ax, 'xy');
    xlim([0 0.004]);

end

