% load SurfaceDataRev2018a.mat
% load CaveDataRev2018a.mat
% Requires QuickPlotGrid, bFishAmpComparo, 


%% Example plots for fish movement versus frequency

    QuickPlotGrid(cave(12), [], [250 450], [240 360]); figure(2); axis([-75, 175, -200, 50]);
    pause(1);
    saveas(1,'caveFreqs.eps', 'epsc'); saveas(2,'cavePositions.eps', 'epsc');
    pause(1);
    QuickPlotGrid(srf(5), [], [250, 450], [120 240]); figure(2); axis([-145, 105, -125, 125]);
    pause(1);
    saveas(1,'surFreqs.eps', 'epsc'); saveas(2,'surfPositions.eps', 'epsc');
    pause(1);

    
%% EOD Amplitudes

stts = bFishAmpComparo(cave, srf);


%% dF versus distance


