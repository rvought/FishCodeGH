function preTExplorer(dFdist, orig)
% ESF
windw = 300;
advanceby = 10;

trialstouse=[]; emptytrials= [];
for j=1:length(dFdist)
    if ~isempty(dFdist(j).dF); trialstouse(end+1) = j; end
    if isempty(dFdist(j).dF); emptytrials(end+1) = j; end
end

fprintf('These are the bad trials: %2d \n', emptytrials);

for j = trialstouse
        
    figure(1); clf;
    
    subplot(411); % plot original frequency data
    hold on;
        plot(orig(dFdist(j).fishnums(1)).tim, orig(dFdist(j).fishnums(1)).EOD, '.b', 'MarkerSize', 2);
        plot(orig(dFdist(j).fishnums(1)).tim, orig(dFdist(j).fishnums(2)).EOD, '.m', 'MarkerSize', 2);
    xlim([0 orig(dFdist(j).fishnums(1)).tim(end)]);
    ylim([200 500]);
    text(100, 450, num2str(dFdist(j).fishnums));
    text(100, 400, num2str(j));

    subplot(412); % plot Distance and dF data
    xlim([0 orig(dFdist(j).fishnums(1)).tim(end)]);
    hold on;
    yyaxis right;
        plot(dFdist(j).tim, dFdist(j).distance, '.');
    ylabel('distance')
    yyaxis left;
        plot(dFdist(j).tim, dFdist(j).dF, '.');
    ylabel('dF')
    
    subplot(425); hold on;
    plot(orig(dFdist(j).fishnums(2)).xy(:,1), orig(dFdist(j).fishnums(2)).xy(:,2), '.', 'MarkerSize', 4, 'Color', [0.7 0.7 0.7]);
    plot(orig(dFdist(j).fishnums(1)).xy(:,1), orig(dFdist(j).fishnums(1)).xy(:,2), '.b', 'MarkerSize', 4);
    xlim([-250, 250]);
    ylim([-250, 250]);
    
    subplot(426); hold on;
    plot(orig(dFdist(j).fishnums(1)).xy(:,1), orig(dFdist(j).fishnums(1)).xy(:,2), '.', 'MarkerSize', 4, 'Color', [0.7 0.7 0.7]);
    plot(orig(dFdist(j).fishnums(2)).xy(:,1), orig(dFdist(j).fishnums(2)).xy(:,2), '.m', 'MarkerSize', 4);
    xlim([-250, 250]);
    ylim([-250, 250]);
    
    subplot(4,1,4); hold on;
        tmpdF = dFdist(j).dF - nanmean(dFdist(j).dF);
        tmpdF = tmpdF/(max(abs(tmpdF)));
        tmpDistance = dFdist(j).distance - nanmean(dFdist(j).distance);
        tmpDistance = tmpDistance / (max(abs(tmpDistance)));
        tims = dFdist(j).tim;
        
        [CC, TT] = slideCorr(tmpdF, tmpDistance, tims, windw, advanceby);
        
        plot([0, max(tims)], [0.7, 0.7], 'r');
        plot([0, max(tims)], [-0.7, -0.7], 'r');
        plot([0, max(tims)], [0.5, 0.5], 'm');
        plot([0, max(tims)], [-0.5, -0.5], 'm');
        plot([0, max(tims)], [0, 0], 'g', 'LineWidth', 6);
        
        plot(TT, CC, '.k-', 'MarkerSize', 4);
        
        xlim([0, max(tims)]); ylim([-1, 1]);
        text(100, 0, num2str(mean(CC)));

%     xa = xcorr(tmpdF, tmpDistance);
%     plot(xa); xlim([0 length(xa)]);
%     [mm, midx] = max(xa);
%     plot([length(xa)/2, length(xa)/2], [0 mm], 'r');
%     plot(midx, mm, 'go');
%     [h,~] = corrcoef(dFdist(j).dF-mean(dFdist(j).dF), dFdist(j).distance-mean(dFdist(j).distance));
%     text(length(xa)/2, 0, num2str(h(2)));
%     ylim([-1000, 1000]);

    endogo = input('foobar 9 to end: ');  
    if endogo == 9; break; end
    
end

%% Embedded function slideCorr
function [currCorr, currTT] = slideCorr(dF, dist, tim, windo, stp)

    length(dF)
    
strts = 0:stp:tim(end)-windo;

 for loopr = 1:length(strts)
    
      curstart = strts(loopr);
      aaa = dF(tim > curstart & tim < curstart+windo);
      bbb = dist(tim > curstart & tim < curstart+windo);
      
        RR = corrcoef(aaa, bbb);
        currCorr(loopr) = RR(2);

      currTT(loopr) = curstart + (windo/2);
                 
 end


end % End of embedded function


end % End of function

    