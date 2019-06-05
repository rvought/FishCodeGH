function out = getCorr(dFdist, orig, windowlength, stepsize)
% ESF

% delays = [50, 60, 70, 80, 90];
% delays = 50;
% delays = [10 40 60 80 120];
%    out.Corr =[]; out.dF = []; out.distance = []; out.tt = [];
for j = length(dFdist):-1:1
    
%     for z = 1:length(delays)
%         [out(j).TE{z}, out(j).tt{z}] = calcTE(dFdist(j), windowlength, stepsize, delays(z));
%     end

[aa,bb,cc,dd] = slideCorr(dFdist(j).dF, dFdist(j).distance, dFdist(j).tim, windowlength, stepsize);
    
    out(j).Corr = aa;
    out(j).meandF = bb;
    out(j).meanDist = cc;
    out(j).tt = dd;
    
    fa = dFdist(j).fishnums(1);
    fb = dFdist(j).fishnums(2);
    
    [aaeod,bbeod] = eodCorr(orig(fa).EOD, orig(fb).EOD, dFdist(j).tim, windowlength, stepsize);

    out(j).eodCorr = aaeod;
    out(j).eodtt = bbeod;
    
end

%% Plot
%     if length(delays) > 1
%     for j=1:length(out)
%         figure(j+10); clf; hold on;
%         for k = 1:length(out(j).TE)
%             plot(out(j).tt{k}, out(j).TE{k});
%         end
%     end
%     end

figure(5); clf; 
        subplot(211); hold on;
            for j=1:length(out)
                plot(out(j).tt, out(j).Corr);
            end
            
            alldists = []; alldFs = []; maxcorrs = []; avgcorrs = [];
            for j=1:length(out)
                alldists = [alldists, out(j).meanDist];
                alldFs = [alldFs, out(j).meandF];
                [~, idx] =  max(abs(out(j).Corr));
                maxcorrs = [maxcorrs, out(j).Corr(idx)];
                avgcorrs = [avgcorrs, median(out(j).Corr)];
            end
            
        subplot(223); hold on;
            plot(maxcorrs, alldFs, 'ro');
            plot(avgcorrs, alldFs, 'co');
        subplot(224); hold on;
            plot(maxcorrs, alldists, 'ro');
            plot(avgcorrs, alldists, 'co');
        
            figure(6); clf; hold on;
            plot3(alldFs, alldists, maxcorrs, 'ro');
            plot3(alldFs, alldists, avgcorrs, 'co');

figure(6); clf; 
        subplot(211); hold on;
            for j=1:length(out)
                plot(out(j).eodtt, out(j).eodCorr);
            end
            
    
    
%% Embedded function slideCorr
function [currCorr, meandF, meanDist, currTT] = slideCorr(dF, dist, tim, windo, stp)

strts = 0:stp:tim(end)-windo;

for loopr = 1:length(strts)
    
      curstart = strts(loopr);
      aaa = dF(tim > curstart & tim < curstart+windo);
      bbb = dist(tim > curstart & tim < curstart+windo);
      
        RR = corrcoef(aaa, bbb);
        currCorr(loopr) = RR(2);
        meandF = mean(aaa);
        meanDist = mean(bbb);

      currTT(loopr) = curstart + (windo/2);
                 
 end


end % End of embedded function

%% Embedded function eodCorr
function [curreodCorr, curreodTT] = eodCorr(eod1, eod2, timtim, wndo, stpstp)

strtseod = 0:stpstp:timtim(end)-wndo;


for loopreod = 1:length(strtseod)
    
      curstarteod = strtseod(loopreod);
      aaaa = fillmissing(eod1(timtim > curstarteod & timtim < curstarteod+wndo), 'linear','EndValues','nearest');
      bbbb = fillmissing(eod2(timtim > curstarteod & timtim < curstarteod+wndo), 'linear','EndValues','nearest');
      
      aaaa = aaaa - min(aaaa); bbbb = bbbb - min(bbbb);
      
        RReod = corrcoef(aaaa, bbbb);
        
        curreodCorr(loopreod) = RReod(2);
        curreodTT(loopreod) = curstarteod + (wndo/2);
                 
 end

curreodCorr = fillmissing(curreodCorr, 'linear','EndValues','nearest');

end % End of eodCorr


end