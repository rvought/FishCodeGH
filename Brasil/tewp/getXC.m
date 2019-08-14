function out = getXC(dFdist, windowlength, stepsize)

% delays = [50, 60, 70, 80, 90];
    
for j = length(dFdist):-1:1
    
    [out(j).CC, out(j).tt] = calcCORR(dFdist(j), windowlength, stepsize);

%     for z = 1:length(delays)
%     [out(j).TE{z}, out(j).tt{z}] = calcCORR(dFdist(j), windowlength, stepsize, delays(z));
%     end
    
end


%% Plot
    for j=1:length(out)
        figure(j+10); clf; hold on;
        for k = 1:length(out(j).TE)
            plot(out(j).tt{k}, out(j).TE{k});
        end
    end

%% Embedded calcCORR function
function [currRR, currTT] = calcCORR(data, windo, stp)

strts = 0:stp:data.tim(end)-windo;

for loopr = 1:length(strts)
    
      curstart = strts(loopr);
      aaa= data.dF(data.tim > curstart & data.tim < curstart+windo);
      bbb = data.distance(data.tim > curstart & data.tim < curstart+windo);
      
      RR = corrcoef(aaa, bbb);
      
      currRR = RR(2);
      currTT(loopr) = curstart + (windo/2);
                 
end


end

end
