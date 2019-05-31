function out = getTE(dFdist, windowlength, stepsize)

delays = [50, 60, 70, 80, 90];
    
for j = length(dFdist):-1:1
    
    for z = 1:length(delays)
    [out(j).TE{z}, out(j).tt{z}] = calcTE(dFdist(j), windowlength, stepsize, delays(z));
    end
    
end


%% Plot
    for j=1:length(out)
        figure(j+10); clf; hold on;
        for k = 1:length(out(j).TE)
            plot(out(j).tt{k}, out(j).TE{k});
        end
    end

%% Embedded calcTE function
function [currTE, currTT] = calcTE(data, windo, stp, kk)

ll = 1; 

currTE = []; currTT = [];

strts = 0:stp:data.tim(end)-windo;


parfor loopr = 1:length(strts)
    
      curstart = strts(loopr);
      aaa= data.dF(data.tim > curstart & data.tim < curstart+windo);
      bbb = data.distance(data.tim > curstart & data.tim < curstart+windo);
      [currTE(loopr),~ ,~] = transferEntropyPartition(aaa(1:2:end), bbb(1:2:end), ll, kk);
      
      [currTE(loopr),~ ,~] = transferEntropyPartition(data.dF(data.tim > curstart & data.tim < curstart+windo), data.distance(data.tim > curstart & data.tim < curstart+windo), ll, kk);      
      %currTE(loopr) = transferEntropyKDE(data.dF(data.tim > curstart & data.tim < curstart+windo), data.distance(data.tim > curstart & data.tim < curstart+windo), ll, kk, 2, 2); 
      %currTE(loopr) = transferEntropyRank(data.dF(data.tim > curstart & data.tim < curstart+windo), data.distance(data.tim > curstart & data.tim < curstart+windo), ll, kk, 2, 2, 10);

      currTT(loopr) = curstart + (windo/2);
                 
end


end

end
