function out = painkiller(data)
for kk = 1:length(data)
    %% calculating the ISI
    ISI = ISIfinder(data(kk).s);
    
    %% creating histograms and heatmaps
    hotstuff = sunburn(data(kk).s);
    
    % in order to call data forms that are made in the creation of the
    % heatmaps and histograms, call out.heathist and you have the options
    % posvel, tmpPV, accvel, tmpAV, pos, vel, and acc
    
     %% calculating DSI and nDSI (must come after the heat and histo because need the output of the heatmaps)
     DSI = killerDSI(hotstuff);
     nDSI = objDSI(hotstuff);  %%this will be done for the size that is selected for the heatmaps and histograms
     
     %% completing the STA
     [spiketa, dat] = sta(hotstuff.spikes, hotstuff.objpos, hotstuff.Fs); 
     
     
end
out.ISI = ISI;
out.hotstuff = hotstuff;
out.DSI = DSI;
out.nDSI = nDSI;
out.sta = [spiketa, dat];
end 