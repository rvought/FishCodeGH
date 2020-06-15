function out = iu_sta(spikes, randspikes, sig, Fs, wid)

tim = 1/Fs:1/Fs:length(sig)/Fs; % Time stamps for the duration of the signal.

% For every spike (starting at the end) get the time "wid" before and after
% the time of the spike.
for idx = length(spikes):-1:1
    if spikes(idx) > wid && spikes(idx) < tim(end)-wid % Make sure that the window does not go before or after the signal.
        temp = interp1(tim, sig, spikes(idx)-wid:1/Fs:spikes(idx)+wid); % Copy the signal 
        sta(idx,:) = temp; % Put the signal into a temporary structure
    end
end

% For every randome spike (starting at the end) get the time "wid" before and after
% the time of the spike.
for idx = length(randspikes):-1:1
    if randspikes(idx) > wid && randspikes(idx) < tim(end)-wid % Make sure that the window does not go before or after the signal.
        temp = interp1(tim, sig, randspikes(idx)-wid:1/Fs:randspikes(idx)+wid);    
        sta_rand(idx,:) = temp;
    end
end

    out.MEAN  = nanmean(sta,1);
    out.STD  = nanstd(sta,0,1);

    out.randMEAN  = nanmean(sta_rand,1);
    out.randSTD  = nanstd(sta_rand,0,1);

    out.time = -wid:1/Fs:wid;

% figure,
% hold on,
% box on,
%     plot([0, 0], [min(out.MEAN), max(out.MEAN)], 'k-', 'LineWidth',1);
%     plot(out.time, out.MEAN, 'b-', 'LineWidth', 3);
%     plot(out.time, out.randMEAN,'r-','LineWidth',3);
%     xlabel('Time (s)')

 
