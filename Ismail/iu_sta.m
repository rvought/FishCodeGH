function out = iu_sta(spikes, randspikes, sig, Fs, wid)
% Function out = iu_sta(spikes, randspikes, sig, Fs, wid)
% spikes are the spike times
% randspikes are shuffled spike times
% sig is the signal (e.g. error_vel) of interest. Behavior...
% Fs is the sample rate (usually 25 for these data, fs = 25
% wid is the width of the spike triggered average in seconds (1 or 2 seconds is good)

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

    out.MEAN  = nanmean(sta,1); % Calculate the mean (which is the STA)
    out.STD  = nanstd(sta,0,1); % Get the standard deviation for each point.

    out.randMEAN  = nanmean(sta_rand,1);
    out.randSTD  = nanstd(sta_rand,0,1);

    out.time = -wid:1/Fs:wid; % Give the user a time base for plotting.

% figure,
% hold on,
% box on,
%     plot([0, 0], [min(out.MEAN), max(out.MEAN)], 'k-', 'LineWidth',1);
%     plot(out.time, out.MEAN, 'b-', 'LineWidth', 3);
%     plot(out.time, out.randMEAN,'r-','LineWidth',3);
%     xlabel('Time (s)')

 
