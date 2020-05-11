%% Setup once per Matlab start

    addpath('~/SparkleShare/github.com/FishCodeGH/General')
    addpath('~/SparkleShare/github.com/FishCodeGH')
    % Relies on fftmachine.m

    % cd into the directory where the original files live
    
%% Setup once per directory

iFiles = dir('GallmanImage*.mat');
eFiles = dir('GallmanElectro*.mat');

Fs = 20000; % Sample rate for EOD in Hz
samlen = 10; % Duration of the sample (60 seconds total)

[b,a] = butter(3, 160 / (Fs/2), 'high'); % Highpass filter to remove 60Hz.
[d,c] = butter(3, 5000 / (Fs/2), 'low'); % Lowpass filter that should not be necessary.

if (length(iFiles) / 6) ~= length(eFiles)
    fprintf('Fuck off motherfucker\n');
end

%% Cycle through the files and build the structure "out"

winwidth = 2; % Duration of the window for analysis
stepsize = winwidth/2; % This is 50% overlap
iidx = 0;
eidx = 1;

while eidx <= length(eFiles)

    figure(1); clf; figure(2); clf;
    eval(['load ' eFiles(eidx).name]); % Load the EOD data
    tmpsigA = filtfilt(b,a,EODonly(:,1)); % Filter both channels
        tmpsigA = filtfilt(d,c,tmpsigA);
    tmpsigB = filtfilt(b,a,EODonly(:,2));
        tmpsigB = filtfilt(d,c,tmpsigB);
    
    fprintf('Entry %i. \n', eidx); % Tell the user where we are    
    
    for j=1:6 % For each 10 second epoch in the 60 second sample
        
        figure(1); 
        eval(['load ' iFiles(iidx+j).name]); % Load the image
        subplot(3,2,j); imshow(vData); % Plot the EODonly
        text(100,100, num2str(j), 'Color', 'g', 'FontSize', 24); % Add the label
        
        figure(2); 
        subplot(3,2,j); hold on;
        tt = find(tim > (j-1)*samlen &  tim <= j*samlen);
        plot(tim(tt(1:4:end)), 2*(tmpsigA(tt(1:4:end)))+0.5, 'b'); 
        plot(tim(tt(1:4:end)), 2*(tmpsigB(tt(1:4:end)))-0.5, 'm');
        text(0.5+((j-1)*10), 0, num2str(j), 'Color', 'k', 'FontSize', 24); % Add the label
        ylim([-1.5 1.5]);
        
        % Get the variance for both
            for zz = 1:(samlen/stepsize) - 1                 
                att = tt(tim(tt) > (tim(tt(1)) + stepsize*(zz-1)) & tim(tt) < (tim(tt(1)) + stepsize*(zz-1) + winwidth));
                
                tmpApeaks = findpeaks(tmpsigA(att) - mean(tmpsigA(att)));
                tmpAvar = var(tmpApeaks(tmpApeaks > 0));
                tmpBpeaks = findpeaks(tmpsigB(att) - mean(tmpsigB(att)));
                tmpBvar = var(tmpBpeaks(tmpBpeaks > 0));
                totalvar(zz) = tmpBvar + tmpAvar;
            end
                
                [minVar(j), minVaridx] = min(totalvar);
                minVartim(j) = tim(tt(1)) + winwidth*(minVaridx-1);
                
                plot([minVartim(j), minVartim(j)], [-1.5, 1.5], 'k-', 'LineWidth', 2);
                plot([minVartim(j)+winwidth, minVartim(j)+winwidth], [-1.5, 1.5], 'k-', 'LineWidth', 2);                
                
    end
    
    drawnow;
    
% Pick frame - OLD VERSION
% %     framNo = input('Best Frame? ');    
% %     if isempty(framNo) % User didn't click a frame    
% %         out(eidx).Ch1 = 0;
% %         out(eidx).Ch2 = 0;        
% %     end    
% %     if ~isempty(framNo) % The fish is in the correct position in these frames                          
% %         out(eidx).Ch1 = tmpsigA(tim > samlen*(framNo-1) & tim <= samlen*framNo);
% %         out(eidx).Ch2 = tmpsigB(tim > samlen*(framNo-1) & tim <= samlen*framNo);            
% %     end

% Click frame - NEW VERSION


[xPos, ~] = ginput(1);

    if ~isempty(xPos)
        out(eidx).Ch1 = tmpsigA(tim > xPos & tim <= xPos + winwidth);
        out(eidx).Ch2 = tmpsigB(tim > xPos & tim <= xPos + winwidth);        
    end
    if isempty(xPos) % User didn't click a frame    
        out(eidx).Ch1 = 0;
        out(eidx).Ch2 = 0;        
    end    


% Light, temp, and time data    
    out(eidx).light = mean(mean(vData(700:800, 500:600)));
    out(eidx).temp = temp;
    out(eidx).fileinfo = eFiles(eidx).name;
    
% Advance our counters
    iidx = iidx+6;
    eidx = eidx+1;
    
% Clear xPos
    clear xPos

end

out(1).Fs = Fs;

%% Analysis

Fs = out(1).Fs;

for j=length(out):-1:1
    
    if sum(out(j).Ch1) ~= 0
        
        tmp = fftmachine(out(j).Ch1, Fs);
        out(j).fftCh1data = tmp.fftdata; out(j).fftCh1freq = tmp.fftfreq;
            [famp, fidx] = max(smooth(out(j).fftCh1data, 10));
        out(j).fftCh1peakfreq = out(j).fftCh1freq(fidx);
        out(j).fftCh1peakamp = famp; 
        out(j).rmsCh1 = rms(out(j).Ch1);
        
        tmp = fftmachine(out(j).Ch2, Fs);
        out(j).fftCh2data = tmp.fftdata; out(j).fftCh2freq = tmp.fftfreq;
            [famp, fidx] = max(smooth(out(j).fftCh2data,10));
        out(j).fftCh2peakfreq = out(j).fftCh2freq(fidx);
        out(j).fftCh2peakamp = famp; 
        out(j).rmsCh2 = rms(out(j).Ch2);

        lightlevels(j) = out(j).light;
        temps(j) = out(j).temp;
        
        out(j).idx = j;
        
    end    
    
    if sum(out(j).Ch1) == 0
        
        out(j).Ch1 = []; % Empty out the 0'ed data 
        out(j).Ch2 = [];

        out(j).idx = [];
        out(j).rmsCh1 = [];
        out(j).rmsCh2 = [];
        out(j).fftCh1freq = []; out(j).fftCh1data = [];
        out(j).fftCh1peakfreq = []; out(j).fftCh1peakamp = [];
        out(j).fftCh2freq = []; out(j).fftCh2data = [];
        out(j).fftCh2peakfreq = []; out(j).fftCh2peakamp = [];
    end
        
end

f3 = figure(3); clf; 
    bb = [0 0 1]; mm = [1 0 1];
    set(f3, 'defaultAxesColorOrder', [bb; mm]);
    ax(1) = subplot(411); 
        yyaxis left; plot([out.idx], [out.fftCh1peakamp], 'b.-', 'MarkerSize', 8); title('FFT amplitude');
        hold on; 
        yyaxis right; plot([out.idx], [out.fftCh2peakamp], 'm.-', 'MarkerSize', 8); 
    ax(2) = subplot(412); 
        yyaxis right; plot([out.idx], [out.rmsCh1], 'b.-', 'MarkerSize', 8);
        hold on; yyaxis left; plot([out.idx], [out.rmsCh2], 'm.-', 'MarkerSize', 8); title('RMS amplitude');
    ax(3) = subplot(413); 
        plot([out.idx], [out.fftCh1peakfreq], '.-', 'MarkerSize', 8); ylim([200 700]); title('EOD Frequency');
        hold on; plot([out.idx], [out.fftCh2peakfreq], '.-', 'MarkerSize', 8); ylim([300 600]);
    ax(4) = subplot(414); 
        yyaxis left; plot(1:length(out), [out.light], '.-', 'MarkerSize', 8); ylim([180 260]);
        title('Light level & Temperature');
        yyaxis right; plot(1:length(out), [out.temp], '.-', 'MarkerSize', 8); ylim([0.7 0.9]);
        
    linkaxes(ax, 'x');
    
    
    
% figure(2); clf; 
%     ax(1) = subplot(411); plot(lightims, fftamp(:,1), '.-', 'MarkerSize', 8);
%     ax(2) = subplot(412); plot(lightims, rmsamp, '.-', 'MarkerSize', 8)
%     ax(3) = subplot(413); plot(lightims, fftamp(:,2), '.-', 'MarkerSize', 8);
%     ax(4) = subplot(414); plot(lightims, lightlevel, '.-', 'MarkerSize', 8)
%     linkaxes(ax, 'x');

%% Other stuff to do

% freeampEODonly = ampEODonly; freeampEODonly(freeampEODonly > 10) = 0;
% % % hdat = real(hilbert(freeampEODonly));
% hdat = envelope(freeampEODonly, 200, 'peak');
% figure(28); clf; plot(timtim(1:3:end), freeampEODonly(1:3:end), timtim(1:3:end), hdat(1:3:end)); ylim([-0.5 0.5]);
% % tt = find(timtim > 920 & timtim < 930);
% % asdf = fftmachine(hdat(tt(1:10:end)), Fs/10);
% % figure(29); hold on; semilogy(asdf.fftfreq, asdf.fftdata); xlim([0 100]);


%% When you are done...

%
% save filname.mat out