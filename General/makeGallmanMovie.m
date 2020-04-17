%% Picking

iFiles = dir('GallmanImage*.mat');
eFiles = dir('GallmanElectro*.mat');

ampdataOne = [];
ampdataTwo = [];
Fs = 20000; % Sample rate for EOD in Hz
samlen = 10; % Duration of the sample (60 seconds total)

[b,a] = butter(3, 160 / (Fs/2), 'high'); % Highpass filter to remove 60Hz.

% out = VideoWriter('video.avi');
% open(out);

% for i=1:length(iFiles)
%     eval(['load ' files(i).name]);
%     imshow(vData); hold on; text(100,100, num2str(i), 'Color', 'r');
%     tmp = getframe(gcf);
%     writeVideo(out, tmp);
%     pause(0.1); clf;
%     vid(i).vData = vData;
% end
% 
% close(out);

if (length(iFiles) / 6) ~= length(eFiles)
    fprintf('Fuck off motherfucker\n');
end

% for i=1:length(files)

iidx = 0;
eidx = 1;

while eidx <= length(eFiles)

    figure(1); clf; figure(2); clf;
    eval(['load ' eFiles(eidx).name]); % Load the EOD data
    tmpsigA = filtfilt(b,a,data(:,1)); % Filter both channles
    tmpsigB = filtfilt(b,a,data(:,2));
    
    fprintf('Entry %i. \n', eidx); % Tell the user where we are    
    
    for j=1:6 % For each 10 second epoch in the 60 second sample
        figure(1); 
        eval(['load ' iFiles(iidx+j).name]); % Load the image
        subplot(3,2,j); imshow(vData); % Plot the data
        text(100,100, num2str(j), 'Color', 'g', 'FontSize', 24); % Add the label
        figure(2); 
        subplot(3,2,j); hold on;
        tt = find(tim > (j-1)*samlen &  tim <= j*samlen);
        plot(tim(tt(1:4:end)), tmpsigA(tt(1:4:end))+0.5); 
        plot(tim(tt(1:4:end)), tmpsigB(tt(1:4:end))-0.5);
        text(0.5+((j-1)*10), 0, num2str(j), 'Color', 'm', 'FontSize', 24); % Add the label
        ylim([-1.5 1.5]);
    end
   
    drawnow;
    framNo = input('Best Frame? ');
    
    if isempty(framNo)
        ampdataOne = [ampdataOne, zeros(1, Fs*10)];
        ampdataTwo = [ampdataTwo, zeros(1, Fs*10)];
    end
    
    if ~isempty(framNo) % The fish is in the correct position in these frames
                           
          ampdataOne = [ampdataOne, tmpsigA(tim > samlen*(framNo-1) & tim <= samlen*framNo)'];
          ampdataTwo = [ampdataTwo, data(tim > samlen*(framNo-1) & tim <= samlen*framNo)'];
            
    end
    
% Mark the EOD data with a value to indicate light or dark    
%     ampdataOne(end+1) = mean(mean(vData(100:200, 500:600)));
%     ampdataTwo(end+1) = mean(mean(vData(100:200, 500:600)));
    ampdataOne(end+1) = mean(mean(vData(700:800, 500:600)));
    ampdataTwo(end+1) = mean(mean(vData(700:800, 500:600)));
    
    iidx = iidx+6;
    eidx = eidx+1;
    
end

%% Analysis

ampdata = ampdataOne - mean(ampdataOne);

timtim = 1/Fs:1/Fs:length(ampdata)/Fs;
lightlevelidxs(1) = 1;
lightlevelidxs = [lightlevelidxs(1) find(ampdata > 10)];

for j=2:length(lightlevelidxs)
    
    tt = lightlevelidxs(j-1)+1:lightlevelidxs(j)-1; % Range between pulses
    lightims(j-1) = timtim(lightlevelidxs(j));

    if sum(ampdata(tt)) ~= 0
        
        sondat = fftmachine(ampdata(tt), Fs, 50);
            [famp, fidx] = max(sondat.fftdata);
        fftamp(j-1,:) = [famp, sondat.fftfreq(fidx)];
        rmsamp(j-1) = rms(ampdata(tt));
        lightlevel(j-1) = ampdata(lightlevelidxs(j));
        
    end    
        
end

gIDX = find(rmsamp ~= 0);

figure(2); clf; 
    ax(1) = subplot(411); plot(lightims(gIDX), fftamp(gIDX,1), '.-', 'MarkerSize', 8);
    ax(2) = subplot(412); plot(lightims(gIDX), rmsamp(gIDX), '.-', 'MarkerSize', 8)
    ax(3) = subplot(413); plot(lightims(gIDX), fftamp(gIDX,2), '.-', 'MarkerSize', 8);
    ax(4) = subplot(414); plot(lightims(gIDX), lightlevel(gIDX), '.-', 'MarkerSize', 8)
    linkaxes(ax, 'x');
% figure(2); clf; 
%     ax(1) = subplot(411); plot(lightims, fftamp(:,1), '-.', 'MarkerSize', 8);
%     ax(2) = subplot(412); plot(lightims, rmsamp, '.', 'MarkerSize', 8)
%     ax(3) = subplot(413); plot(lightims, fftamp(:,2), '.', 'MarkerSize', 8);
%     ax(4) = subplot(414); plot(lightims, lightlevel, '.', 'MarkerSize', 8)
%     linkaxes(ax, 'x');

freeampdata = ampdata; freeampdata(freeampdata > 10) = 0;
% hdat = real(hilbert(freeampdata));
hdat = envelope(freeampdata, 200, 'peak');
figure(28); clf; plot(timtim, freeampdata, timtim, hdat);
tt = find(timtim > 920 & timtim < 930);
asdf = fftmachine(hdat(tt(1:10:end)), Fs/10);
figure(29); hold on; semilogy(asdf.fftfreq, asdf.fftdata); xlim([0 100]);



