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

timtim = 1/Fs:1/Fs:length(ampdataOne)/Fs;
lightlevelidxs(1) = 1;
lightlevelidxs = [lightlevelidxs find(ampdataOne > 10)];

for j=2:length(lightlevelidxs)
    
    tt = lightlevelidxs(j-1)+1:lightlevelidxs(j)-1; % Range between pulses
    
    if sum(ampdataOne(tt)) > 0
        
        sondat = fftmachine(ampdataOne(tt), Fs);
        [famp, fidx] = max(sondat.fftdata);
        fftamp(:,j-1) = [famp, sondat.fftfreq(fidx)];
        rmsamp(j-1) = rms(ampdataOne(tt));
        lightlevel(j-1) = ampdataOne(lightlevelidxs(j));
        lightims(j-1) = timtim(lightlevelidxs(j));
    end    
    
    if sum(ampdataOne(tt)) == 0
        fftamp(:,j-1) = [];
        rmsamp(j-1) = [];
        lightlevel(j-1) = [];
        lightims(j-1) = timtim(lightlevelidxs(j));
    end
    
end

figure(2); clf; 
subplot(311); plot(lightims, fftamp);
subplot(312); plot(lightims, rmsamp)
subplot(312); plot(lightims, lightlevel)

