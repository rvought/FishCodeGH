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
        plot(tmpsigA(tt(1:4:end))+0.5); 
        plot(tmpsigB(tt(1:4:end))-0.5);
    end
   
    drawnow;
    a = input('Best Frame? ');
    
    if isempty(a)
        ampdataOne = [ampdataOne, zeros(1, Fs*10)];
        ampdataTwo = [ampdataTwo, zeros(1, Fs*10)];
    end
    
    if ~isempty(a) % The fish is in the correct position in these frames
                           
          ampdataOne = [ampdataOne, tmpsigA(tim > samlen*(a-1) & tim <= samlen*a)];
          ampdataTwo = [ampdataTwo, data(tim > samlen*(a-1) & tim <= samlen*a)];
            
    end
    
% Mark the EOD data with a value to indicate light or dark    
%     ampdataOne(end+1) = mean(mean(vData(100:200, 500:600)));
%     ampdataTwo(end+1) = mean(mean(vData(100:200, 500:600)));
    ampdataOne(end+1) = mean(mean(vData(700:800, 500:600)));
    ampdataTwo(end+1) = mean(mean(vData(700:800, 500:600)));
    
    iidx = iidx+6;
    eidx = eidx+1;
    
end

