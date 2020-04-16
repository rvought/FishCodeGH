iFiles = dir('GallmanImage*.mat');
eFiles = dir('GallmanElectro*.mat');

ampdataOne = [];
ampdataTwo = [];
Fs = 20000;

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

while eidx < length(eFiles)

    fprintf('Entry %i. \n', eidx);
    
    for j=1:6
        eval(['load ' iFiles(iidx+j).name]);
        subplot(3,2,j); imshow(vData); text(100,100, num2str(j), 'Color', 'g');    
        drawnow;
    end
   
    a = input('Frame? ');
    
    if isempty(a)
        ampdataOne = [ampdataOne, zeros(1, Fs*10)];
        ampdataTwo = [ampdataTwo, zeros(1, Fs*10)];
    end
    
    if ~isempty(a) % The fish is in the correct position in these frames
        
        eval(['load ' eFiles(eidx).name]); % load the electrical data
                   
          ampdataOne = [ampdataOne, data((tim > 10*(a-1) & tim <= 10*a),1)'];
          ampdataTwo = [ampdataTwo, data((tim > 10*(a-1) & tim <= 10*a),2)'];
            
    end
    
%     ampdataOne(end+1) = mean(mean(vData(100:200, 500:600)));
%     ampdataTwo(end+1) = mean(mean(vData(100:200, 500:600)));
    ampdataOne(end+1) = mean(mean(vData(700:800, 500:600)));
    ampdataTwo(end+1) = mean(mean(vData(700:800, 500:600)));
    
    iidx = iidx+6;
    eidx = eidx+1;
    
end

