function out = dblindfishclicker(in)
% Usage out = blindfishclicker(in)

pad = 10;

numcolors = 50;
oclrs = hsv(numcolors); % 20 colors for plotting up to 20 different fishes
clrs = zeros(numcolors, 3);
% figure(3); clf; hold on; for kk=1:20; plot(kk,kk, '*', 'MarkerEdgeColor', oclrs(kk,:)); end;
    %shuff = [11 2 7 13 18 1 16 12 10 17 3 8 14 19 20 4 16 5 15 9];
    shuff = randperm(numcolors);
    for i=1:numcolors; clrs(i,:) = oclrs(shuff(i),:); end; 


figure(1); clf; % XY Plot
hold on; for j=1:length(in); plot(in(j).tx, in(j).ty, '*', 'MarkerEdgeColor', clrs(j,:)); end;

figure(2); clf; % Fake Spectrogram plot
hold on; 
    for j=1:length(in) 
        plot(in(j).tim, in(j).freq, 'o', 'MarkerEdgeColor', clrs(j,:),'MarkerFaceColor', clrs(j,:)); 
        meanfreq(j) = mean(in(j).freq);
        plot(10, meanfreq(j), 'k*');
    end;

figure(2); % Get user clicks

[eventimes, eventfreqs] = ginputc('Color','k','ShowPoints',true,'ConnectPoints',false);

% Cut out the segments
    for j = 2:2:length(eventimes)
        
        ClickFreq = mean([eventfreqs(j-1), eventfreqs(j)]);
        
        diffreq = meanfreq - ClickFreq;
        
        [~, minidx] = min(abs(diffreq)); % This is the idx for the matching frequency
                        
        out(j/2).tim = in(minidx).tim(find(in(minidx).tim > eventimes(j-1)-pad & in(minidx).tim < eventimes(j)+pad));
        out(j/2).freq = in(minidx).freq(find(in(minidx).tim > eventimes(j-1)-pad & in(minidx).tim < eventimes(j)+pad));
        out(j/2).meanfreq = meanfreq(minidx);
        out(j/2).pad = pad;
        out(j/2).clicktimes = eventimes;
        out(j/2).clickfreqs = eventfreqs;
        out(j/2).freqidx = minidx;
    end    
 
figure(27); clf; hold on; 
    for j=1:length(out) 
        plot(out(j).tim, out(j).freq); 
        text(out(j).tim(1), out(j).freq(1), num2str(j));
    end;
    ylim([200 550]);
    
fprintf('Remember to save the data!\n');
