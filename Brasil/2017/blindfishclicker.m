function out = dblindfishclicker(in)
% Usage out = blindfishclicker(in)

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
hold on; for j=1:length(in); plot(in(j).tim, in(j).freq, 'o', 'MarkerEdgeColor', clrs(j,:),'MarkerFaceColor', clrs(j,:)); end;

figure(2); % Get user clicks

[eventimes, eventfreqs] = ginputc('Color','k','ShowPoints',true,'ConnectPoints',false);

    %[eventfreqs, eidx] = sort(eventfreqs);
    %eventimes = eventimes(eidx);
    
    % Find the next time stamp and the associate indices
    for jj = 1:length(eventimes)
        
        rango = 1;
        timrang = [(eventimes(jj) - rango) (eventimes(jj) + rango)];
        
        for zz = 1:length(in)
            freqs(zz) = mean(in(zz).freq(in(zz).tim > timrang(1) & in(zz).tim < timrang(2)));
        end               
        
    end
         
     
    % Get the nearest frequency to the click
    
%     for kk = 1:length(eventimes)
%        
%         
%         
%         diffs = abs(eventfreqs(kk) - in.freq 
%         
%         
%     end
%     
    
    
% Replot the clicks 
    for k = 1:length(eventfreqs)
        figure(2); plot(eventimes(k), eventfreqs(k), 'ko', 'MarkerSize', 10, 'MarkerEdgeColor','y','MarkerFaceColor','k')
        % figure(1); plot(eventimes(k), eventfreqs(k), 'ko', 'MarkerSize', 10, 'MarkerEdgeColor','y','MarkerFaceColor','k')
    end
    
    
    
% Assemble the output structure    
    out.freqs = eventfreqs; % The click frequencies
    out.tims = eventimes; % The click times
 
        
fprintf('Remember to save the data!\n');
