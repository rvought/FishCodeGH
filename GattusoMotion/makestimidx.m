% figure(1); clf;
% 
% Fs = AL(2).s(5).pFs;
% tim = 1/Fs:1/Fs:length(AL(2).s(5).pos)/Fs;
% 
% ax(1) = subplot(212); plot(tim, AL(2).s(5).pos)
% 
% 
% ys = [0 1];
% 
% ax(2) = subplot(211); hold on;
% for j=1:length(AL(2).s(5).st)
%     plot([AL(2).s(5).st(j) AL(2).s(5).st(j)], ys, 'k-');
% end
% linkaxes(ax, 'x');

% Ev(length(AL)).s(1) = AL(length(AL)).s(1);
% 
% for j = 1:length(AL)
% for k = 1:length(AL(j).s)        
%        if ~isnan(AL(j).s(k).pos)           
%            Ev(j).s(end+1) = AL(j).s(k);
%        end    
%         
% end
% end

% 
% fish = al;
% for k = 1:length(al)
%     holder = [];
%     for mm = 1:length(al(k).s)
%         if ~isnan(al(k).s(mm).tub)
%             holder = [holder mm]; 
%         end
%     end
%     if length(holder) > 0
%     for nn = 1:holder(1)
%     fish(k).s(nn).tub = al(k).s(holder(1)).tub; 
%     fish(k).s(nn).date = al(k).s(holder(1)).date;
%     fish(k).s(nn).filename = al(k).s(holder(1)).filename;
%     fish(k).s(nn).USID = al(k).s(holder(1)).USID;
%     end
%     for j = 2:length(holder)
%         for pp = 1:(holder(j)-holder(j-1))
%        fish(k).s((holder(j-1)+pp)).tub = al(k).s(holder(j)).tub;
%        fish(k).s((holder(j-1)+pp)).date = al(k).s(holder(j)).date;
%        fish(k).s((holder(j-1)+pp)).filename= al(k).s(holder(j)).filename;
%        fish(k).s((holder(j-1)+pp)).USID = al(k).s(holder(j)).USID;
%         end
%     end
%     end
% end
%    
fish = al;
howhi = [];
for k = 1:length(al)
    for j = 1:length(al(k).s)
        howhi(end+1) = max(xcorr(al(1).s(1).pos, al(k).s(j).pos));
    end
end
plot(howhi)
[~, y] = ginput(1);
for k = 1:length(al)
    for j = 1:length(al(k).s)
        howhi = max(xcorr(al(1).s(1).pos, al(k).s(j).pos));
        if howhi > y
            fish(k).s(j).USID = 10;
        end
    end
end