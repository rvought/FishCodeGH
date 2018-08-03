% function out = posvelacc(spikes, stim, keys, dur, Fs)
% %Usage: out = postvelacc(spikes, stim, Fs)
% 
% tim = 1/Fs:1/Fs:length(stim)/Fs;8u90
% buff = 0.100;
% 
% [b,a] = butter(3, 30/Fs, 'low'); 
% [d,c] = butter(5, 20/Fs, 'low'); 
% 
% 
% % Make stimuli
% 
% firstorder = filtfilt(b,a,diff(stim));
% secondorder = filtfilt(d,c,diff(firstorder));
% 
% cpos = []; cvel = []; cacc = [];
% 
% % For each epoch...
% for i = length(keys):-1:1;
% 
%     %    st_tt = find(tim > keys(i) & tim < keys(i) + dur);
%     sp_tt = find(spikes > keys(i) & spikes < keys(i) + dur);
% 
% % Get values for each spike
% 
% prelenpos = length(cpos);
% prelenvel = length(cvel);
% prelenacc = length(cacc);
% 
% for ss = length(sp_tt):-1:1;
%     
%     tt = find(tim < spikes(sp_tt(ss)) & tim > spikes(sp_tt(ss)) - buff);
%     cpos(ss+prelenpos) = mean(stim(tt));
%     cvel(ss+prelenvel) = mean(firstorder(tt));
%     cacc(ss+prelenacc) = mean(secondorder(tt));
% 
%     pv(ss+prelenvel,:) = [cpos(ss+prelenpos) cvel(ss+prelenvel)];
%     av(ss+prelenacc,:) = [cacc(ss+prelenacc) cvel(ss+prelenvel)];
%     
% end;
% 
% end;
% 
% 
% out.posvel = hist3(pv,[60 60]);
% out.accvel = hist3(av,[60 60]);
% 
% out.pos = cpos;
% out.vel = cvel;
% out.acc = cacc;
% 
% 
% %% Useful Commands
% 
% % Fs = 1/bg_2015_Ch2.interval; 
% % tim = 1/Fs:1/Fs:bg_2015_Ch2.length/Fs;
% % plot(tim, bg_2015_Ch2.values);
% %
% % oogabuug = posvelacc(bg_2015_Ch405.times, bg_2015_Ch2.values, bg_2015_Ch31.times, 20, Fs);
% % figure(1); subplot(121); surf(oogabuug.posvel'); view(0,90); subplot(122); surf(oogabuug.accvel'); view(0,90);
% % colormap('HOT');
% % caxis([0 10]);
% 
% % figure(2); subplot(121); plot(oogabuug.pos, oogabuug.vel,'*'); subplot(122); plot(oogabuug.acc, oogabuug.vel, '*');
% 
% 
% % save filename.mat oogabuug
% 
