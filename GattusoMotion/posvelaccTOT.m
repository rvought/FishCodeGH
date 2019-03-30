function out = posvelaccTOT(spikechan, stimulus)
%Usage: out = postvelaccTOT(spikes, stim, Fs)

Fs = 1/stimulus.interval;
stim = stimulus.values;

spikes = spikechan.times;

tim = 1/Fs:1/Fs:length(stim)/Fs;

buff = 0.100; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[b,a] = butter(3, 30/Fs, 'low'); 
[d,c] = butter(5, 20/Fs, 'low'); 

% Make stimuli

firstorder = filtfilt(b,a,diff(stim));
secondorder = filtfilt(d,c,diff(firstorder));

cpos = []; cvel = []; cacc = [];


for ss = length(spikes):-1:1;
    
    tt = find(tim < spikes(ss) & tim > spikes(ss) - buff);
    cpos(ss) = mean(stim(tt));
    cvel(ss) = mean(firstorder(tt));
    cacc(ss) = mean(secondorder(tt));

    pv(ss,:) = [cpos(ss) cvel(ss)];
    av(ss,:) = [cacc(ss) cvel(ss)];
    
end;


%out.posvel = hist3(pv,{-6:0.01:6 -0.0012:.0024/128:0.0012});
out.posvel = hist3(pv,[20 20]);
out.accvel = hist3(av,[20 20]);

out.pos = cpos;
out.vel = cvel;
out.acc = cacc;

%surf(out.posvel');view(0,90); caxis([0 5]);
%xlim([0 120]);
%ylim([0 102]);
%% Useful Commands

% Fs = 1/bg_2015_Ch2.interval; 
% tim = 1/Fs:1/Fs:bg_2015_Ch2.length/Fs;
% plot(tim, bg_2015_Ch2.values);
%
% oogabuug = posvelacc(bg_2015_Ch405.times, bg_2015_Ch2.values, bg_2015_Ch31.times, 20, Fs);
% figure(1); subplot(121); surf(oogabuug.posvel'); view(0,90); subplot(122); surf(oogabuug.accvel'); view(0,90);
% colormap('HOT');
% caxis([0 10]);

% figure(2); subplot(121); plot(oogabuug.pos, oogabuug.vel,'*'); subplot(122); plot(oogabuug.acc, oogabuug.vel, '*');


% save filename.mat oogabuug

