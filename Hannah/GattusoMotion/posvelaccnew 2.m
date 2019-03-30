function out = posvelaccnew(spikes, stim, keys, dur)
%Usage: out = postvelacc(spikes, stim, Fs)
Fs = 1/stim.interval; 
tim = 1/Fs:1/Fs:length(stim.values)/Fs;
buff = 0.400;

[b,a] = butter(3, 30/Fs, 'low'); 
[d,c] = butter(5, 20/Fs, 'low'); 


% Make stimuli

firstorder = filtfilt(b,a,diff(stim.values));
secondorder = filtfilt(d,c,diff(firstorder));

cpos = []; cvel = []; cacc = [];

% For each epoch...
for i = length(keys):-1:1;

    %    st_tt = find(tim > keys(i) & tim < keys(i) + dur);
    sp_tt = find(spikes > keys(i) & spikes < keys(i) + dur);

% Get values for each spike

prelenpos = length(cpos);
prelenvel = length(cvel);
prelenacc = length(cacc);

for ss = length(sp_tt):-1:1;
    
    tt = find(tim < spikes(sp_tt(ss)) & tim > spikes(sp_tt(ss)) - buff);
    cpos(ss+prelenpos) = mean(stim.values(tt));
    cvel(ss+prelenvel) = mean(firstorder(tt));
    cacc(ss+prelenacc) = mean(secondorder(tt));

    pv(ss+prelenvel,:) = [cpos(ss+prelenpos) cvel(ss+prelenvel)];
    av(ss+prelenacc,:) = [cacc(ss+prelenacc) cvel(ss+prelenvel)];
    
end;

end;


out.posvel = hist3(pv,[25 25]);

out.accvel = hist3(av,[25 25]);

out.pos = cpos;
out.vel = cvel;
out.acc = cacc;

cmaxPV = max(max(out.posvel));
cmaxAV = max(max(out.accvel));
thresh = 0.0001;
pp = find(abs(out.acc) < thresh); %eliminates peak at zero velocity and zero acceleration

figure; clf;
subplot(221); surf(out.posvel'); view(0,90); caxis([0 cmaxPV]); xlabel('Position'); ylabel('Velocity'); 
subplot(222); surf(out.accvel'); view(0,90); caxis([0 cmaxAV]); xlabel('Acceleration'); ylabel('Velocity');
colormap('HOT');
subplot(223); plot(out.pos(pp), out.vel(pp),'.'); xlabel('Position'); ylabel('Velocity');
subplot(224); plot(out.acc(pp), out.vel(pp), '.'); xlabel('Acceleration'); ylabel('Velocity');

figure;
velsteps = 0.0005;
accsteps = 0.0000015;
posedges = -5:0.1:5;
veledges = min(out.vel):velsteps:max(out.vel);
accedges = min(out.acc):accsteps:max(out.acc);
subplot(131); a = histcounts(out.pos, posedges); plot(posedges(1:end-1), a/sum(a), 'b'), 
subplot(132); b = histcounts(out.vel, veledges); plot(veledges(1:end-1), b/sum(b), 'b'),
subplot(133); c = histcounts(out.acc, accedges); plot(accedges(1:end-1), c/sum(c), 'b'), 


stimedges = -5:0.1:5;
foedges = min(firstorder):velsteps: max(firstorder);
soedges = min(secondorder):accsteps:max(secondorder);

subplot(131); hold on, d = histcounts(stim.values, stimedges); plot(stimedges(1:end-1), d/ (sum(d)), 'r'), xlabel('Position'); ylabel('Time Percentge'); legend('Spikes', 'Stimulus')
subplot(132); hold on, f = histcounts(firstorder, foedges); plot(foedges(1:end-1), f / (sum(f)), 'r'), xlabel('Velocity'); ylabel('Time Percentage'); legend('Spikes', 'Stimulus')
subplot(133); hold on, g = histcounts(secondorder, soedges); plot(soedges(1:end-1), g/ (sum(g)), 'r'), xlabel('Acceleration'); ylabel('Time Percentage'); legend('Spikes', 'Stimulus');
% axis([-5*10^(-5), 5*10^(-5), 0, 0.005])
xlim([-5*10^(-5), 5*10^(-5)]);




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

