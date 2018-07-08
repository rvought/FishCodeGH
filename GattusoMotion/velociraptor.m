function out = velociraptor(in, idx)
% Usage: out = velociraptor(sz, find(typeNsize(:,2)==#));

Fs = in(1).pFs;

stim =[];
tim = [];
spikes = [];

for j = 1:length(idx)
    stim = [stim in(idx(j)).pos'];
    if j == 1 
        nextim = 0; 
    else 
        nextim = tim(end); 
    end
    
    curtim = (1/Fs:1/Fs:length(in(idx(j)).pos)/Fs) + nextim;
    tim = [tim curtim];
    
    spikes = [spikes, (in(j).spiketimes' + nextim)];
    
end

buff = 0.100; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
[b,a] = butter(3, 30/Fs, 'low'); 
[d,c] = butter(5, 20/Fs, 'low'); 

% Make stimuli

firstorder = filtfilt(b,a,diff(stim));
secondorder = filtfilt(d,c,diff(firstorder));
%firstorder = diff(stim);
%secondorder = diff(firstorder);

cpos = []; cvel = []; cacc = [];


for ss = length(spikes):-1:1;    
    tt = find(tim(1:end-2) < spikes(ss) & tim(1:end-2) > spikes(ss) - buff);
    cpos(ss) = mean(stim(tt));
    cvel(ss) = mean(firstorder(tt)); %does not like when you select an even number of options*****
    cacc(ss) = mean(secondorder(tt));

    pv(ss,:) = [cpos(ss) cvel(ss)];
    av(ss,:) = [cacc(ss) cvel(ss)];
    
end;


thresh = 0.0001;
goodpoints = find(abs(av(:,1)) < thresh); 
pv = pv(goodpoints,:);
av = av(goodpoints,:);

%out.posvel = hist3(pv,{-6:0.01:6 -0.0012:.0024/128:0.0012});
out.posvel = hist3(pv,[50 50]);
    out.tmpPV = medfilt2(out.posvel, [3 3]); cmaxPV = max(max(out.tmpPV)); 
out.accvel = hist3(av,[50 50]);
    out.tmpAV = medfilt2(out.accvel, [3 3]); cmaxAV =  max(max(out.tmpAV));

out.pos = cpos;
out.vel = cvel;
out.acc = cacc;

% surf(out.posvel');view(0,90); caxis([0 5]);
% xlim([0 120]);
% ylim([0 102]);
%% Useful Commands

% Fs = 1/bg_2015_Ch2.interval; 
% tim = 1/Fs:1/Fs:bg_2015_Ch2.length/Fs;
% plot(tim, bg_2015_Ch2.values);
%
% oogabuug = posvelacc(bg_2015_Ch405.times, bg_2015_Ch2.values, bg_2015_Ch31.times, 20, Fs);

pp = find(abs(out.acc) < thresh); 

%figure(1); clf;
figure; clf;
subplot(221); surf(out.posvel', 'EdgeColor', 'none'); view(0,90); caxis([0 cmaxPV]); xlabel('Position'); ylabel('Velocity'); 
subplot(222); surf(out.accvel', 'EdgeColor', 'none'); view(0,90); caxis([0 cmaxAV]); xlabel('Acceleration'); ylabel('Velocity');
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

subplot(131); hold on, d = histcounts(stim, stimedges); plot(stimedges(1:end-1), d/ (sum(d)), 'r'), xlabel('Position'); ylabel('Time Percentge'); legend('Spikes', 'Stimulus')
subplot(132); hold on, f = histcounts(firstorder, foedges); plot(foedges(1:end-1), f / (sum(f)), 'r'), xlabel('Velocity'); ylabel('Time Percentage'); legend('Spikes', 'Stimulus')
subplot(133); hold on, g = histcounts(secondorder, soedges); plot(soedges(1:end-1), g/ (sum(g)), 'r'), 
xlabel('Acceleration'); ylabel('Time Percentage'); legend('Spikes', 'Stimulus'); % axis([-5*10^(-5), 5*10^(-5), 0, 0.005])


% length(out.pos)
% length(stim)



% figure; clf;
% subplot(311); plot(posratio, r, '*')
% subplot(312); plot(velratio, b, '*')
% subplot(313); plot(accratio, g, '*')





% save filename.mat oogabuug

