function plotvelacc(stim, Fs)

tim = 1/Fs:1/Fs:length(stim)/Fs;

subplot(211); plot(tim, stim);

vel = diff(stim);

subplot(223);

plot(stim(2:end), vel);

ax = axis;

acc = diff(vel);


subplot(224);

plot(stim(2:end-1),acc);

axis([ax(1) ax(2) ax(3)/100 ax(4)/100]);

