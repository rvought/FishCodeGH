function [out] = squigglesquaggle(highfreq, lowfreq, stepsize)



Fs = 1000;

tim001 = 1/Fs:1/Fs:0.1;
tim002 = 1/Fs:1/Fs:0.2;
tim005 = 1/Fs:1/Fs:0.5;
tim010 = 1/Fs:1/Fs:1.0;
tim020 = 1/Fs:1/Fs:2.0;
tim040 = 1/Fs:1/Fs:4.0;
tim050 = 1/Fs:1/Fs:5.0;
tim100 = 1/Fs:1/Fs:10.0;
tim200 = 1/Fs:1/Fs:20.0;

%500 msec of zeros (for beginning and end)
zipper = zeros (1, length (tim010));
pause = ones(1, length(tim010));

% Initial offset and reset from/to zero
        cvinddn = 0:-1/length(tim010):-stepsize;
%         cvinddn = cvinddn(1:length(tim010));
        rcvinddn = 0:1/length(tim010): 1;
        cvendup = cvinddn(end:-1:1);
        retreat = -4: 1/length(tim005):0;
        finish = 4: -1/length(tim005):0;
        
% Top and bottom for CV stimuli        

    top = ones(1, length(tim002));
    bottom = -top;
    
%% slow sinusoid that changes location
squiggle = -1*sin(tim200*2*pi*lowfreq);
begin = [squiggle zipper];
forward = [];
for i = stepsize:stepsize:5;
    forward = [forward (cvinddn - (i-stepsize)) (squiggle - i) -1*pause*i];
end 
backward = [];
backsquiggle = 1*sin(tim200*2*pi*lowfreq);
for i = stepsize:stepsize:5;
    backward = [backward (-1*(cvinddn - (i-stepsize))) pause*i (backsquiggle + i)  pause*i];
end
slowstim(:,1) = [zipper begin forward retreat backward finish];
slowstim = (1/5)*slowstim;
stimtim = 1/Fs:1/Fs:length(slowstim)/Fs;
figure(1);
plot(stimtim,slowstim)

%% fast sinusoid that changes location 
fastsquiggle = sin(tim100*2*pi*highfreq);

fastbegin = [fastsquiggle fastsquiggle zipper];
fastforward = [];
for i = stepsize:stepsize:5;
    fastforward = [fastforward (cvinddn - (i-stepsize)) -1*pause*i (fastsquiggle - i) (fastsquiggle-i) -1*pause*i];
end 
fastbackward = [];
fastbacksquiggle = 1*sin(tim100*2*pi*highfreq);
for i = stepsize:stepsize:5;
    fastbackward = [fastbackward (-1*(cvinddn - (i-stepsize))) pause*i (fastbacksquiggle + i) (fastbacksquiggle + i) pause*i];
end
faststim(:,1) = [zipper fastbegin fastforward retreat fastbackward finish];
faststim = (1/5)*faststim;
faststimtim = 1/Fs:1/Fs:length(faststim)/Fs;
figure(2);
plot(faststimtim,faststim)


%% fast on slow sinusoid7
combo = sin(tim100*2*pi*lowfreq) + sin(tim100*2*pi*highfreq);
combo = (1/5)*combo;
combostim(:,1)=[zipper (1/5)*squiggle zipper (1/5)*fastsquiggle (1/5)*fastsquiggle zipper combo -1*combo zipper];
combotim = 1/Fs:1/Fs:length(combostim)/Fs;
figure(3);
plot(combotim, combostim)
% hold on
% plot(testtim, test, 'r')
% hold on
% plot(test2tim, test2, 'g')


%% tuberous 

%constant amplitude
tstim = (1/5)*ones(1, length(slowstim));
tstimfast = (1/5)*ones(1, length(faststim));
tstimcombo = (1/5)*ones(1, length(combostim));
%oscillating tuberous
AMFreq = 20;
ttim = 1/Fs:1/Fs:length(slowstim)/Fs;
ttimfast = 1/Fs:1/Fs:length(faststim)/Fs;
ttimcombo = 1/Fs:1/Fs:length(combostim)/Fs;
tstimosc = (1/40)* (sin(2*pi*ttim*AMFreq)) + .2;
tstimoscfast = (1/40)* (sin(2*pi*ttimfast*AMFreq)) + .2;
tstimosccombo = (1/40)* (sin(2*pi*ttimcombo*AMFreq)) +.2;
figure(4)
plot(tstimosc)
%% writing waveforms
% slowstim(:,2) = tstim;
% audiowrite('slowstimconstant.wav', slowstim, Fs);

faststim(:,2) = tstimfast;
audiowrite('faststimconstant.wav', faststim, Fs);

combostim(:,2) = tstimcombo;
audiowrite('combostimconstant.wav', combostim, Fs);


% slowstim(:,2) = tstimosc;
% audiowrite('slowstimoscillating.wav', slowstim, Fs);

faststim(:,2) = tstimoscfast;
audiowrite('faststimoscillating.wav', faststim, Fs);

combostim(:,2) = tstimosccombo;
audiowrite('combostimoscillating.wav', combostim, Fs);




end

