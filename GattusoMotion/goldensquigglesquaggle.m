function [out] = goldensquigglesquaggle(lowfreq, stepsize)



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
tim400 = 1/Fs:1/Fs:40.0;
tim800 = 1/Fs:1/Fs:80.0;

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
    highfreq = (2*(1+sqrt(5))/2)*lowfreq; %multiplies lower frequency by the golden ratio and an integer
    
%% slow sinusoid that changes location
squiggle = (-1/2)*sin(tim800*2*pi*lowfreq);
begin = [squiggle zipper];
forward = [];
for i = stepsize:stepsize:5;
    forward = [forward (cvinddn - (i-stepsize)) (squiggle - i) -1*pause*i];
end 
backward = [];
backsquiggle = (1/2)*sin(tim800*2*pi*lowfreq);
for i = stepsize:stepsize:5;
    backward = [backward (-1*(cvinddn - (i-stepsize))) pause*i (backsquiggle + i)  pause*i];
end
slowstim(:,1) = [zipper begin forward retreat backward finish];
slowstim = (1/5)*slowstim;
stimtim = 1/Fs:1/Fs:length(slowstim)/Fs;
% figure(1);
% plot(stimtim,slowstim)

%% fast sinusoid that changes location 
fastsquiggle = sin(tim400*2*pi*highfreq);

fastbegin = [fastsquiggle zipper];
fastforward = [];
for i = stepsize:stepsize:5;
    fastforward = [fastforward (cvinddn - (i-stepsize)) -1*pause*i (fastsquiggle - i) -1*pause*i];
end 
fastbackward = [];
fastbacksquiggle = 1*sin(tim400*2*pi*highfreq);
for i = stepsize:stepsize:5;
    fastbackward = [fastbackward (-1*(cvinddn - (i-stepsize))) pause*i (fastbacksquiggle + i) pause*i];
end
faststim(:,1) = [zipper fastbegin fastforward retreat fastbackward finish];
faststim = (1/5)*faststim;
faststimtim = 1/Fs:1/Fs:length(faststim)/Fs;
% figure(2);
% plot(faststimtim,faststim)


%% fast on slow sinusoid
combo = (1/2)*sin(tim800*2*pi*lowfreq) + sin(tim800*2*pi*highfreq);
combo = (1/5)*combo;
combostim(:,1)=[zipper (1/5)*squiggle (1/5)*squiggle zipper (1/5)*fastsquiggle (1/5)*fastsquiggle  zipper combo combo zipper];
combotim = 1/Fs:1/Fs:length(combostim)/Fs;
% figure(3);
% plot(combotim, combostim)
% hold on
% plot(testtim, test, 'r')
% hold on
% plot(test2tim, test2, 'g')

%plot to test if covers full range
combovel = diff(combo);
comboacc = diff(combovel);
combovel(:,end) = [];
% figure(5)
% plot(combovel, comboacc)

%% finding the receptive field (written with 1's because it will scale up when put into spike2)
rampdown = 0:-1/length(tim005):-5;
rampup = 5:-1/length(tim005):0; 
% freqprofile = [1 2 4 8 16 32]; 
% freqtim = 1/Fs: 1/Fs:(1/6)*length(go)/Fs;
freqtim = 1/Fs: 1/Fs: length(tim010)/Fs;
oscstim1 = sin(2*pi*freqtim*1);
oscstim2 = sin(2*pi*freqtim*2);
oscstim3 = sin(2*pi*freqtim*4);
oscstim4 = sin(2*pi*freqtim*8);
oscstim5 = sin(2*pi*freqtim*16);
oscstim6 = sin(2*pi*freqtim*32);
figure(8)
oscstim = (1/5)*[oscstim1 oscstim1 oscstim2 oscstim3 oscstim4 oscstim5 oscstim6];
go = ones(1,length(oscstim));
plot(oscstim)
workbb=[];
for i= 0:1:9;
    move = -5+i:1/length(tim002):-4+i;
    workbb = [workbb go*(-5+i) move];
end 
halt = 5*[go];
fieldtest = [rampdown workbb halt rampup];
fieldtest = 0.2*fieldtest;
figure(6)
plot(fieldtest)
silence1 = zeros(1, length(rampdown));
silence2 = zeros(1, length(move));
silence3 = zeros(1, length(rampup));
testtub = [silence1 oscstim silence2 oscstim silence2 oscstim silence2 oscstim silence2 oscstim silence2 oscstim silence2 oscstim silence2 oscstim silence2 oscstim silence2 oscstim silence2 oscstim silence3];
figure(6)
hold on
plot(testtub, 'r')
receptivestim(:,1) = fieldtest;
receptivestim(:,2) = testtub;

% backward = [];
% backsquiggle = (1/2)*sin(tim800*2*pi*lowfreq);
% for i = stepsize:stepsize:5;
%     backward = [backward (-1*(cvinddn - (i-stepsize))) pause*i (backsquiggle + i)  pause*i];
% end

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
% figure(4)
% plot(tstimosc)
%% writing waveforms
slowstim(:,2) = tstim;
audiowrite('slowstimconstant.wav', slowstim, Fs);

faststim(:,2) = tstimfast;
audiowrite('faststimconstant.wav', faststim, Fs);

combostim(:,2) = tstimcombo;
audiowrite('combostimconstant.wav', combostim, Fs);


slowstim(:,2) = tstimosc;
audiowrite('slowstimoscillating.wav', slowstim, Fs);

faststim(:,2) = tstimoscfast;
audiowrite('faststimoscillating.wav', faststim, Fs);

combostim(:,2) = tstimosccombo;
audiowrite('combostimoscillating.wav', combostim, Fs);

audiowrite('receptivestim.wav', receptivestim, Fs);





end

