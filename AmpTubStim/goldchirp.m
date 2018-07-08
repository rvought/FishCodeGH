function out = goldchirp(centerfreq, amplitudes, Fs, dur, filename)
% out = goldchirp(centerfreq, amplitudes, Fs, dur, filename)
%centerfreq: ex (1,2,4,8, or 16)
%amplitudes: ex [0.1 0.2]
%Fs: ex 10000
%dur: duration of stimulis in seconds (32, 16, 8, 4, or 2)
%filename: ('filename.wav')

scalr = 0.2; % This is "one Volt" when imported into Spike2
amplitudes = amplitudes * scalr;

goldenRatio = (1+sqrt(5))/2;

chirp1Freq = centerfreq;
chirp2Freq = centerfreq*goldenRatio;

sampleFreq = Fs;
nSamples = dur*Fs;
tim = (1:nSamples)/sampleFreq; 

out.amp = sin(2*pi*chirp1Freq*tim)* amplitudes(1);
out.tub = sin(2*pi*chirp2Freq*tim)* amplitudes(2);

plot(out.amp,out.tub);
figure(2)
hold off
plot(tim,out.amp)
hold on
plot(tim,out.tub,'r')

O = scalr * ones(1,floor(Fs*0.1));
Z = zeros(1,floor(Fs*0.1));

out.Fs = Fs;
out.tim = tim;

damp = diff(diff(out.amp));
dtub = diff(diff(out.tub));

out.tub = scalr + (out.tub); %0.2 = 1 when importing into Spike2

%out.amp = [Z Z Z Z Z Z Z Z Z Z Z Z Z Z Z Z Z Z out.amp Z Z Z Z Z Z Z out.amp Z Z Z Z Z Z Z out.amp Z Z Z Z Z Z Z Z Z Z Z Z Z Z Z Z Z Z]; %for CF = 32
%out.tub = [O O O O O O O O O O O O O O O O O O out.tub O O O O O O O out.tub O O O O O O O out.tub O O O O O O O O O O O O O O O O O O]; %for CF = 32

out.amp = [Z Z Z Z Z Z Z Z Z Z Z Z Z Z Z Z Z Z out.amp Z Z Z Z Z Z Z Z Z Z Z Z Z Z Z Z Z Z];
out.tub = [O O O O O O O O O O O O O O O O O O out.tub O O O O O O O O O O O O O O O O O O];

figure(3); plot(damp,dtub, 'm');
% max(abs(out.tub))
% wavwrite(transpose([out.amp; out.tub]), Fs, 16, filename);
end