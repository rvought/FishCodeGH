function out = SumSines( freqAmp, freqTub, amplitudes, Fs, dur, filename1, filename2, filename3, filename4)
% out = SumSines(freqAmp, freqTub, amplitudes, Fs, dur, filename);
%KEY: make sure to make one of the frequences common between amp and tub
%and put that at the end of your string. Make all the other numbers golden
%ratio (multiply by (1 + sqrt(5))/2  ) 
%freqAmp, pick five frequencies for the amp ex [1 6.796 10.6 29.448 12.5]
%freqTub, pick five frequencies for the tub ex  [1.62 4.2 17.151 18.2 12.5]
%amplitudes in volts [ampullary tuberous]. Ampullary should be less than 0.1 V (probably much
%less). Tuberous we need to work out - since that is modulation of the S1.
%Fs: ex 10000
%dur: duration of stimulis in seconds
%filename1: ('filename.wav') for regular Ampullary and... [tub both amp]
%filename2: for inverted ampullary.
%filename3: for [amp both tub]
%filename4: for inverted amp in [amp both tub]
%
%
%once done running, make sure to save out: ex) save 'SS_12.5.mat out'
%inside out is each individual sine wave used to generate out.amp and
%out.tub, which are found in: out.ampsingle and out.tubsingle
%to get these sine waves, type out.ampsingle(1,:)

scalr = 0.2; % This is "one Volt" when imported into Spike2
amplitudes = amplitudes * scalr;

sampleFreq = Fs;
nSamples = dur*Fs;
tim = (1:nSamples)/sampleFreq; 
sstub = zeros(1,length(tim));
ssamp = zeros(1,length(tim));

%making zeros/ones where only playing tub or amp
noTub = scalr * ones(1,length(tim));
noAmp = zeros(1,length(tim));

% End business
O = scalr * ones(1,floor(Fs*0.1));
Z = zeros(1,floor(Fs*0.1));

%Adding the Ampullary sines together, randomly picking if they should be
%neg or not
for k = length(freqAmp)-1:-1:1;
    
    %making 'random' either 1 or -1
    random = randi(2);
    if random == 2
        random = 1;
    else
        random = -1;
    end
        
    a(k,:) = random*sin(2*pi*freqAmp(k)*tim) * amplitudes(1);
    ssamp = ssamp + a(k,:);
    
    %to save the sine waves that were used to make out.amp
    out.ampsingle(k,:) = [noAmp a(k,:) a(k,:) Z];   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end;

%keeping the freq that is the same between Amp and Tub positive
last = length(freqAmp);
a(last,:) = sin(2*pi*freqAmp(last)*tim) * amplitudes(1);
ssamp = ssamp + a(last,:);
%to save the sine waves that were used to make out.amp
out.ampsingle(last,:)= [noAmp a(last,:) a(last,:) Z];  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
%Adding the Tuberous sines together
for k = length(freqTub):-1:1;
    t(k,:) = sin(2*pi*freqTub(k)*tim) * amplitudes(2);
    sstub = sstub + t(k,:);
    %to save the sine waves that were used to make out.tub
    out.tubsingle(k,:) = [noTub t(k,:) t(k,:) O];   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end;


    
%saving them
out.amp = ssamp;
out.tub = sstub;
out.Fs = Fs;
out.tim = tim;

%out.amp = out.amp * amplitudes(1);
out.tub = scalr + out.tub;

if max(abs(out.tub)) > 1; error('Lower the voltage'); end;

%making the sines end at 0 and 1 
out.amp = [Z Z Z Z Z Z Z Z Z Z Z Z Z Z Z Z Z Z noAmp out.amp out.amp Z Z Z Z Z Z Z Z Z Z Z Z Z Z Z Z Z Z];
out.tub = [O O O O O O O O O O O O O O O O O O out.tub out.tub noTub O O O O O O O O O O O O O O O O O O];

out.ampback = [Z Z Z Z Z Z Z Z Z Z Z Z Z Z Z Z Z Z out.amp out.amp noAmp Z Z Z Z Z Z Z Z Z Z Z Z Z Z Z Z Z Z];
out.tubback = [O O O O O O O O O O O O O O O O O O noTub out.tub out.tub O O O O O O O O O O O O O O O O O O];

%plot
figure(1)
newtim = 1/Fs:1/Fs:length(out.amp)/Fs;
    clf;
    plot(newtim,out.amp, 'g')
    hold on
    plot(newtim,out.tub,'r')

wavwrite(transpose([out.amp; out.tub]), Fs, 16, filename1);
wavwrite(transpose([-out.amp; out.tub]), Fs, 16, filename2);
wavwrite(transpose([out.ampback; out.tubback]), Fs, 16, filename3);
wavwrite(transpose([-out.ampback; out.tubback]), Fs, 16, filename4);
end


