function out = ccxd(spikes, stimulus, Fs, plotme)
% Usage out = ccxd(spikes, stimulus, Fs, plotme)
% spikes are the spike times
% Stimulus is the stimulus waveform
% Fs is the samplerate for the stimulus
% If you want to plot, include the plot window number for plotme
% Depends on fftmachine

%% Setup

stepreduction = 50; % An integer to reduce the sample rate to focus on low frequencies

stimulus = stimulus(1:stepreduction:end);
Fs = Fs/stepreduction;

out = fftmachine(stimulus, Fs, 3);

tim = 1/Fs:1/Fs:length(stimulus)/Fs;

spikedur = 2; % spike duration in milliseconds
    spikedur = spikedur / 1000;

%% ALPS
% The alps weren't built in a day. This is our alpha function for the spike
% times.
alps = zeros(1, length(stimulus));
    for jj = 1:length(spikes)-1;
    
        alps(tim >= spikes(jj) & tim < spikes(jj)+spikedur) = 1;    
    
    end;
[b,a] = butter(3,2*100/Fs, 'low');
alps = filtfilt(b,a,alps);


[out.cpPxy, out.cpF] = cpsd(stimulus,alps,hamming(2048),2000,2048,Fs);
[out.mscoPxy,out.mscoF] = mscohere(stimulus,alps,hamming(2048),2000,2048,Fs);

%% Plot
if nargin > 3;
    figure(plotme); clf;
    subplot(311); plot(out.fftfreq, out.fftdata); xlim([0 50]);
    title('Stimulus');
    subplot(312); plot(out.cpF, abs(out.cpPxy), '*-'); xlim([0 50]);
    title('Cross Spectrum Phase');
    subplot(313); plot(out.mscoF, out.mscoPxy, '*-'); xlim([0 50]); ylim([0 1]);
    title('Coherence');
end;

