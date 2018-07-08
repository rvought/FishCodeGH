function [ampullary tuberous combo] = amptub(varargin)
% sinnoise(PARAMS ...);  % Generates wav file signal with random amp mods

sf = 10000;          % sampling frequency (Hz)
nbits = 16;          % signifies 16 bit wave, don't change
seed = 234324;       % seed for random number generator
dees = 272727;       % seed for the second noisy signal
modfreq = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDIT THESE

%%%%%%%%%%%%%%% S1 frequency
freq = 300;          % Frequency of carrier sine wave (Hz)

%%%%%%%%%%%%%%% Noise 
lofnoise = 0.1;      % Low Cutoff frequency for band-pass noise (Hz)
hifnoise = 10;       % High Cutoff frequency for band-pass noise (Hz)

%%%%%%%%%%%%%%% Relative Amplitudes 
A = .9;             % Amplitude (fraction of maxumum = 1)
amp_weight = 0.01;   % Relative size of ampulary signal
dam = 20;           % depth of modulation (percent)

%%%%%%%%%%%%%%% Timing
tmax = 20;          % length of wav file (seconds)
tshift = 5;         % time shift for second noise (ampulary)

%%%%%%%%%%%%%%% Filename(s)
combofile = 's300lh01_10t09a01dam20_ts5';   % output filename (string)
tubfile = 'tubout';             % envelope of tuberous
shiftfile = 'ampout';           % shifted noise (ampulary)

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h = 1/sf;
Wn(1) = lofnoise*2/sf;
Wn(2) = hifnoise*2/sf;
[b,a] = butter(3,Wn(1),'high');
[d,c] = butter(3,Wn(2),'low');
time = h:h:tmax;
lt = length(time);

% Noise that will be the basis for our signals.
	randn('state',seed);
	rseq = randn(1,lt);

% Noise that will be added to uncohere our signals
	randn('state',dees);
	qesr = randn(1,lt);

% Low pass filter
noise = filter(d,c,rseq);
esion = filter(d,c,qesr);

% High pass filter
noise = filter(b,a,noise);
esion = filter(b,a,esion);

of = std(noise)*sqrt(2)*200/dam - max(noise);
fo = std(esion)*sqrt(2)*200/dam - max(esion);

% These are our basic noisy signals
noise = (noise + of)/(max(noise) + of);
esion = (esion + fo)/(max(esion) + fo);

% Construct the noise that goes in and out of coherence 

sinCOH    = (sin(2*pi*modfreq*time)+1)/2;
cosINCOH  = (cos(2*pi*modfreq*time)+1)/2;

cnoise = noise - mean(noise);
nnoise = esion - mean(esion);

cnoise = cnoise .* sinCOH;
nnoise = nnoise .* cosINCOH;

cohincoh = cnoise + nnoise; % this goes above and below zero.  

cohincoh = cohincoh - min(cohincoh); % this is positive only.

%figure; hold on; plot(1+qwer*2,'b'); plot(1+cohincoh*2,'r'); plot(cohincoh-qwer,'g');

% Take our Noise Signal and make the AMPULLARY
rangenoise = max(cohincoh) - min(cohincoh);
savenoise = 1.98*cohincoh/(rangenoise);
savenoise = savenoise + 0.99 - max(savenoise);
ampullary = savenoise;

% Make our TUBEROUS signal
y = A*sin(2*pi*freq*time).*noise;
tuberous = y;

% Construct our signals to add together.

siglen = length(ampullary);
nulsig = zeros(1,siglen);

% Make combined Tuberous signal
tubsig = [ tuberous tuberous nulsig ];

% Make combined Ampullary signal
ampsig = [ nulsig ampullary ampullary ];

wampsig = ampsig*amp_weight;
wtubsig = tubsig*(1-amp_weight);

combo = (wampsig + wtubsig)/2;

% wavwrite(combo,sf,nbits,combofile);

