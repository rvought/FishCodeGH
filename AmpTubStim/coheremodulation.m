clear all
clf
fs = 10000;      % 10KHz sampling
Ts = 1/fs;       % 10KHz sampling
fcuttoff = 10;   % 50Hz
Nduration = 2^16;
v = randn(Nduration,1);
V = fft(v);

[b,a]=butter(2,fcuttoff/fs,'low');
%vfilt = filtfilt(b,a,v);

G = tf(b,a,-1);
fir = impulse(G,[1:Nduration]);
vfilt = fftfilt(fir,v);

FIR = fft(fir);
VFILT = FIR.*V;
vfilt = ifft(VFILT);

t = [0:Nduration-1]';

plot(t,[vfilt+2])


t1 = t;
t2 = [0:0.8:Nduration-1]';

v1 = interp1(t,vfilt,t1);
v2 = interp1(t,vfilt,t2);

v1input = repmat(v1,40,1);
v2input = repmat(v2,36,1);

m = min(length(v1input),length(v2input));
v1input = v1input(1:m);
v2input = v2input(1:m);

xc = xcorr(v1input,v2input);

%%
clear coh 
for i = 1:36;
    kstart = round((i-1) * m / 36 + 1);
    kend = kstart + round(m/36) -1;
    a=v1input(kstart:kend);
    b=v2input(kstart:kend);
    coh(i) = sum(mscohere(a,b));
end

plot(coh);
    


%clf
%plot([v2; v2])

