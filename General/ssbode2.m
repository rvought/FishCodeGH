function [mag, phase, G] = ssbode2(data,f)
% function [mag, phase, G] = ssbode(data,f)
% INPUT:
% data - iddata
% N - number of frequency peaks
% OUTPUT:
% mag - gain in dB
% phase - phase in rad
% f - same as input
% s - T.F

datf = fft(data);

%%
U = datf.u;
Y = datf.y;
magU = abs(U);

%% Find the index of input frequencies

N = length(f);
for k = 1:N
    [~,peak_inds(k)] = min(abs(datf.Freq-f(k)));
end


%% Calculate parameters

mag = abs(Y(peak_inds))./abs(U(peak_inds));
phase = (angle(Y(peak_inds))-angle(U(peak_inds)));
offset = round(phase(1)/(2*pi))*2*pi;
phase = unwrap(phase-offset);
G = Y(peak_inds)./U(peak_inds);