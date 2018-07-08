function [mag, phase, f, s] = ssbode(data,N)
% function [mag, phase, f, s] = ssbode(data,N)
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

% Identify peaks of input
peak_inds = localmax(magU);
% Only use the 'real' peaks
peak_inds = peak_inds(peak_inds .* magU(peak_inds)>4);

%[s,idx] = sort(peak_inds,2,'descend');

% Sort
[~,idx] = sort(magU(peak_inds),1,'descend');
peak_inds = peak_inds(idx);

% Need only the strongest N peaks
if N <= length(peak_inds)
    peak_inds = peak_inds(1:N);
else
    error('Not enough peaks');
end
%%
% if N~=length(peak_inds)
%     N
%     length(peak_inds)
% figure(777)
% plot(datf.Frequency, magU, 'r', datf.Frequency(peak_inds),magU(peak_inds),'ko'); pause;
% xlim([0 26])
% figure(666)
% end

% peak_inds = sort(peak_inds(idx(1:N)));
% if 0
%     clf
%     plot(datf)
%     hold on
%     plot(datf.Frequency(peak_inds),magU(peak_inds),'ko');
%     xlim([0 20])
% end

%peak_inds = peak_inds(2);

mag = abs(Y(peak_inds))./abs(U(peak_inds));
phase = (angle(Y(peak_inds))-angle(U(peak_inds)));
offset = round(phase(1)/(2*pi))*2*pi;
phase = unwrap(phase-offset);
s = Y(peak_inds)./U(peak_inds);
f = datf.Frequency(peak_inds)/(2*pi);