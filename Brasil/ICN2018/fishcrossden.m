function out = fishcrossden(in, idx, idxs)
% Run fishdist or fishdistS first

figure(1); clf; hold on;

if nargin == 2
    idxs = 1:length(in(idx).pair);
end

for j = idxs

di = in(idx).pair(j).descartes - mean(in(idx).pair(j).descartes);
di = di / max(abs(di));                                        
fr = in(idx).pair(j).dF - mean(in(idx).pair(j).dF);
fr = fr / max(abs(fr));                    

Fs = median(diff(in(idx).pair(j).sharedtims));

[out(j).p, out(j).w] = cpsd(di(1:20:end), fr(1:20:end), [], [], [], 0.05*Fs); 
    figure(1); plot(out(j).w, real(out(j).p), '-*'); 
    figure(1+j); clf; hold on; 
    plot(in(idx).pair(j).sharedtims, di, '*'); 
    plot(in(idx).pair(j).sharedtims, fr, '*');

end

