function out = fishcrossden(in, idx, idxs)
% Run fishdist or fishdistS first
% e.g. in = fishdist(f, [100 1100])


figure(1); clf; hold on;

numfish = length(in(idx).fish);
for k=1:numfish
    out(j).fishnum(k).meandist = [];
    out(j).fishnum(k).meandF = [];
end

if nargin == 2
    idxs = 1:length(in(idx).pair);
end

for j = idxs

di = in(idx).pair(j).descartes - mean(in(idx).pair(j).descartes);
    out(j).fishnum(in(idx).pair(j).fishnums(1)).meandist(end+1) = mean(di);
    out(j).fishnum(in(idx).pair(j).fishnums(2)).meandist(end+1) = mean(di);    
di = di / max(abs(di));                                        

fr = in(idx).pair(j).dF - mean(in(idx).pair(j).dF);
    out(j).fishnum(in(idx).pair(j).fishnums(1)).meandF(end+1) = mean(di);
    out(j).fishnum(in(idx).pair(j).fishnums(2)).meandF(end+1) = mean(di);
fr = fr / max(abs(fr));                    

Fs = median(diff(in(idx).pair(j).sharedtims));

[out(j).p, out(j).w] = cpsd(di(1:20:end), fr(1:20:end), [], [], [], 0.05*Fs); 
    figure(1); plot(out(j).w, real(out(j).p), '-*'); 
    figure(1+j); clf; hold on; 
    plot(in(idx).pair(j).sharedtims, di, '*'); 
    plot(in(idx).pair(j).sharedtims, fr, '*');

end

