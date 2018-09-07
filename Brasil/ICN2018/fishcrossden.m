function out = fishcrossden(in, idx)


figure(1); clf; hold on;

for j=1:length(in(idx).pair)

di = in(idx).pair(j).descartes - mean(in(idx).pair(j).descartes);
di = di / max(abs(di));                                        
fr = in(idx).pair(j).dF - mean(in(idx).pair(j).dF);
fr = fr / max(abs(fr));                    

Fs = median(diff(out(idx).pair(j).sharedtims));

[p, w] = cpsd(di(1:10:end), fr(1:10:end), [], [], [], 0.1*Fs); 
    figure(1); plot(w, real(p), '-*'); 
    figure(1+j); clf; hold on; 
    plot(out(idx).pair(j).sharedtims, di, '*'); 
    plot(out(idx).pair(j).sharedtims, fr, '*');

end

>> [p, w] = cpsd(di(1:10:end), -fr(1:10:end), [], [], [], 0.1/median(diff(out(3).pair(1).sharedtims))); figure(2); clf; plot(w,real(p), '-*');
>> [p, w] = cpsd(di(1:10:end), -fr(1:10:end), [], [], [], 0.1/median(diff(out(3).pair(1).sharedtims))); figure(2); clf; plot(w,real(p), '-*');
>> 
>> figure(1); clf; plot(out(3).pair(1).descartes, di, '*'); hold on; plot(out(3).pair(1).descartes, fr, '*');
>> figure(1); clf; plot(out(3).pair(1).sharedtims, di, '*'); hold on; plot(out(3).pair(1).sharedtims, fr, '*');
>> [p, w] = cpsd(di(1:10:end), fr(1:10:end), [], [], [], 0.1/median(diff(out(3).pair(1).sharedtims))); figure(2); clf; plot(w,real(p), '-*');
