function out = squeezestim(spikes, stim, Fs)
% Usage out = sequeezestim(spikes, stim, Fs);

newstim = stim(1:10:end);

z = zeros(1,length(newstim));

tt = find(abs(newstim) < 0.005);

z(tt) = 1;

starts = find(diff(z) == 1);
ends = find(diff(z) == -1);


plot(stim); hold on; plot(stim(starts), 'go'); plot(stim(ends), 'ro'):
