function out = phaseplot(spiketimes, stimtims, bins)
% out = phaseplot(spiketimes, stimtims, bins)

foo = [];

for i = 2:length(stimtims);
    out.trial{i} = spiketimes(find(spiketimes > stimtims(i-1) & spiketimes < stimtims(i))) - stimtims(i-1);
    out.dur(i) = stimtims(i) - stimtims(i-1);
    
    out.deg{i} = 360 * (out.trial{i} ./ out.dur(i));
    foo = [foo out.deg{i}'];
end;

hist(foo,bins); xlim([0 360]);