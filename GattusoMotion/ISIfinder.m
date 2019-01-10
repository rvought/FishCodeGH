function ISI = ISIfinder(neuron)
ISI = [];
size = input('size for ISI');
binsize = input('binsize');
cutoff = input('ISI cutoff');
for kk = 1:length(neuron)
    for jj = 1:length(neuron(kk).s)
        if neuron(kk).s(jj).sizeDX == size
        ISI = [ISI diff(neuron(kk).s(jj).st)']; %for later data remove "prime"
        end
    end
end
fredrico = 0:binsize:cutoff;
figure;
a = histogram(ISI, 'BinEdges', fredrico);
xlabel('ISI');
ylabel('Frequency');
end
