function ISI = ISIfinder(neuron)
ISI = [];
size = input('size for ISI');
for jj = 1:length(neuron)
    if neuron(jj).sizeDX == size
    ISI = [ISI diff(neuron(jj).st)']; 
    end
end
binsize = input('binsize');
cutoff = input('ISI cutoff');
fredrico = 0:binsize:cutoff;
figure;
a = histogram(ISI, 'BinEdges', fredrico);
xlabel('ISI');
ylabel('Frequency');
end