function TE = fishtransferEntropy(dat, met, ind)

idx = (met.dataset == ind);
listofpairs = unique(met.pair(idx));

ll = 1;
kk = 15;
tt = 100;
ww = 1;

figure(2); clf; hold on;

TE{length(listofpairs)} = [];

for j = 1:length(listofpairs)
    
   pairidx = find(met.dataset == ind & met.pair == listofpairs(j));
   
   
   for k = 1:length(pairidx)
       
      [TE{j}(end+1),~,~] = transferEntropyPartition(dat.trial{pairidx(k)}(1,:), dat.trial{pairidx(k)}(2,:), ll, kk); 
      %TE{j}(end+1) = transferEntropyKDE(dat.trial{pairidx(k)}(1,:), dat.trial{pairidx(k)}(2,:), ll, kk, tt, ww); 
       
   end
    
   plot(TE{j});
   
end