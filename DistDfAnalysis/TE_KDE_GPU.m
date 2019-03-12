function TE = TE_KDE_GPU(X,Y,Xrange,Yrange,t,N)

TE = zeros(size(t));
for i = 1:length(t)
    TE(i) = transferEntropyKDE_range(X,Y,Xrange,Yrange,t(i),1,N,1); 
end