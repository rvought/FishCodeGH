function p = compareDistributions(x1,x2)
% Compare distributions using the Kolmogorov-Smirnov test when sample
% sizes are potentially unequal

x1 = x1(:);
x2 = x2(:);

n1 = length(x1);
n2 = length(x2);

L = round(min(n1,n2)/2);
N = 10000;

h = zeros(N,1);
for n = 1:N
    r1 = datasample(x1,L,'Replace',false);
    r2 = datasample(x2,L,'Replace',false);
    
    h(n) = kstest2(r1,r2,'Alpha',0.01);   
end

p = 1 - (sum(h)/N);
