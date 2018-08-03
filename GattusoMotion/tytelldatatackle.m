function out = tytelldatatackle(st, tub, pos, keys) 

tFs = 1/tub.interval;           %establish time
pFs = 1/pos.interval;

ttim = 1/tFs:1/tFs:tub.length/tFs;
ptim = 1/pFs:1/pFs:pos.length/pFs;

firstorder = diff(pos.values);
secondorder = diff(firstorder);

for ss = length(st):-1:1;    
    tt = find(ptim(1:end-2) < st(ss) & ptim(1:end-2) > st(ss) - buff);
    cpos(ss) = mean(stim(tt));
    cvel(ss) = mean(firstorder(tt)); %does not like when you select an even number of options*****
    cacc(ss) = mean(secondorder(tt));

    pv(ss,:) = [cpos(ss) cvel(ss)];
    av(ss,:) = [cacc(ss) cvel(ss)];
    
end;

thresh = 0.0001;
goodpoints = find(abs(av(:,1)) < thresh); 
pv = pv(goodpoints,:);
av = av(goodpoints,:);

out.posvel = hist3(pv,[50 50]);
    out.tmpPV = medfilt2(out.posvel, [3 3]); cmaxPV = max(max(out.tmpPV)); 
out.accvel = hist3(av,[50 50]);
    out.tmpAV = medfilt2(out.accvel, [3 3]); cmaxAV =  max(max(out.tmpAV));
