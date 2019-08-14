
cavedFDistCC = [];
caveDistance = [];
cavedF = [];

figure(1); clf;

%Cave 3
[dFdist, orig] = preTE(cave(3));
    out = getCorr(dFdist, orig, 300, 10, [1 3]);
    
for j=1:length(out) 
    cavedFDistCC = [cavedFDistCC out.dfdistCC];
    caveDistance = [caveDistance out.meanDist];
    cavedF = [cavedF out.meandF];
end
    figure(1);
    subplot(121); plot(caveDistance, cavedFDistCC, 'b*'); hold on;
    subplot(122); plot(cavedF, cavedFDistCC, 'b*'); hold on;

% Cave 4
[dFdist, orig] = preTE(cave(4));
    out = getCorr(dFdist, orig, 300, 10, [1 4]);
    
    figure(1);
    subplot(121); plot(caveDistance, cavedFDistCC, 'b*'); 
    subplot(122); plot(cavedF, cavedFDistCC, 'b*'); 
        
for j=1:length(out) 
    cavedFDistCC = [cavedFDistCC out.dfdistCC];
    caveDistance = [caveDistance out.meanDist];
    cavedF = [cavedF out.meandF];
end

% Cave 5
[dFdist, orig] = preTE(cave(5));
    out = getCorr(dFdist, orig, 300, 10, [1 5]);
    
    figure(1);
    subplot(121); plot(caveDistance, cavedFDistCC, 'b*'); 
    subplot(122); plot(cavedF, cavedFDistCC, 'b*'); 
        
for j=1:length(out) 
    cavedFDistCC = [cavedFDistCC out.dfdistCC];
    caveDistance = [caveDistance out.meanDist];
    cavedF = [cavedF out.meandF];
end

% Cave 7
[dFdist, orig] = preTE(cave(7));
    out = getCorr(dFdist, orig, 300, 10, [1 7]);
    
    figure(1);
    subplot(121); plot(caveDistance, cavedFDistCC, 'b*'); 
    subplot(122); plot(cavedF, cavedFDistCC, 'b*'); 

    
for j=1:length(out) 
    cavedFDistCC = [cavedFDistCC out.dfdistCC];
    caveDistance = [caveDistance out.meanDist];
    cavedF = [cavedF out.meandF];
end





