ss = [];
for k = 18:23
    for j=1:length(outs(k).pair)
        if ~isempty(outs(k).pair(j).dFmean)
            ss(end+1) = outs(k).pair(j).dFmean;
        end
    end
end

cc = [];
for k = 3:17
    for j=1:length(outc(k).pair)
        if ~isempty(outc(k).pair(j).dFmean)
            cc(end+1) = outc(k).pair(j).dFmean;
        end
    end
end

ssd = [];
for k = 18:23
    for j=1:length(outs(k).pair)
        if ~isempty(outs(k).pair(j).meanDist)
            ssd(end+1) = outs(k).pair(j).meanDist;
        end
    end
end

ccd = [];
for k = 3:17
    for j=1:length(outc(k).pair)
        if ~isempty(outc(k).pair(j).meanDist)
            ccd(end+1) = outc(k).pair(j).meanDist;
        end
    end
end

%ssd = ssd(ssd < 250);
%ccd = ccd(ccd < 250);

% Plot dF distribution
figure; clf; 
aa = 0:5:250;              
    ax(1) = subplot(211); hist(cc, aa);
    ax(2) = subplot(212); hist(ss, aa);
    linkaxes(ax, 'x');
xlim([0 250])

figure; clf;
aa = logspace(0.1,3, 50);
cch = hist(cc, aa);
ssh = hist(ss, aa);
clf; semilogx(aa,cch/max(cch), '*-');
hold on; semilogx(aa, ssh/max(ssh), '*-');


% Plot distance distribution
% aa = 0:5:250;              
%     ax(1) = subplot(211); hist(cc, aa);
%     ax(2) = subplot(212); hist(ss, aa);
%     linkaxes(ax, 'x');
% xlim([0 250])

% Find out distance of dF's lower than 10 Hz
dists = []; dFs = [];
for kk=18:22
   for jj = 1:length(outs(kk).pair)
    dists = [dists outs(kk).pair(jj).descartes(outs(kk).pair(jj).dF < 10)];
    dFs = [dFs outs(kk).pair(jj).dF(outs(kk).pair(jj).dF < 10)];
   end    
end

%% Covariance between distance and dF for outc
posneg = []; fidx = [];
sigsig = []; pairdx = [];
dist = []; mdF = [];

for fn = 1:17
    if length(outc(fn).pair) > 1
        for xx = 1:length(outc(fn).pair)
            %figure(xx); clf; 
            
        if length(outc(fn).pair(xx).dF) > 50 % make sure that we have data in there
            if outc(fn).pair(xx).covDistdFpval(1,2) < 0.0000001; sigsig(end+1) = 1; end
            if outc(fn).pair(xx).covDistdFpval(1,2) > 0.0000001; sigsig(end+1) = 0; end
            if outc(fn).pair(xx).covDistdF(1,2) > 0; posneg(end+1) = 1; end
            if outc(fn).pair(xx).covDistdF(1,2) < 0; posneg(end+1) = -1; end
            dist(end+1) = outc(fn).pair(xx).meanDist;   
            mdF(end+1) = outc(fn).pair(xx).dFmean;   
            fidx(end+1) = fn;
            pairdx(end+1) = xx;
        end
        end
    end
end
sigposidx = find(posneg == 1 & sigsig ==1);
signegidx = find(posneg == -1 & sigsig == 1);
nsigidx = find(sigsig == 0);
%% Covariance between distance and dF for outs
posnegS = []; fidxS = [];
sigsigS = []; pairdxS = [];
distS = []; mdFS = [];

for fn = 18:23
    if length(outs(fn).pair) > 1
        for xx = 1:length(outs(fn).pair)
            %figure(xx); clf; 
            
        if length(outs(fn).pair(xx).dF) > 50 % make sure that we have data in there
            if outs(fn).pair(xx).covDistdFpval(1,2) < 0.0000001; sigsigS(end+1) = 1; end
            if outs(fn).pair(xx).covDistdFpval(1,2) > 0.0000001; sigsigS(end+1) = 0; end
            if outs(fn).pair(xx).covDistdF(1,2) > 0; posnegS(end+1) = 1; end
            if outs(fn).pair(xx).covDistdF(1,2) < 0; posnegS(end+1) = -1; end
            distS(end+1) = outs(fn).pair(xx).meanDist;   
            mdFS(end+1) = outs(fn).pair(xx).dFmean;   
            fidxS(end+1) = fn;
            pairdxS(end+1) = xx;
        end
        end
    end
end
sigposidxS = find(posnegS == 1 & sigsigS ==1);
signegidxS = find(posnegS == -1 & sigsigS == 1);
nsigidxS = find(sigsigS == 0);

%%
    plot(5*(outc(fn).pair(xx).dF - mean(outc(fn).pair(xx).dF)), '.')   
        hold on; 
    plot(cvcv*(outc(fn).pair(xx).descartes - mean(outc(fn).pair(xx).descartes)), '.'); 
    text(100,10,num2str(outc(fn).pair(xx).covDistdF(1,2)));
    text(100,5,num2str(outc(fn).pair(xx).covDistdFpval(1,2)));
    
% figure;
% yyy = xcorr(outc(fn).pair(xx).descartes - mean(outc(fn).pair(xx).descartes), -(outc(fn).pair(xx).dF - mean(outc(fn).pair(xx).dF)));
% plot(yyy);

%% Plot examples for poster

figure(27); clf;
subplot(211); plot((outc(3).pair(1).dF-mean(outc(3).pair(1).dF))*10, '.');
    hold on; plot(-(outc(3).pair(1).descartes-mean(outc(3).pair(1).descartes)), '.');         
figure(27); subplot(212); plot((outc(8).pair(8).dF-mean(outc(8).pair(8).dF))*10, '.');    
    hold on; plot((outc(8).pair(8).descartes-mean(outc(8).pair(8).descartes)), '.');     

    
%% Plot all trajectories


figure(27); clf; hold on;
for fn = 18:23
        for xx = 1:length(f(fn).s)
            plot(f(fn).s(xx).tx, f(fn).s(xx).ty, '.');
        end
end

figure(28); clf; hold on;
for fn = 1:17
        for xx = 1:length(f(fn).c)
            plot(f(fn).c(xx).tx, f(fn).c(xx).ty, '.');
        end
end
            
            

    
    
    
    
    
    
    