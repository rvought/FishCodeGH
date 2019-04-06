%% Load Tracked Data
dataFolderName = '~/Downloads';
% dataFolderName = '/Users/Shared/Data/Brasil';
load(fullfile(dataFolderName,'CaveDataRev2018a.mat'));
load(fullfile(dataFolderName,'SurfaceDataRev2018a.mat'));

nCave = length(cave);
nSrf = length(srf);

mmnorm = @(x) ( x - repmat(min(x),size(x,1),1) ) ./ ( repmat(max(x),size(x,1),1) - repmat(min(x),size(x,1),1) );

%% Compile data and create windowed data for TRENTOOL

data.trial = {};
data.time = {};
[meta.dataset,meta.pair,meta.fish1,meta.fish2,...
    meta.window,meta.df,meta.ddf,meta.ddf_abs,...
    meta.dist,meta.ddist,meta.ddist_abs] = deal([]);

[freq_all,dist_all,df_all,ddist_all,ddf_all,ddist_abs_all,ddf_abs_all] = deal(cell(nCave+nSrf,1));

nFish_all = zeros(nCave+nSrf,1);

for j = 1:(nCave + nSrf)
    % Use both cave and surface data
    if j <= nCave
        dat = cave(j);
        dataType = 'cave';
    else
        dat = srf(j-nCave);
        dataType = 'srf';
    end
    
    nFish = length(dat.fish);
    nFish_all(j) = nFish;
    
    for k = 1:nFish
        freq_all{j} = [dat.fish.freq];
        freq_all{j} = freq_all{j}(:,2:2:end);
    end

    if nFish>1  % analyze if there is at least one pair
        t = dat.t;
        nTimes = length(t);
        dt = mean(diff(t));

        % Pairwise states
        C = nchoosek(1:nFish,2);
        nPairs = size(C,1);
        [dist,df] = deal(zeros(nTimes,nPairs));
        for k = 1:nPairs
            dist(:,k) = sqrt((dat.fish(C(k,1)).x - dat.fish(C(k,2)).x).^2 + (dat.fish(C(k,1)).y - dat.fish(C(k,2)).y).^2);
            df(:,k) = abs(dat.fish(C(k,1)).freq(:,2) - dat.fish(C(k,2)).freq(:,2));
        end
       
        nanIdx = isnan(df);
        dist(nanIdx) = NaN;

        % Smooth distance and df
        smoothInterval = 30; % seconds
        smoothLength = round(smoothInterval/dt);
        dist = smoothdata(dist,'movmean',smoothLength,'omitnan');
        df = smoothdata(df,'movmean',smoothLength,'omitnan');

        dist_all{j} = dist;
        df_all{j} = df;
        
        % Compute and smooth derivatives
        ddist = diff(dist);
        ddf = diff(df);
      
        ddist_abs = abs(ddist);
        ddf_abs = abs(ddf);
        
        ddist = smoothdata(ddist,'movmean',smoothLength,'omitnan');
        ddf = smoothdata(ddf,'movmean',smoothLength,'omitnan');
        
        ddist_all{j} = ddist;
        ddf_all{j} = ddf;
        
        ddist_abs = smoothdata(ddist_abs,'movmean',smoothLength,'omitnan');
        ddf_abs = smoothdata(ddf_abs,'movmean',smoothLength,'omitnan');
        
        ddist_abs_all{j} = ddist_abs;
        ddf_abs_all{j} = ddf_abs;
        
        % Window and create data pairs
        windowInterval = 200; % seconds
        windowLength = round(windowInterval/dt);
        percOverlap = 90;
        overlapLength = floor(windowLength*percOverlap/100);
        
%         T = buffer(t(1:end-1),windowLength,overlapLength,'nodelay');
        T = (1:windowLength)*dt;
        for k = 1:nPairs
            DST = buffer(dist(:,k),windowLength,overlapLength,'nodelay');
            DF = buffer(df(:,k),windowLength,overlapLength,'nodelay');
            
            DDST = buffer(ddist(:,k),windowLength,overlapLength,'nodelay');
            DDF = buffer(ddf(:,k),windowLength,overlapLength,'nodelay');
            
            DDST_ABS = buffer(ddist_abs(:,k),windowLength,overlapLength,'nodelay');
            DDF_ABS = buffer(ddf_abs(:,k),windowLength,overlapLength,'nodelay');
            
            for w = 1:(size(DDST,2)-1)
%                 data.time = [data.time,T(:,w)'];
                data.time = [data.time,T];
                data.trial = [data.trial,[DST(:,w),DF(:,w),DDST(:,w),DDF(:,w),DDST_ABS(:,w),DDF_ABS(:,w)]'];
                
                meta.dataset = [meta.dataset,j];
                meta.pair = [meta.pair,k];
                meta.fish1 = [meta.fish1,C(k,1)];
                meta.fish2 = [meta.fish2,C(k,2)];
                meta.window = [meta.window,w];
                
                meta.dist = [meta.dist,nanmean(DST(:,w))];
                meta.ddist = [meta.ddist,nanmean(DDST(:,w))];
                meta.ddist_abs = [meta.ddist_abs,nanmean(DDST_ABS(:,w))];
                
                meta.df = [meta.df,nanmean(DF(:,w))];
                meta.ddf = [meta.ddf,nanmean(DDF(:,w))];
                meta.ddf_abs = [meta.ddf_abs,nanmean(DDF_ABS(:,w))];
                
                if nanmean(DDF_ABS(:,w))<0
                    error;
                end
            end
        end
    end
end

%% Compute sensitivity and correlation values

[sens,sens_abs,maxr,maxr_abs,lag,lag_abs] = deal(zeros(length(data.trial),1));
for k = 1:length(data.trial)
    sens(k) = nanmean(data.trial{k}(3,:)./data.trial{k}(4,:));
    sens_abs(k) = nanmean(data.trial{k}(5,:)./data.trial{k}(6,:));
    
    [r,l] = xcorr(data.trial{k}(3,:),data.trial{k}(4,:));
    [maxr(k),idx] = max(abs(r));
    lag(k) = l(idx);
    
    [r,l] = xcorr(data.trial{k}(5,:),data.trial{k}(6,:));
    [maxr_abs(k),idx] = max(abs(r));
    lag_abs(k) = l(idx);
end

%% Example of low-df behavior, Cave fish
% There is only one pair where the df is kinda small (<5) and that is a crossing pair 
idx = meta.df<5 & meta.dataset<=nCave;
dataset = unique(meta.dataset(idx));
pair = unique(meta.pair(idx));
time = [0,60];

dat = cave(dataset);
C = nchoosek(1:dat.nFish,2);

timeIdx = dat.t>time(1) & dat.t<time(2);

figure(1); clf, hold on;

plot(dat.t(timeIdx), dat.fish(C(pair,1)).freq(timeIdx,2), 'LineWidth', 2);
plot(dat.t(timeIdx), dat.fish(C(pair,2)).freq(timeIdx,2), 'LineWidth', 2);

title('Single low-df example in cave fish, crossing pair');
xlabel('Time (s)');
ylabel('Frequency (Hz)');

ylim([250, 300]);
hold off;

%% Example of low-df behavior, Surface fish

idx = meta.df<3 & meta.dataset>nCave;

% Choose the pair with the largest number of windows
longestIdx = mode(meta.dataset(idx) + 1i*meta.pair(idx));
dataset = real(longestIdx);
pair = imag(longestIdx);
time = [0,dat.t(end)];

dat = srf(dataset-nCave);
C = nchoosek(1:dat.nFish,2);

timeIdx = dat.t>time(1) & dat.t<time(2);

figure(2); clf, hold on;

plot(dat.t(timeIdx),dat.fish(C(pair,1)).freq(timeIdx,2), 'LineWidth', 2);
plot(dat.t(timeIdx),dat.fish(C(pair,2)).freq(timeIdx,2), 'LineWidth', 2);

title('Longest low-df example in surface fish');
xlabel('Time (s)');
ylabel('Frequency (Hz)');
ylim([320, 370]);
hold off;

%% Compute mean values for each pair

dist_pair_cave = cellfun(@(x) nanmean(x),dist_all(1:nCave),'UniformOutput',false);
dist_pair_cave = [dist_pair_cave{:}];

dist_pair_srf = cellfun(@(x) nanmean(x),dist_all(nCave+1:end),'UniformOutput',false);
dist_pair_srf = [dist_pair_srf{:}];

df_pair_cave = cellfun(@(x) nanmean(x),df_all(1:nCave),'UniformOutput',false);
df_pair_cave = [df_pair_cave{:}];

df_pair_srf = cellfun(@(x) nanmean(x),df_all(nCave+1:end),'UniformOutput',false);
df_pair_srf = [df_pair_srf{:}];

ddist_pair_cave = cellfun(@(x) nanmean(x),ddist_all(1:nCave),'UniformOutput',false);
ddist_pair_cave = [ddist_pair_cave{:}];

ddist_pair_srf = cellfun(@(x) nanmean(x),ddist_all(nCave+1:end),'UniformOutput',false);
ddist_pair_srf = [ddist_pair_srf{:}];

ddf_pair_cave = cellfun(@(x) nanmean(x),ddf_all(1:nCave),'UniformOutput',false);
ddf_pair_cave = [ddf_pair_cave{:}];

ddf_pair_srf = cellfun(@(x) nanmean(x),ddf_all(nCave+1:end),'UniformOutput',false);
ddf_pair_srf = [ddf_pair_srf{:}];

ddist_abs_pair_cave = cellfun(@(x) nanmean(x),ddist_abs_all(1:nCave),'UniformOutput',false);
ddist_abs_pair_cave = [ddist_abs_pair_cave{:}];

ddist_abs_pair_srf = cellfun(@(x) nanmean(x),ddist_abs_all(nCave+1:end),'UniformOutput',false);
ddist_abs_pair_srf = [ddist_abs_pair_srf{:}];

ddf_abs_pair_cave = cellfun(@(x) nanmean(x),ddf_abs_all(1:nCave),'UniformOutput',false);
ddf_abs_pair_cave = [ddf_abs_pair_cave{:}];

ddf_abs_pair_srf = cellfun(@(x) nanmean(x),ddf_abs_all(nCave+1:end),'UniformOutput',false);
ddf_abs_pair_srf = [ddf_abs_pair_srf{:}];

%% For each pair - Distance

max_val = max([dist_pair_cave,dist_pair_srf]);
edges = linspace(0,max_val,20);

figure(3); clf;
hold on;

ch = histogram(dist_pair_cave,edges);
    chcounts = ch.BinCounts; chx = ch.BinEdges(2:end);
sh = histogram(dist_pair_srf,edges);
    shcounts = sh.BinCounts; shx = sh.BinEdges(2:end);

    legend('Cave','Surface');
    xlabel('Distance');
    ylabel('Pairs')
    title('Distance - cave vs surface');
    hold off;

figure(4); clf; hold on;
plot(chx, chcounts/max(chcounts)); plot(shx, shcounts/max(shcounts));    

    legend('Cave','Surface');
    xlabel('Distance');
    ylabel('Pairs')
    title('Distance - cave vs surface');
    hold off;

p = compareDistributions(dist_pair_cave,dist_pair_srf);
fprintf('\nProbability that the samples are from the same distribution: %.2f\n',p);

%% For each pair - df

max_val = max([df_pair_cave,df_pair_srf]);
edges = linspace(0,max_val,40);

figure(5); clf; subplot(211); hold on;

yyaxis left;
dFc = histogram(df_pair_cave,edges);
    dFcounts = dFc.BinCounts; dFcx = dFc.BinEdges(2:end);
yyaxis right;
dFs = histogram(df_pair_srf,edges);
    dFscounts = dFs.BinCounts; dFsx = dFs.BinEdges(2:end);

    legend('Cave','Surface');
    xlabel('Df (Hz)');
    ylabel('Pairs')
    title('Df - cave vs surface');
    hold off;

figure(5); subplot(212); hold on;
plot(dFcx, dFcounts/max(dFcounts)); plot(dFsx, dFscounts/max(dFscounts));    

    legend('Cave','Surface');
    xlabel('Df (Hz)');
    ylabel('Pairs')
    title('Df - cave vs surface');
    hold off;


p = compareDistributions(df_pair_cave,df_pair_srf);
fprintf('\nProbability that the samples are from the same distribution: %.2f\n',p);

%% For each pair - change in distance

max_val = max([ddist_pair_cave,ddist_pair_srf]);
min_val = min([ddist_pair_cave,ddist_pair_srf]);
edges = linspace(min_val,max_val,100);

figure(6); clf;
hold on;

histogram(ddist_pair_cave,edges);
histogram(ddist_pair_srf,edges);

legend('Cave','Surface');
xlabel('Change in Distance');
ylabel('Pairs')
title('Change in Distance - cave vs surface');
hold off;

p = compareDistributions(ddist_pair_cave,ddist_pair_srf);
fprintf('\nProbability that the samples are from the same distribution: %.2f\n',p);

%% For each pair - change in df

max_val = max([ddf_pair_cave,ddf_pair_srf]);
min_val = min([ddf_pair_cave,ddf_pair_srf]);
edges = linspace(min_val,max_val,100);

clf;
hold on;

histogram(ddf_pair_cave,edges);
histogram(ddf_pair_srf,edges);

legend('Cave','Surface');
xlabel('Change in Df (Hz)');
ylabel('Pairs')
title('Change in Df - cave vs surface');
hold off;

p = compareDistributions(ddf_pair_cave,ddf_pair_srf);
fprintf('\nProbability that the samples are from the same distribution: %.2f\n',p);

%% For each pair - relative speed

max_val = max([ddist_abs_pair_cave,ddist_abs_pair_srf]);
edges = linspace(0,max_val,100);

clf;
hold on;

histogram(ddist_abs_pair_cave,edges);
histogram(ddist_abs_pair_srf,edges);

legend('Cave','Surface');
xlabel('Relative speed');
ylabel('Pairs')
title('Relative speed - cave vs surface');
hold off;

p = compareDistributions(ddist_abs_pair_cave,ddist_abs_pair_srf);
fprintf('\nProbability that the samples are from the same distribution: %.2f\n',p);

%% For each pair - relative frequency speed

max_val = max([ddf_abs_pair_cave,ddf_abs_pair_srf]);
edges = linspace(0,max_val,100);

clf;
hold on;

histogram(ddf_abs_pair_cave,edges);
histogram(ddf_abs_pair_srf,edges);

legend('Cave','Surface');
xlabel('Relative frequency speed (Hz)');
ylabel('Pairs')
title('Relative frequency speed - cave vs surface');
hold off;

p = compareDistributions(ddf_abs_pair_cave,ddf_abs_pair_srf);
fprintf('\nProbability that the samples are from the same distribution: %.2f\n',p);

%% For each pair - relative speed split by direction

max_val = max([ddist_abs_pair_cave,ddist_abs_pair_srf]);
edges = linspace(0,max_val,100);
fulledges = linspace(-2,2,0.10);

figure(7); clf; plot(ddist_pair_cave*5); hold on; plot(ddist_abs_pair_cave);
%histogram(ddist_pair_cave, fulledges,'FaceColor','b');

figure(6); clf;
subplot(2,1,1)
hold on;

histogram(ddist_abs_pair_cave(ddist_pair_cave<0),edges,'FaceColor','b');
histogram(ddist_abs_pair_cave(ddist_pair_cave>0),edges,'FaceColor','c');

legend('Moving closer','Moving further');
ylabel('Pairs');
title('Cave');
hold off;


subplot(2,1,2);
hold on;
histogram(ddist_abs_pair_srf(ddist_pair_srf<0),edges,'FaceColor','r');
histogram(ddist_abs_pair_srf(ddist_pair_srf>0),edges,'FaceColor','m');

title('Surface');
legend('Moving closer','Moving further');
xlabel('Relative speed');
ylabel('Pairs')
hold off;

p = compareDistributions(ddist_abs_pair_cave(ddist_pair_cave<0),ddist_abs_pair_cave(ddist_pair_cave>0));
fprintf('\nProbability that the samples are from the same distribution: %.2f\n',p);

p = compareDistributions(ddist_abs_pair_srf(ddist_pair_srf<0),ddist_abs_pair_srf(ddist_pair_srf>0));
fprintf('\nProbability that the samples are from the same distribution: %.2f\n',p);

%% For each pair - relative speed split by direction

max_val = max([ddf_abs_pair_cave,ddf_abs_pair_srf]);
edges = linspace(0,max_val,100);

clf;
subplot(2,1,1)
hold on;

histogram(ddf_abs_pair_cave(ddist_pair_cave<0),edges,'FaceColor','b');
histogram(ddf_abs_pair_cave(ddist_pair_cave>0),edges,'FaceColor','c');

legend('Moving closer','Moving further');
ylabel('Pairs');
title('Cave');
hold off;


subplot(2,1,2);
hold on;
histogram(ddf_abs_pair_srf(ddist_pair_srf<0),edges,'FaceColor','r');
histogram(ddf_abs_pair_srf(ddist_pair_srf>0),edges,'FaceColor','m');

legend('Moving closer','Moving further');
xlabel('Relative frquency speed');
ylabel('Pairs')
title('Surface');
hold off;


p = compareDistributions(ddf_abs_pair_cave(ddist_pair_cave<0),ddf_abs_pair_cave(ddist_pair_cave>0));
fprintf('\nProbability that the samples are from the same distribution: %.2f\n',p);

p = compareDistributions(ddf_abs_pair_srf(ddist_pair_srf<0),ddf_abs_pair_srf(ddist_pair_srf>0));
fprintf('\nProbability that the samples are from the same distribution: %.2f\n',p);

%% Compute mean values for each fish

[dist_fish_cave,dist_fish_srf,...
df_fish_cave,df_fish_srf,...
ddist_fish_cave,ddist_fish_srf,...
ddf_fish_cave,ddf_fish_srf,...
ddist_abs_fish_cave,ddist_abs_fish_srf,...
ddf_abs_fish_cave,ddf_abs_fish_srf] = deal([]);

for j = 1:(nCave+nSrf)
    nFish = nFish_all(j);
  
    if nFish>2
        C = nchoosek(1:nFish,2);

        for f = 1:nFish
            fishIdx = find(sum(C==f,2));

            if j<=nCave
                dist_fish_cave = [dist_fish_cave,nanmean(nanmean(dist_all{j}(:,fishIdx)))];
                df_fish_cave = [df_fish_cave,nanmean(nanmean(df_all{j}(:,fishIdx)))];
                ddist_fish_cave = [ddist_fish_cave,nanmean(nanmean(ddist_all{j}(:,fishIdx)))];
                ddf_fish_cave = [ddf_fish_cave,nanmean(nanmean(ddf_all{j}(:,fishIdx)))];
                ddist_abs_fish_cave = [ddist_abs_fish_cave,nanmean(nanmean(ddist_abs_all{j}(:,fishIdx)))];
                ddf_abs_fish_cave = [ddf_abs_fish_cave,nanmean(nanmean(ddf_abs_all{j}(:,fishIdx)))];
            else
                dist_fish_srf = [dist_fish_srf,nanmean(nanmean(dist_all{j}(:,fishIdx)))];
                df_fish_srf = [df_fish_srf,nanmean(nanmean(df_all{j}(:,fishIdx)))];
                ddist_fish_srf = [ddist_fish_srf,nanmean(nanmean(ddist_all{j}(:,fishIdx)))];
                ddf_fish_srf = [ddf_fish_srf,nanmean(nanmean(ddf_all{j}(:,fishIdx)))];
                ddist_abs_fish_srf = [ddist_abs_fish_srf,nanmean(nanmean(ddist_abs_all{j}(:,fishIdx)))];
                ddf_abs_fish_srf = [ddf_abs_fish_srf,nanmean(nanmean(ddf_abs_all{j}(:,fishIdx)))];
            end
        end
    end
end

%% For each fish - distance

max_val = max([dist_fish_cave,dist_fish_srf]);
edges = linspace(0,max_val,30);

figure(7); clf;
hold on;

histogram(dist_fish_cave,edges);
histogram(dist_fish_srf,edges);

legend('Cave','Surface');
xlabel('Distance');
ylabel('Pairs')
title('Distance - cave vs surface');
hold off;

p = compareDistributions(dist_fish_cave,dist_fish_srf);
fprintf('\nProbability that the samples are from the same distribution: %.2f\n',p);

%% For each fish - df

max_val = max([df_fish_cave,df_fish_srf]);
edges = linspace(0,max_val,30);

figure(8); clf;
hold on;

histogram(df_fish_cave,edges);
histogram(df_fish_srf,edges);

legend('Cave','Surface');
xlabel('Df');
ylabel('Pairs')
title('Df - cave vs surface');
hold off;

p = compareDistributions(df_fish_cave,df_fish_srf);
fprintf('\nProbability that the samples are from the same distribution: %.2f\n',p);

%% For each fish - change in distance

max_val = max([ddist_fish_cave,ddist_fish_srf]);
edges = linspace(0,max_val,50);

clf;
hold on;

histogram(ddist_fish_cave,edges);
histogram(ddist_fish_srf,edges);

legend('Cave','Surface');
xlabel('Change in Distance');
ylabel('Fish')
title('Change in Distance - cave vs surface');
hold off;

p = compareDistributions(ddist_fish_cave,ddist_fish_srf);
fprintf('\nProbability that the samples are from the same distribution: %.2f\n',p);

%% For each fish - change in df

max_val = max([ddf_fish_cave,ddf_fish_srf]);
edges = linspace(0,max_val,50);

clf;
hold on;

histogram(ddf_fish_cave,edges);
histogram(ddf_fish_srf,edges);

legend('Cave','Surface');
xlabel('Change in Df');
ylabel('Fish')
title('Change in Df - cave vs surface');
hold off;

p = compareDistributions(ddf_fish_cave,ddf_fish_srf);
fprintf('\nProbability that the samples are from the same distribution: %.2f\n',p);

%% For each fish - relative speed

max_val = max([ddist_abs_fish_cave,ddist_abs_fish_srf]);
edges = linspace(0,max_val,50);

clf;
hold on;

histogram(ddist_abs_fish_cave,edges);
histogram(ddist_abs_fish_srf,edges);

legend('Cave','Surface');
xlabel('Relative speed');
ylabel('Fish')
title('Relative speed - cave vs surface');
hold off;

p = compareDistributions(ddist_abs_fish_cave,ddist_abs_fish_srf);
fprintf('\nProbability that the samples are from the same distribution: %.2f\n',p);

%% For each fish - relative frequency speed

max_val = max([ddf_abs_fish_cave,ddf_abs_fish_srf]);
edges = linspace(0,max_val,50);

clf;
hold on;

histogram(ddf_abs_fish_cave,edges);
histogram(ddf_abs_fish_srf,edges);

legend('Cave','Surface');
xlabel('Relative frequency speed (Hz)');
ylabel('Fish')
title('Relative frequency speed - cave vs surface');
hold off;

p = compareDistributions(ddf_abs_fish_cave,ddf_abs_fish_srf);
fprintf('\nProbability that the samples are from the same distribution: %.2f\n',p);

%% Distance vs df plots

figure(9); clf, hold on;
plot(dist_fish_cave,df_fish_cave,'.b','MarkerSize',10)
plot(dist_fish_srf,df_fish_srf,'.r','MarkerSize',10)

grid on;

xlabel('Distance');
ylabel('Df');
hold off;

%% Distance vs df plots - change

figure(10); clf, hold on;
plot(ddist_fish_cave,ddf_fish_cave,'.b','MarkerSize',10)
plot(ddist_fish_srf,ddf_fish_srf,'.r','MarkerSize',10)

grid on;

xlabel('Change in distance');
ylabel('Change in frequency');
hold off;

%% Distance vs df plot - abs change DeltaDistanceVSDeltadF
clf, hold on;

% plot(ddist_abs_fish_cave,ddf_abs_fish_cave,'.b','MarkerSize',10)

plot(ddist_abs_fish_cave(ddist_fish_cave<0),ddf_abs_fish_cave(ddist_fish_cave<0),'.b','MarkerSize',10)
plot(ddist_abs_fish_cave(ddist_fish_cave>0),ddf_abs_fish_cave(ddist_fish_cave>0),'ob','MarkerSize',3)


% px = linspace(0,max(ddist_abs_fish_cave),10);
% p = polyfit(ddist_abs_fish_cave,ddf_abs_fish_cave,1);
% py = polyval(p,px);
% [R_ddist_df_abs_cave,pVal_ddist_df_abs_cave] = corrcoef(ddist_abs_fish_cave,ddf_abs_fish_cave);
% fprintf('\nr = %.2f, n = %d, df = %d, p = %.3e',R_ddist_df_abs_cave(2,1),length(ddist_abs_fish_cave),length(ddist_abs_fish_cave)-2,pVal_ddist_df_abs_cave(2,1));
% plot(px,py,'-b');


% plot(ddist_abs_fish_srf,ddf_abs_fish_srf,'.r','MarkerSize',10)

plot(ddist_abs_fish_srf(ddist_fish_srf<0),ddf_abs_fish_srf(ddist_fish_srf<0),'.r','MarkerSize',10)
plot(ddist_abs_fish_srf(ddist_fish_srf>0),ddf_abs_fish_srf(ddist_fish_srf>0),'or','MarkerSize',3)

% px = linspace(0,max(ddist_abs_fish_srf),10);
% p = polyfit(ddist_abs_fish_srf,ddf_abs_fish_srf,1);
% py = polyval(p,px);
% [R_ddist_df_abs_srf,pVal_ddist_df_abs_srf] = corrcoef(ddist_abs_fish_srf,ddf_abs_fish_srf);
% fprintf('\nr = %.2f, n = %d, df = %d, p = %.3e',R_ddist_df_abs_srf(2,1),length(ddist_abs_fish_srf),length(ddist_abs_fish_srf)-2,pVal_ddist_df_abs_srf(2,1));
% plot(px,py,'-r');


grid on;
xlabel('Relative speed');
ylabel('Relative frequency speed');

hold off;

%% Shuffling statistics between r values

X = [ddist_abs_fish_cave,ddist_abs_fish_srf];
Y = [ddf_abs_fish_cave,ddf_abs_fish_srf];

plot(X,Y,'.');

l1 = length(ddist_abs_fish_cave);
l2 = length(ddist_abs_fish_srf);
L = l1+l2;

N = 10000;
R_rand = zeros(N,1);
for n = 1:N
    idx = randsample(L,l1);
    
    R = corrcoef(X(idx),Y(idx));
    R_rand(n) = R(2,1);
end

clf;
subplot(2,1,1);
hold on;

histogram(R_rand,100);
% plot(R_ddist_df_abs_cave(2,1),0,'*r'); %% ESF ERROR CHANGE

legend('','Cave R value');
title('Shuffled R values for Cave fish');
xlabel('R values');
ylabel('Counts');

hold off;

N = 10000;
R_rand = zeros(N,1);
for n = 1:N
    idx = randsample(L,l2);
    
    R = corrcoef(X(idx),Y(idx));
    R_rand(n) = R(2,1);
end

subplot(2,1,2);
hold on;

histogram(R_rand,100);
% plot(R_ddist_df_abs_cave(2,1),0,'*r');  %% ESF ERROR CHANGE

legend('','Surface R value');
title('Shuffled R values for Surface fish');
xlabel('R values');
ylabel('Counts');

%% Frequency 

[freq_fish_cave,freq_fish_srf] = deal([]);
for j = 1:(nCave+nSrf)
    nFish = nFish_all(j);
  
    if nFish>2
        if j<=nCave
            freq_fish_cave = [freq_fish_cave,nanmean(freq_all{j})];
        else
            freq_fish_srf = [freq_fish_srf,nanmean(freq_all{j})];
        end        
    end
end

max_val = max([freq_fish_cave,freq_fish_srf]);
edges = linspace(0,max_val,50);

figure(11); clf;
hold on;

histogram(freq_fish_cave,edges);
histogram(freq_fish_srf,edges);

legend('Cave','Surface');
xlabel('Frequency (Hz)');
ylabel('Fish')
title('Frequency - cave vs surface');
hold off;

p = compareDistributions(freq_fish_cave,freq_fish_srf);
fprintf('\nProbability that the samples are from the same distribution: %.2f\n',p);

%% Frequency vs distance

clf, hold on;

% plot(freq_fish_cave,dist_fish_cave,'.b');
% plot(freq_fish_srf,dist_fish_srf,'.r');
plot(dist_fish_cave,freq_fish_cave,'.b');
plot(dist_fish_srf,freq_fish_srf, '.r');

hold off;

%% Frequency vs Df

clf, hold on;

plot(freq_fish_cave,df_fish_cave,'.b');
plot(freq_fish_srf,df_fish_srf,'.r');

hold off;


%% Demonstration that Freq vs Df in not a real phenomenon

N = 50;
foo = randn(N,1)*100 + 400;

C = nchoosek(1:N,2);

comb = abs(foo(C(:,1))-foo(C(:,2)));

bar = zeros(N,1);
for k = 1:N
    idx = find(sum(C==k,2));
    bar(k) = mean(comb(idx));
end

plot(foo,bar,'.');


%% Run TRENTOOL analysis

%%
% resultFolderName = '~/Google Drive/trentool_results';
resultFolderName = '/Users/eric/Downloads/trentool_results';
%resultFolderName = '/Users/Shared/Data/Brasil/trentool_results';

load(fullfile(resultFolderName,'all_channels_windowed.mat'));
trentool_result = [trentool_result{:}];

%%

TEprepare = [trentool_result.TEprepare];
u_in_ms = [TEprepare.u_in_ms];
u_in_ms = u_in_ms(1,:);
n_u = length(u_in_ms);

channelPairs = TEprepare(1).channelcombilabel;

nChannelPairs = size(channelPairs,1);

nWindows =  size(TEprepare(1).ACT,3);

TE = nan(n_u,nChannelPairs,nWindows);

for u = 1:n_u
    for c = 1:nChannelPairs
        TE(u,c,TEprepare(u).trials{c,1}) = trentool_result(u).TEmat(c,1:TEprepare(u).nrtrials(c,1));        
    end
end


channelLabels = cell(nChannelPairs,1);
for k = 1:nChannelPairs
    channelLabels{k} = [channelPairs{k,1},' -> ',channelPairs{k,2}];
end

%%

[TEmax,TEmaxIdx] = max(TE,[],1);

TEmax = squeeze(TEmax)';
TEmaxIdx = squeeze(TEmaxIdx)';

idx = TEmax<=0;
TEmax(idx) = NaN;
TEmaxIdx(idx) = NaN;


%% Examples of high-distance, high TE interactions
idx = find(TEmax(:,3)>prctile(TEmax(:,3),95) & meta.dist'>prctile(meta.dist,95));

for k = 1:length(idx)
    clf;
    
    subplot(2,1,1);
    plot(data.trial{idx(k)}(1,:), '.-');
    title('Distance (cm)');
    
    subplot(2,1,2);
    plot(data.trial{idx(k)}(2,:), '.-');
    title('Df (Hz)');
    pause;
end


%% Examples of high-df, high TE interactions
idx = find(TEmax(:,3)>prctile(TEmax(:,3),95) & meta.df'>prctile(meta.df,95));

for k = 1:length(idx)
    clf;
    
    subplot(2,1,1);
    plot(data.trial{idx(k)}(1,:), '.-');
    title('Distance (cm)');
    
    subplot(2,1,2);
    plot(data.trial{idx(k)}(2,:), '.-');
    title('Df (Hz)');
    pause;
end


%% TE comparison Cave <-> Surface pairs

[TE_pair_cave,TE_pair_srf,TE_pair_cave_std,TE_pair_srf_std] = deal([]);

for j = 1:(nCave+nSrf)
    nFish = nFish_all(j);
  
    if nFish>2
        C = nchoosek(1:nFish,2);
        
        nPairs = size(C,1);
        
        for p = 1:nPairs
            idx = meta.dataset==j & meta.pair==p;
            
            if j<=nCave
                TE_pair_cave = [TE_pair_cave;nanmean(TEmax(idx,:))];
                TE_pair_cave_std = [TE_pair_cave_std;nanstd(TEmax(idx,:))];
            else
                TE_pair_srf = [TE_pair_srf;nanmean(TEmax(idx,:))];
                TE_pair_srf_std = [TE_pair_srf_std;nanstd(TEmax(idx,:))];
            end
        end
    end
end

clf;
for k = 1:nChannelPairs
    subplot(nChannelPairs,1,k);
    hold on    
    
    max_val = max([TE_pair_cave(:,k);TE_pair_srf(:,k)]);
    edges = linspace(0,max_val,50);

    histogram(TE_pair_cave(:,k),edges);
    histogram(TE_pair_srf(:,k),edges);

    legend('Cave','Surface');
    xlabel('TE');
    ylabel('Pairs')
    title(channelLabels{k},'Interpreter','none');
    hold off;
end

%% Proportion of interacting windows for each fish
p = prctile(TEmax(:,3),50);

prop = [];
for j = 1:(nCave+nSrf)
    nFish = nFish_all(j);
  
    if nFish>2
        C = nchoosek(1:nFish,2);
        
        nPairs = size(C,1);
        
        for f = 1:nFish
            fishIdx = find(sum(C==f,2));
            
            idx = meta.dataset==j & ismember(meta.pair,fishIdx);
            
            prop = [prop,sum(TEmax(idx,3)>p)/sum(idx)];

        end
    end
end 

%% TE comparison Cave <-> Surface Fish

[TE_fish_cave,TE_fish_srf,...
    TE_fish_cave_std,TE_fish_srf_std] = deal([]);

for j = 1:(nCave+nSrf)
    nFish = nFish_all(j);
  
    if nFish>2
        C = nchoosek(1:nFish,2);
        
        nPairs = size(C,1);
        
        for f = 1:nFish
            fishIdx = find(sum(C==f,2));
            
            idx = meta.dataset==j & ismember(meta.pair,fishIdx);
            
            if j<=nCave
                TE_fish_cave = [TE_fish_cave;nanmean(TEmax(idx,:))];
                TE_fish_cave_std = [TE_fish_cave_std;nanstd(TEmax(idx,:))];
            else
                TE_fish_srf = [TE_fish_srf;nanmean(TEmax(idx,:))];
                TE_fish_srf_std = [TE_fish_srf_std;nanstd(TEmax(idx,:))];
            end
        end
    end
end

clf;
for k = 1:nChannelPairs
    subplot(nChannelPairs,1,k);
    hold on    
    
    max_val = max([TE_fish_cave(:,k);TE_fish_srf(:,k)]);
    edges = linspace(0,max_val,50);

    histogram(TE_fish_cave(:,k),edges);
    histogram(TE_fish_srf(:,k),edges);

    legend('Cave','Surface');
    xlabel('TE');
    ylabel('Fish')
    title(channelLabels{k},'Interpreter','none');
    hold off;
end

%% TE vs Distance

clf;
for k = 1:nChannelPairs
    subplot(nChannelPairs,1,k);
    hold on 
    
    plot(dist_fish_cave,TE_fish_cave(:,k),'.b');
    plot(dist_fish_srf,TE_fish_srf(:,k),'.r');
    
    xlabel('Distance');
    ylabel('TE')
    title(channelLabels{k},'Interpreter','none');
    
    hold off 
end

%% TE vs DF

clf;
for k = 1:nChannelPairs
    subplot(nChannelPairs,1,k);
    hold on 
    
    plot(df_fish_cave,TE_fish_cave(:,k),'.b');
    plot(df_fish_srf,TE_fish_srf(:,k),'.r');
    
    xlabel('DF');
    ylabel('TE')
    title(channelLabels{k},'Interpreter','none');
    
    hold off 
end

%% TE vs Change in Distance

clf;
for k = 1:nChannelPairs
    subplot(nChannelPairs,1,k);
    hold on 
    
    plot(ddist_fish_cave,TE_fish_cave(:,k),'.b');
    plot(ddist_fish_srf,TE_fish_srf(:,k),'.r');
    
    xlabel('Change in Distance');
    ylabel('TE')
    title(channelLabels{k},'Interpreter','none');
    
    hold off 
end

%% TE vs Change in DF

clf;
for k = 1:nChannelPairs
    subplot(nChannelPairs,1,k);
    hold on 
    
    plot(ddf_fish_cave,TE_fish_cave(:,k),'.b');
    plot(ddf_fish_srf,TE_fish_srf(:,k),'.r');
    
    xlabel('Change in DF');
    ylabel('TE')
    title(channelLabels{k},'Interpreter','none');
    
    hold off 
end

%% TE vs Relative Speed

clf;
for k = 1:nChannelPairs
    subplot(nChannelPairs,1,k);
    hold on 
    
    plot(ddist_abs_fish_cave,TE_fish_cave(:,k),'.b');
    plot(ddist_abs_fish_srf,TE_fish_srf(:,k),'.r');
    
    xlabel('Relative Speed');
    ylabel('TE')
    title(channelLabels{k},'Interpreter','none');
    
    hold off 
end

%% TE vs Relative frequency speed

clf;
for k = 1:nChannelPairs
    subplot(nChannelPairs,1,k);
    hold on 
    
    plot(ddf_abs_fish_cave,TE_fish_cave(:,k),'.b');
    plot(ddf_abs_fish_srf,TE_fish_srf(:,k),'.r');
    
    if k==nChannelPairs
        xlabel('Relative Frequency Speed');
    end
    ylabel('TE')
    title(channelLabels{k},'Interpreter','none');
    
    hold off 
end

%% TE vs Frequency

clf;
for k = 1:nChannelPairs
    subplot(nChannelPairs,1,k);
    hold on 
    
    plot(freq_fish_cave,TE_fish_cave(:,k),'.b');
    plot(freq_fish_srf,TE_fish_srf(:,k),'.r');
    
    if k==nChannelPairs
        xlabel('Relative Frequency Speed');
    end
    ylabel('TE')
    title(channelLabels{k},'Interpreter','none');
    
    hold off 
end

%% All variables vs. TE of ddist_abs->ddf_abs

vars = {'dist_fish','df_fish','ddist_fish','ddf_fish',...
    'ddist_abs_fish','ddf_abs_fish','freq_fish'};

clf;
for k = 1:length(vars)
    subplot(3,3,k);
    hold on;

    x = eval([vars{k},'_cave']);
    x = x(:);
    y = TE_fish_cave(:,3);
    plot(x,y,'.b')

    nanIdx = ~isnan(x) & ~isnan(y);
    x = x(nanIdx);
    y = y(nanIdx);

    p = polyfit(x,y,1);
    R_cave = corrcoef(x,y);
    
    plot(x,polyval(p,x),'-b')
    
    x = eval([vars{k},'_srf']);
    x = x(:);
    y = TE_fish_srf(:,3);
    plot(x,y,'.r')

    nanIdx = ~isnan(x) & ~isnan(y);
    x = x(nanIdx);
    y = y(nanIdx);

    p = polyfit(x,y,1);
    R_srf = corrcoef(x,y);
    
    plot(x,polyval(p,x),'-r')
    title(sprintf('%s, R_c=%.2f, R_s=%.2f',vars{k},R_cave(2,1),R_srf(2,1)),'Interpreter','none');
    hold off;
end


%% Poisson testing of TE (in progress)

TE_cave = TEmax(meta.dataset<=nCave,3);
TE_srf = TEmax(meta.dataset>nCave,3);

N = 10000;
L1 = size(TE_cave,1);
L2 = size(TE_srf,1);
l = 1000;
[m1,s1,m2,s2] = deal(zeros(N,1));

for n = 1:N
    idx1 = randsample(L1,l);
    m1(n) = nanmean(TE_cave(idx1));
    s1(n) = nanstd(TE_cave(idx1));
	
    idx2 = randsample(L2,l);
    m2(n) = nanmean(TE_srf(idx2));
    s2(n) = nanstd(TE_srf(idx2));
end

clf;

edges = linspace(0.6,1.2,30);
subplot(2,1,1);
hold on;

histogram(m1./s1,edges,'Normalization','pdf');
histogram(TE_fish_cave./TE_fish_cave_std,edges,'Normalization','pdf');


hold off;


%% All three variables

clf, hold on;
plot3(ddist_abs_fish_cave,ddf_abs_fish_cave,TE_fish_cave(:,3),'.b','MarkerSize',5)
plot3(ddist_abs_fish_srf,ddf_abs_fish_srf,TE_fish_srf(:,3),'.r','MarkerSize',5)

xlabel('Relative speed');
ylabel('Relative frequency speed');
zlabel('TE');
grid on;
hold off;

%% TE bias - dist->df - df->dist

TE_dir_cave = TE_fish_cave(:,1)-TE_fish_cave(:,2);
TE_dir_srf = TE_fish_srf(:,1)-TE_fish_srf(:,2);

edges = linspace(min([TE_dir_cave;TE_dir_srf]),max([TE_dir_cave;TE_dir_srf]),20);

clf, hold on;
histogram(TE_dir_cave,edges);
histogram(TE_dir_srf,edges);

legend('Cave','Surface');
xlabel('TE bias in the dist -> df direction');
ylabel('Fish');

hold off;

%% TE bias - dist_abs->df_abs - df_abs->dist_abs

TE_dir_abs_cave = TE_fish_cave(:,3)-TE_fish_cave(:,4);
TE_dir_abs_srf = TE_fish_srf(:,3)-TE_fish_srf(:,4);

edges = linspace(min([TE_dir_abs_cave;TE_dir_abs_srf]),max([TE_dir_abs_cave;TE_dir_abs_srf]),20);

clf, hold on;
histogram(TE_dir_abs_cave,edges);
histogram(TE_dir_abs_srf,edges);

legend('Cave','Surface');
xlabel('TE bias in the dist_abs -> df_abs direction');
ylabel('Fish');

hold off;


%% JUNK BELOW THIS

d = 3;
p = 4;

idx = meta.dataset==d & meta.pair==p;


[TEmax,TEidx] = deal(zeros(sum(idx),nChannelPairs));
for k = 1:nChannelPairs
    TEpair = squeeze(TE(:,k,idx));
    
    [TEmax(:,k),TEidx(:,k)] = nanmax(TEpair);
end

negIdx = TEmax<0;
TEmax(negIdx) = NaN;
TEidx(negIdx) = NaN;

clf;

subplot(3,3,1);
hold on
plot(mmnorm(meta.dist(idx)));
plot(mmnorm(meta.df(idx)));
l = legend('dist','df');
hold off

subplot(3,3,4);
hold on
plot(TEmax);
legend(channelLabels,'Interpreter','none');
hold off

subplot(3,3,5);
hold on
plot(TEmax(:,1) - TEmax(:,2));
plot(TEmax(:,3) - TEmax(:,4));
hold off;

subplot(3,3,6);
hold on
plot(TEmax(:,3) - TEmax(:,1));
plot(TEmax(:,4) - TEmax(:,2));
hold off;


subplot(3,3,7);
hold on;
plot(TEidx/n_u);
% legend(channelLabels,'Interpreter','none');
hold off;


subplot(3,3,8);
hold on;
plot(TEidx(:,1) - TEidx(:,2));
plot(TEidx(:,3) - TEidx(:,4));
% legend(channelLabels,'Interpreter','none');
hold off;