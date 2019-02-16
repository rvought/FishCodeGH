% addpath('TRENTOOL3/');
% addpath(genpath('fieldtrip-20190202/'))
% rmpath(fullfile('fieldtrip-20190202','external','signal'));
% rmpath(genpath(fullfile('fieldtrip-20190202','compat')));

addpath('tewp');

load CaveDataRev2018a.mat

mmnorm = @(x) ( x - repmat(min(x),size(x,1),1) ) ./ ( repmat(max(x),size(x,1),1) - repmat(min(x),size(x,1),1) );


%% Compute dfs and distances

c = 5;

nFish = length(cave(c).fish);
t = cave(c).t;
nTimes = length(t);
dt = mean(diff(t));

% Pairwise states
C = nchoosek(1:nFish,2);
nPairs = size(C,1);
[dist,df] = deal(zeros(nTimes,nPairs));
for k = 1:nPairs
    dist(:,k) = sqrt((cave(c).fish(C(k,1)).x - cave(c).fish(C(k,2)).x).^2 + (cave(c).fish(C(k,1)).y - cave(c).fish(C(k,2)).y).^2);
    df(:,k) = abs(cave(c).fish(C(k,1)).freq(:,2) - cave(c).fish(C(k,2)).freq(:,2));
end
dist(isnan(df)) = NaN;

dist = fillmissing(dist,'spline');
df = fillmissing(df,'spline');

smoothInterval = 15; % seconds
smoothLength = round(smoothInterval/dt);

dist = smoothdata(dist,'movmean',smoothLength,'omitnan');
df = smoothdata(df,'movmean',smoothLength,'omitnan');

% dist = highpass(smoothdata(dist,'movmean',smoothLength,'omitnan'),0.001,1/dt,'ImpulseResponse','iir');
% df = highpass(smoothdata(df,'movmean',smoothLength,'omitnan'),0.001,1/dt,'ImpulseResponse','iir');
%%
windowInterval = 60; % seconds
windowLength = round(windowInterval/dt);
percOverlap = 95;
overlapLength = floor(windowLength*percOverlap/100);

[DST,DF] = deal([]);
for k = 1:nPairs
    DST(:,:,k) = buffer(dist(:,k),windowLength,overlapLength,'nodelay');
    DF(:,:,k) = buffer(df(:,k),windowLength,overlapLength,'nodelay');
end
T = buffer(t,windowLength,overlapLength,'nodelay');

nWindows = size(DST,2)-1;
maxDelay = round(windowLength*0.9);

%%


[TE_DST_DF,TE_DF_DST] = deal(zeros(nWindows,maxDelay,nPairs));

parfor c = 1:nPairs
    for k = 1:nWindows
        for j = 1:maxDelay
    %         TE(k,j) = transferEntropyKDE(squeeze(DST(:,k,c)),squeeze(DF(:,k,c)),j,0,20,1);
            TE_DST_DF(k,j,c) = transferEntropyPartition(squeeze(DST(:,k,c)),squeeze(DF(:,k,c)),0,j);
            TE_DF_DST(k,j,c) = transferEntropyPartition(squeeze(DF(:,k,c)),squeeze(DST(:,k,c)),0,j);
        end
    end
end

save te_results2 TE_DST_DF TE_DF_DST

%%

% smoothInterval = 30; % seconds
% smoothLength = round(smoothInterval/dt);
% mean_dist = highpass(smoothdata(dist,'movmean',smoothLength,'omitnan'),0.001,1/dt,'ImpulseResponse','iir');
% mean_df = highpass(smoothdata(df,'movmean',smoothLength,'omitnan'),0.001,1/dt,'ImpulseResponse','iir');

[M,D] = deal(zeros(nWindows,nPairs));
for c = 1:nPairs
%     [M(:,c),idx] = max(smoothdata(squeeze(TE_DST_DF(:,:,c))','movmean',smoothLength/4));
%     D(:,c) = idx*dt;
    [M(:,c),idx] = max(smoothdata(squeeze(TE_DF_DST(:,:,c))','movmean',smoothLength/4));
    D(:,c) = idx*dt;
end

for c = 1:nPairs
    clf, hold on;
    title(sprintf('Pair %d: %d<->%d',c,C(c,1),C(c,2)));

    plot(t,mmnorm(dist(:,c)));
    plot(t,mmnorm(df(:,c)));

    plot(mean(T(:,1:end-1)),mmnorm(D(:,c)),'.');
    plot(mean(T(:,1:end-1)),mmnorm(M(:,c)));
    
%     plot(mean(T(:,1:end-1)),mmnorm(D2(:,c)),'.');
%     plot(mean(T(:,1:end-1)),mmnorm(M2(:,c)));

    legend('Dist','DF','delay1','Max_TE1','delay2','Max_TE2');
    hold off;
    pause;
end

%%

% crange = [-1,1];
crange = [min(M(:)),max(M(:))];
% crange = [min(D(:)),max(D(:))*0.75];

for j = 1:nWindows
    clf, hold on;
    for k = 1:nPairs
        rectangle('Position',[C(k,1),C(k,2),1,1],'FaceColor',vals2colormap(M(j,k),'jet',crange))
    end
    hold off;
    pause;
end

%%

% crange = [-1,1];
% crange = [min(M(:)),max(M(:))];
crange = [min(D(:)),10];

clf, hold on;
for j = 1:100%nActualWindows-1
    for k = 1:nPairs
        plotcube([1,1,1],[C(k,1),C(k,2),j],0.2,vals2colormap(D(j,k),'jet',crange))
    end
end
hold off;


%% define cfg for TEprepare.m

c = 1;
data.trial = {};
data.time = {};
for k = 1:size(DST,2)-1
    data.time = [data.time,T(:,k)'];
    data.trial = [data.trial,[squeeze(DST(:,k,c))';squeeze(DF(:,k,c))']];
end

data.label = {'dist','df'};
data.fsample = 1/dt;


cfgTEP = [];


% data
cfgTEP.toi                 = [min(data.time{1,1}),max(data.time{1,1})]; % time of interest
cfgTEP.sgncmb              = {'dist' 'df'};  % channels to be analyzed

% scanning of interaction delays u
cfgTEP.predicttimemin_u    = 0;      % minimum u to be scanned
cfgTEP.predicttimemax_u    = 10000;	  % maximum u to be scanned
cfgTEP.predicttimestepsize = 100; 	  % time steps between u's to be scanned

% estimator
cfgTEP.TEcalctype  = 'VW_ds'; % use the new TE estimator (Wibral, 2013)

% ACT estimation and constraints on allowed ACT(autocorelation time)
cfgTEP.actthrvalue = 100;   % threshold for ACT
cfgTEP.maxlag      = 100;
cfgTEP.minnrtrials = 1;   % minimum acceptable number of trials

% optimizing embedding
cfgTEP.optimizemethod ='ragwitz';  % criterion used
cfgTEP.ragdim         = 2:6;       % criterion dimension
cfgTEP.ragtaurange    = [0.2 0.4]; % range for tau
cfgTEP.ragtausteps    = 5;        % steps for ragwitz tau steps
cfgTEP.repPred        = 100;         % size(data.trial{1,1},2)*(3/4);

cfgTEP.trialselect    = 'no';   % all trials
% cfgTEP.ensemblemethod = 'yes';

% kernel-based TE estimation
cfgTEP.flagNei = 'Mass' ;           % neigbour analyse type
cfgTEP.sizeNei = 4;                 % neigbours to analyse

%% define cfg for TEsurrogatestats_ensemble.m

cfgTESS = [];

% use individual dimensions for embedding
cfgTESS.optdimusage = 'indivdim';

% statistical and shift testing
cfgTESS.tail           = 1;
cfgTESS.numpermutation = 5e4;
cfgTESS.shifttesttype  ='TEshift>TE';
cfgTESS.surrogatetype  = 'trialshuffling';

% results file name
OutputDataPath = './trentool_results/';
cfgTESS.fileidout  = strcat(OutputDataPath,'Lorenzdata_1->2_');

%% calculation - scan over specified values for u

TGA_results = InteractionDelayReconstruction_calculate(cfgTEP,cfgTESS,data);

%%

mmnorm = @(x) ( x - repmat(min(x),size(x,1),1) ) ./ ( repmat(max(x),size(x,1),1) - repmat(min(x),size(x,1),1) );

smoothInterval = 30; % seconds
smoothLength = round(smoothInterval/dt);
mean_dist = highpass(smoothdata(dist,'movmean',smoothLength,'omitnan'),0.001,1/dt,'ImpulseResponse','iir');
mean_df = highpass(smoothdata(df,'movmean',smoothLength,'omitnan'),0.001,1/dt,'ImpulseResponse','iir');

clf, hold on;

plot(t,mmnorm(mean_dist(:,c)));
plot(t,mmnorm(mean_df(:,c)));
plot(mean(T(:,1:end-1)),TGA_results.TEmat);
plot(mean(T(:,1:end-1)),TGA_results.MImat);

legend('Dist','DF','TE','MI');