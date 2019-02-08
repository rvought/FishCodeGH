% addpath('TRENTOOL3/');
% addpath(genpath('fieldtrip-20190202/'))
% rmpath(fullfile('fieldtrip-20190202','external','signal'));
% rmpath(genpath(fullfile('fieldtrip-20190202','compat')));
% 
% load CaveDataRev2018a.mat

%% Compute dfs and distances

c = 5;

nFish = length(cave(c).fish);
t = cave(c).t;
nTimes = length(t);

% Pairwise states
C = nchoosek(1:nFish,2);
nPairs = size(C,1);
[dist,df] = deal(zeros(nTimes,nPairs));
for k = 1:nPairs
    dist(:,k) = sqrt((cave(c).fish(C(k,1)).x - cave(c).fish(C(k,2)).x).^2 + (cave(c).fish(C(k,1)).y - cave(c).fish(C(k,2)).y).^2);
    df(:,k) = cave(c).fish(C(k,1)).freq(:,2) - cave(c).fish(C(k,2)).freq(:,2);
end
dist(isnan(df)) = NaN;

dist = fillmissing(dist,'spline');
df = fillmissing(df,'spline');


nData = size(dist,1);
nChunks =  10;
nWindow = ceil(nData/nChunks);
percOverlap = 48;
nOverlap = floor(nWindow*percOverlap/100);

[DST,DF] = deal([]);
for k = 1:size(dist,2)
    DST(:,:,k) = buffer(dist(:,k),nWindow,nOverlap,'nodelay');
    DF(:,:,k) = buffer(df(:,k),nWindow,nOverlap,'nodelay');
end
T = buffer(t,nWindow,nOverlap,'nodelay');
    

c = 1;
data.trial = {};
data.time = {};
for k = 1:size(DST,2)-1
    data.time = [data.time,T(:,k)'];
    data.trial = [data.trial,[squeeze(DST(:,k,c))';squeeze(DF(:,k,c))']];
end

data.label = {'dist','df'};

dt = mean(diff(t));
data.fsample = 1/dt;

%% define cfg for TEprepare.m

cfgTEP = [];


% data
cfgTEP.toi                 = [min(data.time{1,1}),max(data.time{1,1})]; % time of interest
cfgTEP.sgncmb              = {'dist' 'df'};  % channels to be analyzed

% scanning of interaction delays u
cfgTEP.predicttimemin_u    = 40;      % minimum u to be scanned
cfgTEP.predicttimemax_u    = 50;	  % maximum u to be scanned
cfgTEP.predicttimestepsize = 1; 	  % time steps between u's to be scanned

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

winInterval = 30; % seconds
winLength = round(winInterval/dt);
mean_dist = highpass(smoothdata(dist,'movmean',winLength,'omitnan'),0.001,1/dt,'ImpulseResponse','iir');
mean_df = highpass(smoothdata(df,'movmean',winLength,'omitnan'),0.001,1/dt,'ImpulseResponse','iir');

clf, hold on;

plot(t,mmnorm(mean_dist(:,c)));
plot(t,mmnorm(mean_df(:,c)));
plot(mean(T(:,1:end-1)),TGA_results.TEmat);