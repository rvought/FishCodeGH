addpath('TRENTOOL3/');
addpath(genpath('fieldtrip-20190202/'))

load CaveDataRev2018a.mat

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

%%

data.time = {t};
data.trial = {[dist(:,1),df(:,1)]'};
data.label = {'dist','df'};
data.fsample = 1/mean(diff(t));

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
cfgTEP.maxlag      = 1000;
cfgTEP.minnrtrials = 15;   % minimum acceptable number of trials

% optimizing embedding
cfgTEP.optimizemethod ='ragwitz';  % criterion used
cfgTEP.ragdim         = 2:9;       % criterion dimension
cfgTEP.ragtaurange    = [0.2 0.4]; % range for tau
cfgTEP.ragtausteps    = 5;        % steps for ragwitz tau steps
cfgTEP.repPred        = 100;      % size(data.trial{1,1},2)*(3/4);

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
cfgTESS.fileidout  = strcat(OutputDataPath,'Lorenzdata_1->2_');

%% calculation - scan over specified values for u

TGA_results = InteractionDelayReconstruction_calculate(cfgTEP,cfgTESS,data);