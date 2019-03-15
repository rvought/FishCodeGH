addpath('TRENTOOL3/');
addpath(genpath('fieldtrip-20190202/'))
rmpath(fullfile('fieldtrip-20190202','external','signal'));
rmpath(genpath(fullfile('fieldtrip-20190202','compat')));
ft_defaults;

load CaveDataRev2018a.mat
load SurfaceDataRev2018a.mat

nCave = length(cave);
nSrf = length(srf);

mmnorm = @(x) ( x - repmat(min(x),size(x,1),1) ) ./ ( repmat(max(x),size(x,1),1) - repmat(min(x),size(x,1),1) );

%%

data.trial = {};
data.time = {};
[meta.dataset,meta.pair,meta.window,meta.df,meta.ddf,meta.dist,meta.ddist] = deal([]);

for j = 1:(nCave + nSrf)
    if j <= nCave
        dat = cave(j);
        dataType = 'cave';
    else
        dat = srf(j-nCave);
        dataType = 'srf';
    end
    
    nFish = length(dat.fish);
    if nFish>1
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

        smoothInterval = 30; % seconds
        smoothLength = round(smoothInterval/dt);

        dist = smoothdata(dist,'movmean',smoothLength,'omitnan');
        df = smoothdata(df,'movmean',smoothLength,'omitnan');

        ddist = diff(dist);
        ddf = diff(df);
        
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
            
            for w = 1:(size(DDST,2)-1)
%                 data.time = [data.time,T(:,w)'];
                data.time = [data.time,T];
                data.trial = [data.trial,[DDST(:,w),DDF(:,w)]'];
                
                meta.dataset = [meta.dataset,j];
                meta.pair = [meta.pair,k];
                meta.window = [meta.window,w];
                meta.dist = [meta.dist,nanmean(DST(:,w))];
                meta.ddist = [meta.ddist,nanmean(DDST(:,w))];
                meta.df = [meta.df,nanmean(DF(:,w))];
                meta.ddf = [meta.ddf,nanmean(DDF(:,w))];
                
%                 [r,lag] = xcorr(mmnorm(DST(:,w)),mmnorm(DF(:,30)),'coeff');
%                 [~,i] = max(abs(r));
            end
        end
    end
end

data.label = {'ddist','ddf'};
data.fsample = 1/mean(diff(dat.t));


%%  TRENTOOL - CPU method


OutputDataPath = './results/';


% define cfg for TEprepare.m
cfgTEP = [];


% data
cfgTEP.toi                 = [min(data.time{1,1}),max(data.time{1,1})]; % time of interest
cfgTEP.sgncmb              = {'ddist' 'ddf'};  % channels to be analyzed

% scanning of interaction delays u
cfgTEP.predicttimemin_u    = 0;      % minimum u to be scanned (ms)
cfgTEP.predicttimemax_u    = 30000;	  % maximum u to be scanned  (ms)
cfgTEP.predicttimestepsize = 1000; 	  % time steps between u's to be scanned

% estimator
cfgTEP.TEcalctype  = 'VW_ds'; % use the new TE estimator (Wibral, 2013)

% ACT estimation and constraints on allowed ACT(autocorelation time)
cfgTEP.actthrvalue = 100;   % threshold for ACT
cfgTEP.maxlag      = round(100/dt); % samples
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

% use ensemble method
cfgTEP.ensemblemethod = 'no';

% define cfg for TEsurrogatestats_ensemble.m
cfgTESS = [];

% use individual dimensions for embedding
cfgTESS.optdimusage = 'indivdim';

% statistical and shift testing
cfgTESS.tail           = 1;
cfgTESS.numpermutation = 5e4;
cfgTESS.shifttesttype  ='TEshift>TE';
cfgTESS.surrogatetype  = 'trialshuffling';

% results file name
cfgTESS.fileidout  = strcat(OutputDataPath,'all_windowed');

% calculation - scan over specified values for u

TGA_results = InteractionDelayReconstruction_calculate(cfgTEP,cfgTESS,data);


%%

load TGA_results;

%%
act = squeeze(TGA_results.ACT.actvalue(:,2,:));
teIdx = act>=0 & act<=100;

TE = NaN(size(teIdx));
TE(teIdx) = TGA_results.TEmat;

d = 3;
p = 10;

idx = meta.dataset==d & meta.pair==p;
time = data.time(idx);
trial = data.trial(idx);

% t = cellfun(@(x) mean(x),time);
ddist = cellfun(@(x) mean(x(1,:)),trial);
ddf = cellfun(@(x) mean(x(2,:)),trial);
t = (1:length(ddist))*(windowLength-overlapLength)*dt;

clf, hold on
plot(t,mmnorm(meta.dist(idx)));
plot(t,mmnorm(meta.df(idx)));
plot(t,TE(idx));
hold off;

%%

clf, hold on
for d = 1:(nCave + nSrf)
    if d <= nCave
        dat = cave(d);
    else
        dat = srf(d-nCave);
    end
    
    nFish = dat.nFish;
    freq = nanmean([dat.fish.freq]);
    freq = freq(2:2:end);
    
    vel = nanmean(sqrt(diff([dat.fish.x]).^2 + diff([dat.fish.y]).^2 + diff([dat.fish.z]).^2));
    if nFish>1
        C = nchoosek(1:nFish,2);
        nPairs = size(C,1);
        
        te_pair = zeros(nPairs,1);
        for p = 1:nPairs
            idx = meta.dataset==d & meta.pair==p;
            te_pair(p) = nanmean(TE(idx));
        end

        te_fish = zeros(nFish,1);
        te_fish_var = zeros(nFish,1);
        for f = 1:nFish
            fishIdx = find(sum(C==f,2));
            te_fish(f) = nanmean(te_pair(fishIdx));
            te_fish_var(f) = nanstd(te_pair(fishIdx));
            
%             plot(freq(f),te_pair(fishIdx),'.');
        end
        
        if d<=nCave
            col = 'r';
        else
            col = 'b';
        end
        
        subplot(3,1,1), hold on;
            plot(freq,te_fish,'.','Color',col);
        subplot(3,1,2), hold on;
            plot(freq,te_fish_var,'.','Color',col);
        subplot(3,1,3), hold on;
            plot(freq,te_fish./te_fish_var,'.','Color',col);

%         plot(freq,vel,'.','Color',col);
%         plot(vel,te_fish,'.','Color',col);
%         plot(vel,te_fish_var,'.','Color',col);
%         plot(vel,te_fish./te_fish_var,'.','Color',col);
    end
end

% subplot(3,1,1), grid on;
% subplot(3,1,2), grid on;
% subplot(3,1,3), grid on;

hold off;