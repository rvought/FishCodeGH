addpath('tewp');

load CaveDataRev2018a.mat
load SurfaceDataRev2018a.mat

nCave = length(cave);
nSrf = length(srf);

%% Compute dfs and distances

data = cell(nCave+nSrf,1);

dt =  mean(diff(cave(1).t));

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

    %     dist = fillmissing(dist,'nearest');
    %     df = fillmissing(df,'nearest');

        data{j}.type = dataType;
        data{j}.nFish = nFish;
        data{j}.nTimes = nTimes;
        data{j}.nPairs = nPairs;
        data{j}.pairs = C;
        data{j}.dist = dist;
        data{j}.df = df;
    end
end

data = [data{:}];

%%

DST = arrayfun(@(x) diff(x.dist),data,'UniformOutput',false);
DST = cellfun(@(x) x(:),DST,'UniformOutput',false);
DST = vertcat(DST{:});
DST = DST(~isnan(DST));
DST = DST(DST>(mean(DST)-3*std(DST)) & DST<(mean(DST)+3*std(DST)));

DF = arrayfun(@(x) diff(x.df),data,'UniformOutput',false);
DF = cellfun(@(x) x(:),DF,'UniformOutput',false);
DF = vertcat(DF{:});
DF = DF(~isnan(DF));
DF = DF(DF>(mean(DF)-3*std(DF)) & DF<(mean(DF)+3*std(DF)));



%%


N = 20;
% Xpati=linspace(min(DST)-0.1*range(DST),max(DST)+0.1*range(DST),N);
% Ypati=linspace(min(DF)-0.1*range(DF),max(DF)+0.1*range(DF),N);
% Yti = Ypati;
% pdf=zeros(N,N,N);


X = diff(data(1).dist(:,1));
Y = diff(data(1).df(:,1));
X=X(:);
Y=Y(:);
N = 20;

Xrange = [min(DST)-0.1*range(DST),max(DST)+0.1*range(DST)];
Yrange = [min(DF)-0.1*range(DF),max(DF)+0.1*range(DF)];

%%

windowInterval = 200; % seconds
windowLength = round(windowInterval/dt);
percOverlap = 90;
overlapLength = floor(windowLength*percOverlap/100);

XB = buffer(X,windowLength,overlapLength,'nodelay');
YB = buffer(Y,windowLength,overlapLength,'nodelay');
T = buffer(cave(1).t,windowLength,overlapLength,'nodelay');

nWindows = size(XB,2)-1;

maxDelay = round(windowLength*0.5);
delayInt = round(1/dt); % Every 1 s;
delays = (1:delayInt:maxDelay);
nDelays = length(delays);

%%
TE = [];
for b = 1:nWindows 
    b
    tic;
    TE{b} = TE_KDE_GPU(XB(:,b),YB(:,b),Xrange,Yrange,delays,30);
    toc;
end
    
%%
TE = [];
for t = 1:10:1000
    t
    T = transferEntropyPartition_mex(X,Y,t,1);
    TE = [TE,T];
end


%%
TE = [];
for t = 1:10:1000
    T = transferEntropyKDE_range(X,Y,Xrange,Yrange,t,1,N,1);
    TE = [TE,T];
end

%%

% fix block lengths at 1
l=1; k=1;

w = 1;
TE = [];
for t = 1:10:100
    t
    
    % go through the time series X and Y, and populate Xpat, Ypat, and Yt
    Xpat=[]; Ypat=[]; Yt=[];
    for i=max([l+t k+w]):1:min([length(X) length(Y)])
        Xpat=[Xpat; X(i-l-t+1:i-t)];
        Ypat=[Ypat; Y(i-k-w+1:i-w)];
        Yt=[Yt; Y(i)];    
    end

    bw_coeff = 1;
    tic;
    for i=1:length(Xpati)
        for j=1:length(Ypati)
            for k=1:length(Yti)
                pdf(i,j,k)=mdKDE([Xpat Ypat Yt],[Xpati(i) Ypati(j) Yti(k)],bw_coeff);            
            end
        end
    end
    toc;

    pdf = pdf./sum(sum(sum(pdf)));    % normalize

    A = pdf;
    B = repmat(sum(pdf,3),[1,1,size(pdf,3)]);
    C = repmat(sum(pdf,1),[size(pdf,1),1,1]);
    D = repmat(sum(sum(pdf,1),3),[size(pdf,1),1,size(pdf,3)]);
    T = A.*log2((A.*D)./(B.*C));
    T = sum(T(~isnan(T) & ~isinf(T)));
    
    TE = [TE,T];
end

