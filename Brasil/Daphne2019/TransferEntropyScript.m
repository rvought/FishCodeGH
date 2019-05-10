%% Load Tracked Data
% dataFolderName = '~/Downloads';
% % dataFolderName = '/Users/Shared/Data/Brasil';
% load(fullfile(dataFolderName,'CaveDataRev2018a.mat'));
% load(fullfile(dataFolderName,'SurfaceDataRev2018a.mat'));

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



%% Does this go here?
[dist_fish_cave,dist_fish_srf,...
df_fish_cave,df_fish_srf,...
ddist_fish_cave,ddist_fish_srf,...
ddf_fish_cave,ddf_fish_srf,...
ddist_abs_fish_cave,ddist_abs_fish_srf,...
ddf_abs_fish_cave,ddf_abs_fish_srf] = deal([]);








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



[TEmax,TEmaxIdx] = max(TE,[],1);

TEmax = squeeze(TEmax)';
TEmaxIdx = squeeze(TEmaxIdx)';

idx = TEmax<=0;
TEmax(idx) = NaN;
TEmaxIdx(idx) = NaN;


%% Examples of high-distance, high TE interactions
idx = find(TEmax(:,3) > prctile(TEmax(:,3),85) & meta.dist' > prctile(meta.dist,85));

for k = 1:length(idx)
    clf;
    idx(k)
    subplot(3,1,1);
    plot(data.trial{idx(k)}(1,:), '.-');
    title('Distance (cm)');
    
    subplot(3,1,2);
    plot(data.trial{idx(k)}(2,:), '.-');
    title('Df (Hz)');
    
    subplot(3,1,3); 
    yyaxis left;
    plot(data.trial{idx(k)}(1,:), 'b.-');
    yyaxis right; hold on;
    plot(data.trial{idx(k)}(2,:), 'r.-');    
    title('Df (Hz)');
    
    
    pause;
end

%% Our paper examples 

% Positive (unexpected) correlation
figure(28); clf; yyaxis left; plot(data.trial{944}(1,:), 'b.-'); yyaxis right; plot(data.trial{944}(2,:), 'r.-');
figure(28); yyaxis left; ylim([60 200]);
figure(28); yyaxis right; ylim([13 16]);

% Negative correlation, at a SPOOKY distance
figure(27); clf; yyaxis left; plot(data.trial{879}(1,:), 'b.-'); yyaxis right; plot(data.trial{879}(2,:), 'r.-');
figure(27); yyaxis left; ylim([60 200]);
figure(27); yyaxis right; ylim([83.5 86.5]);




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
    
    plot(dist_fish_cave, TE_fish_cave(:,k),'.b');
    plot(dist_fish_srf, TE_fish_srf(:,k),'.r');
    
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
