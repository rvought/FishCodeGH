% Load Tracked Data
dataFolderName = '.';

load(fullfile(dataFolderName,'CaveDataRev2018a.mat'));
load(fullfile(dataFolderName,'SurfaceDataRev2018a.mat'));

nCave = length(cave);
nSrf = length(srf);

mmnorm = @(x) ( x - repmat(min(x),size(x,1),1) ) ./ ( repmat(max(x),size(x,1),1) - repmat(min(x),size(x,1),1) );

%% Compile data and create windowed data for TRENTOOL

data.trial = {};
data.time = {};
[meta.dataset,meta.pair,meta.window,meta.df,meta.ddf,meta.dist,meta.ddist] = deal([]);

[dist_all,df_all,ddist_all,ddf_all,ddist_abs_all,ddf_abs_all] = deal(cell(nCave+nSrf,1));

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
                meta.window = [meta.window,w];
                
                meta.dist = [meta.dist,nanmean(DST(:,w))];
                meta.ddist = [meta.ddist,nanmean(DDST(:,w))];
                meta.ddist_abs = [meta.ddist,nanmean(DDST_ABS(:,w))];
                
                meta.df = [meta.df,nanmean(DF(:,w))];
                meta.ddf = [meta.ddf,nanmean(DDF(:,w))];
                meta.ddf_abs = [meta.ddf,nanmean(DDF_ABS(:,w))];
            end
        end
    end
end

%% Comparison of computed variables between cave and srf

%% For each pair - distance
dist_pair_cave = cellfun(@(x) nanmean(x),dist_all(1:nCave),'UniformOutput',false);
dist_pair_cave = [dist_pair_cave{:}];

dist_pair_srf = cellfun(@(x) nanmean(x),dist_all(nCave+1:end),'UniformOutput',false);
dist_pair_srf = [dist_pair_srf{:}];

max_val = max([dist_pair_cave,dist_pair_srf]);
edges = linspace(0,max_val,100);

clf;
hold on;

histogram(dist_pair_cave,edges);
histogram(dist_pair_srf,edges);

legend('Cave','Surface');
xlabel('Distance');
ylabel('Pairs')
title('Distance - cave vs surface');
hold off;

[h_dist,p_dist,ks2stat_dist] = kstest2(dist_pair_cave,dist_pair_srf);
fprintf('\nProbability of %f that the distances are from the same distribution',p_dist);

%% For each pair - df
df_pair_cave = cellfun(@(x) nanmean(x),df_all(1:nCave),'UniformOutput',false);
df_pair_cave = [df_pair_cave{:}];

df_pair_srf = cellfun(@(x) nanmean(x),df_all(nCave+1:end),'UniformOutput',false);
df_pair_srf = [df_pair_srf{:}];

max_val = max([df_pair_cave,df_pair_srf]);
edges = linspace(0,max_val,100);

clf;
hold on;

histogram(df_pair_cave,edges);
histogram(df_pair_srf,edges);

legend('Cave','Surface');
xlabel('Df (Hz)');
ylabel('Pairs')
title('Df - cave vs surface');
hold off;

[h_df,p_df,ks2stat_df] = kstest2(df_pair_cave,df_pair_srf);
fprintf('\nProbability of %f that the dfs are from the same distribution',p_df);

%% For each pair - change in distance
ddist_pair_cave = cellfun(@(x) nanmean(x),ddist_all(1:nCave),'UniformOutput',false);
ddist_pair_cave = [ddist_pair_cave{:}];

ddist_pair_srf = cellfun(@(x) nanmean(x),ddist_all(nCave+1:end),'UniformOutput',false);
ddist_pair_srf = [ddist_pair_srf{:}];

max_val = max([ddist_pair_cave,ddist_pair_srf]);
edges = linspace(0,max_val,100);

clf;
hold on;

histogram(ddist_pair_cave,edges);
histogram(ddist_pair_srf,edges);

legend('Cave','Surface');
xlabel('Change in Distance');
ylabel('Pairs')
title('Change in Distance - cave vs surface');
hold off;

[h_ddist,p_ddist,ks2stat_ddist] = kstest2(ddist_pair_cave,ddist_pair_srf);
fprintf('\nProbability of %f that the change in distances are from the same distribution',p_ddist);

%% For each pair - change in df
ddf_pair_cave = cellfun(@(x) nanmean(x),ddf_all(1:nCave),'UniformOutput',false);
ddf_pair_cave = [ddf_pair_cave{:}];

ddf_pair_srf = cellfun(@(x) nanmean(x),ddf_all(nCave+1:end),'UniformOutput',false);
ddf_pair_srf = [ddf_pair_srf{:}];

max_val = max([ddf_pair_cave,ddf_pair_srf]);
edges = linspace(0,max_val,100);

clf;
hold on;

histogram(ddf_pair_cave,edges);
histogram(ddf_pair_srf,edges);

legend('Cave','Surface');
xlabel('Change in Df (Hz)');
ylabel('Pairs')
title('Change in Df - cave vs surface');
hold off;

[h_ddf,p_ddf,ks2stat_ddf] = kstest2(ddf_pair_cave,ddf_pair_srf);
fprintf('\nProbability of %f that the change in dfs are from the same distribution',p_ddf);

%% For each pair - relative speed
ddist_abs_pair_cave = cellfun(@(x) nanmean(x),ddist_abs_all(1:nCave),'UniformOutput',false);
ddist_abs_pair_cave = [ddist_abs_pair_cave{:}];

ddist_abs_pair_srf = cellfun(@(x) nanmean(x),ddist_abs_all(nCave+1:end),'UniformOutput',false);
ddist_abs_pair_srf = [ddist_abs_pair_srf{:}];

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

[h_ddist_abs,p_ddist_abs,ks2stat_ddist_abs] = kstest2(ddist_abs_pair_cave,ddist_abs_pair_srf);
fprintf('\nProbability of %f that the relative speeds are from the same distribution',p_ddist_abs);

%% For each pair - relative frequency speed
ddf_abs_pair_cave = cellfun(@(x) nanmean(x),ddf_abs_all(1:nCave),'UniformOutput',false);
ddf_abs_pair_cave = [ddf_abs_pair_cave{:}];

ddf_pair_srf = cellfun(@(x) nanmean(x),ddf_abs_all(nCave+1:end),'UniformOutput',false);
ddf_pair_srf = [ddf_pair_srf{:}];

max_val = max([ddf_abs_pair_cave,ddf_pair_srf]);
edges = linspace(0,max_val,100);

clf;
hold on;

histogram(ddf_abs_pair_cave,edges);
histogram(ddf_pair_srf,edges);

legend('Cave','Surface');
xlabel('Relative frequency speed (Hz)');
ylabel('Pairs')
title('Relative frequency speed - cave vs surface');
hold off;

[h_ddf_abs,p_ddf_abs,ks2stat_ddf_abs] = kstest2(ddf_abs_pair_cave,ddf_pair_srf);
fprintf('\nProbability of %f that the relative frequency speeds are from the same distribution',p_ddf_abs);


%% STATS FOR EACH FISH

%% For each fish - distance

[dist_fish_cave,dist_fish_srf] = deal([]);
for j = 1:(nCave+nSrf)
    nFish = nFish_all(j);
  
    if nFish>2
        C = nchoosek(1:nFish,2);

        for f = 1:nFish
            fishIdx = find(sum(C==f,2));

            if j<=nCave
                dist_fish_cave = [dist_fish_cave,nanmean(nanmean(dist_all{j}(:,fishIdx)))];
            else
                dist_fish_srf = [dist_fish_srf,nanmean(nanmean(dist_all{j}(:,fishIdx)))];
            end
        end
    end
end

max_val = max([dist_fish_cave,dist_fish_srf]);
edges = linspace(0,max_val,50);

clf;
hold on;

histogram(dist_fish_cave,edges);
histogram(dist_fish_srf,edges);

legend('Cave','Surface');
xlabel('Distance');
ylabel('Pairs')
title('Distance - cave vs surface');
hold off;

[h_fish_dist,p_fish_dist,ks2stat_fish_dist] = kstest2(dist_fish_cave,dist_fish_srf);
fprintf('\nProbability of %f that the distances are from the same distribution',p_fish_dist);

%% For each fish - df

[df_fish_cave,df_fish_srf] = deal([]);
for j = 1:(nCave+nSrf)
    nFish = nFish_all(j);
  
    if nFish>2
        C = nchoosek(1:nFish,2);

        for f = 1:nFish
            fishIdx = find(sum(C==f,2));

            if j<=nCave
                df_fish_cave = [df_fish_cave,nanmean(nanmean(df_all{j}(:,fishIdx)))];
            else
                df_fish_srf = [df_fish_srf,nanmean(nanmean(df_all{j}(:,fishIdx)))];
            end
        end
    end
end

max_val = max([df_fish_cave,df_fish_srf]);
edges = linspace(0,max_val,50);

clf;
hold on;

histogram(df_fish_cave,edges);
histogram(df_fish_srf,edges);

legend('Cave','Surface');
xlabel('Df');
ylabel('Pairs')
title('Df - cave vs surface');
hold off;

[h_fish_df,p_fish_df,ks2stat_fish_df] = kstest2(df_fish_cave,df_fish_srf);
fprintf('\nProbability of %f that the dfs are from the same distribution',p_fish_df);

%% For each fish - change in distance

[ddist_fish_cave,ddist_fish_srf] = deal([]);
for j = 1:(nCave+nSrf)
    nFish = nFish_all(j);
  
    if nFish>2
        C = nchoosek(1:nFish,2);

        for f = 1:nFish
            fishIdx = find(sum(C==f,2));

            if j<=nCave
                ddist_fish_cave = [ddist_fish_cave,nanmean(nanmean(ddist_all{j}(:,fishIdx)))];
            else
                ddist_fish_srf = [ddist_fish_srf,nanmean(nanmean(ddist_all{j}(:,fishIdx)))];
            end
        end
    end
end

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

[h_fish_ddist,p_fish_ddist,ks2stat_fish_ddist] = kstest2(ddist_fish_cave,ddist_fish_srf);
fprintf('\nProbability of %f that the change in distances are from the same distribution',p_fish_ddist);

%% For each fish - change in df

[ddf_fish_cave,ddf_fish_srf] = deal([]);
for j = 1:(nCave+nSrf)
    nFish = nFish_all(j);
  
    if nFish>2
        C = nchoosek(1:nFish,2);

        for f = 1:nFish
            fishIdx = find(sum(C==f,2));

            if j<=nCave
                ddf_fish_cave = [ddf_fish_cave,nanmean(nanmean(ddf_all{j}(:,fishIdx)))];
            else
                ddf_fish_srf = [ddf_fish_srf,nanmean(nanmean(ddf_all{j}(:,fishIdx)))];
            end
        end
    end
end

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

[h_fish_ddf,p_fish_ddf,ks2stat_fish_ddf] = kstest2(ddf_fish_cave,ddf_fish_srf);
fprintf('\nProbability of %f that the dfs are from the same distribution',p_fish_ddf);

%% For each fish - relative speed

[ddist_abs_fish_cave,ddist_abs_fish_srf] = deal([]);
for j = 1:(nCave+nSrf)
    nFish = nFish_all(j);
  
    if nFish>2
        C = nchoosek(1:nFish,2);

        for f = 1:nFish
            fishIdx = find(sum(C==f,2));

            if j<=nCave
                ddist_abs_fish_cave = [ddist_abs_fish_cave,nanmean(nanmean(ddist_abs_all{j}(:,fishIdx)))];
            else
                ddist_abs_fish_srf = [ddist_abs_fish_srf,nanmean(nanmean(ddist_abs_all{j}(:,fishIdx)))];
            end
        end
    end
end

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

[h_fish_ddist_abs,p_fish_ddist_abs,ks2stat_fish_ddist_abs] = kstest2(ddist_abs_fish_cave,ddist_abs_fish_srf);
fprintf('\nProbability of %f that the change in distances are from the same distribution',p_fish_ddist_abs);

%% For each fish - relative frequency speed

[ddf_abs_fish_cave,ddf_abs_fish_srf] = deal([]);
for j = 1:(nCave+nSrf)
    nFish = nFish_all(j);
  
    if nFish>2
        C = nchoosek(1:nFish,2);

        for f = 1:nFish
            fishIdx = find(sum(C==f,2));

            if j<=nCave
                ddf_abs_fish_cave = [ddf_abs_fish_cave,nanmean(nanmean(ddf_abs_all{j}(:,fishIdx)))];
            else
                ddf_abs_fish_srf = [ddf_abs_fish_srf,nanmean(nanmean(ddf_abs_all{j}(:,fishIdx)))];
            end
        end
    end
end

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

[h_fish_ddf_abs,p_fish_ddf_abs,ks2stat_fish_ddf_abs] = kstest2(ddf_abs_fish_cave,ddf_abs_fish_srf);
fprintf('\nProbability of %f that the relative frequency speeds are from the same distribution',p_fish_ddf_abs);

%% Distance vs df plots

clf, hold on;
plot(ddist_fish_cave,ddf_fish_cave,'.b')
plot(ddist_fish_srf,ddf_fish_srf,'.r')

grid on;

hold off;

%%
clf, hold on;

plot(ddist_abs_fish_cave,ddf_abs_fish_cave,'.b')
plot(ddist_abs_fish_srf,ddf_abs_fish_srf,'.r')


hold off;

%% Run TRENTOOL analysis

%%
resultFolderName = '~/Google Drive/trentool_results';

load(fullfile(resultFolderName,'all_channels_windowed.mat'));
trentool_result = [trentool_result{:}];


%%

TEprepare = [trentool_result.TEprepare];
u_in_ms = [TEprepare.u_in_ms];
u_in_ms = u_in_ms(1,:);
n_u = length(u_in_ms);

channelPairs = TEprepare(1).channelcombilabel;

% channelPairs = channelPairs(:,1:2)'; % REMOVE LATER
nChannelPairs = size(channelPairs,1);

nWindows =  size(TEprepare(1).ACT,3);

TE = nan(n_u,nChannelPairs,nWindows);

for u = 1:n_u
    for c = 1:nChannelPairs
        TE(u,c,TEprepare(u).trials{c,1}) = trentool_result(u).TEmat(c,1:TEprepare(u).nrtrials(c,1));
        
%         TE(u,c,TEprepare(u).trials{c,1}) = TEpermtest.TEmat_sur(c,1:TEprepare(u).nrtrials(c,1));
    end
end


channelLabels = cell(nChannelPairs,1);
for k = 1:nChannelPairs
    channelLabels{k} = [channelPairs{k,1},' -> ',channelPairs{k,2}];
end

%%


% TEmax = zeros(


%%

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

%%
plot(squeeze(max(TE(:,2,idx))));


% plot(t,TE(idx));
% plot(t,TE2(idx));
% plot(t,mmnorm(meta.sens(idx)));
% plot(t(1:end-1),mmnorm(abs(diff(meta.dist(idx)))./abs(diff(meta.df(idx)))));
% legend('Dist','Df','TE','TE2')
hold off;

%%

act = squeeze(TGA_results.ACT.actvalue(:,2,:));
teIdx = act>=0 & act<=100;
TE = NaN(size(teIdx));
TE(teIdx) = TGA_results.TEmat;


% act = squeeze(TGA_results2.ACT.actvalue(:,2,:));
% teIdx = act>=0 & act<=100;
% TE2 = NaN(size(teIdx));
% TE2(teIdx) = TGA_results2.TEmat;

d = 4;
p = 9;

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
% plot(t,TE2(idx));
% plot(t,mmnorm(meta.sens(idx)));
% plot(t(1:end-1),mmnorm(abs(diff(meta.dist(idx)))./abs(diff(meta.df(idx)))));
legend('Dist','Df','TE','TE2')
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
%         
% %         subplot(3,1,1), hold on;
%             plot(freq,te_fish,'.','Color',col);
%         subplot(3,1,2), hold on;
%             plot(freq,te_fish_var,'.','Color',col);
%         subplot(3,1,3), hold on;
%             plot(freq,te_fish./te_fish_var,'.','Color',col);
%             plot(te_fish,te_fish_var,'.','Color',col);

%         plot(freq,vel,'.','Color',col);
%         plot(vel,te_fish,'.','Color',col);
%         plot(vel,te_fish_var,'.','Color',col);
%         plot(vel,te_fish./te_fish_var,'.','Color',col);
    end
end

% subplot(3,1,1), grid on;
% subplot(3,1,2), grid on;
% subplot(3,1,3), grid on;
grid on;
hold off;