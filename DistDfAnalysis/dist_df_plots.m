load results/all_channels_windowed.mat
trentool_result = [trentool_result{:}];

load CaveDataRev2018a.mat
load SurfaceDataRev2018a.mat

mmnorm = @(x) ( x - repmat(min(x),size(x,1),1) ) ./ ( repmat(max(x),size(x,1),1) - repmat(min(x),size(x,1),1) );


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
        TE(u,c,TEprepare(u).trials{c,1}) = trentool_result(1).TEmat(c,1:TEprepare(u).nrtrials(c,1));
    end
end

%%

d = 3;
p = 6;

idx = meta.dataset==d & meta.pair==p;

clf, hold on
plot(mmnorm(meta.dist(idx)));
plot(mmnorm(meta.df(idx)));

plot(squeeze(max(TE(:,1,idx))));
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