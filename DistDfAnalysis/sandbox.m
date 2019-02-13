load CaveDataRev2018a.mat
load SurfaceDataRev2018a.mat
%%

addpath('MIDERv2');

%%
nCaveTrials = length(cave);
nSurfaceTrials = length(srf);

%% Plot example data

% c = 3;
clf, hold on;
for c = 1:nCaveTrials
    nFish = length(cave(c).fish);
    if nFish>1
        nTimes = length(cave(c).t);
        col = jet(nFish);
        % clf
        % for f = 1:nFish
        %     subplot(1,2,1), hold on;
        %     plot3(cave(c).fish(f).x,cave(c).fish(f).y,cave(c).fish(f).z,'Color',col(f,:));
        %     hold off;
        %     
        %     subplot(1,2,2), hold on;
        %     plot(cave(c).t,cave(c).fish(f).freq(:,2),'Color',col(f,:));
        %     hold off;
        % end
    
        % Pairwise states

        C = nchoosek(1:nFish,2);
        nPairs = size(C,1);

        [dist,df] = deal(zeros(nTimes,nPairs));
        for k = 1:nPairs
            dist(:,k) = sqrt((cave(c).fish(C(k,1)).x - cave(c).fish(C(k,2)).x).^2 + (cave(c).fish(C(k,1)).y - cave(c).fish(C(k,2)).y).^2);
            df(:,k) = cave(c).fish(C(k,1)).freq(:,2) - cave(c).fish(C(k,2)).freq(:,2); 
        end

        mmnorm = @(x) ( x - repmat(min(x),size(x,1),1) ) ./ ( repmat(max(x),size(x,1),1) - repmat(min(x),size(x,1),1) );

        ddist = mmnorm(smoothdata(diff(smoothdata(dist,'movmean',100)),'movmean',10));
        ddf = mmnorm(smoothdata(diff(smoothdata(df,'movmean',100)),'movmean',10 ));
%         clf;
        [R,L] = deal(zeros(nPairs,1));
        for k = 1:nPairs
%             r = corrcoef(ddist(:,k),ddf(:,k),'Rows','complete');
%             R(k) = r(1,2);

            [r,lag] = xcorr(ddist(:,k),ddf(:,k),'coeff');
            [R(k),maxIdx] = max(abs(r));
            L(k) = lag(maxIdx);

        %     plot(r), pause;

        %     clf, hold on;
        %     plot(ddist(:,k),'b');
        %     plot(ddf(:,k),'r');
        %     hold off;
        %     pause;

        end

        plot(mean(dist),abs(R),'.')
    end
end

%% Linear response matrix

[G_dist,G_df] = deal(zeros(nPairs,nPairs,nTimes-1));

F = [cave(c).fish.x];
Y = [cave(c).fish.y];
F = [cave(c).fish.freq];
F = F(:,2:2:end);

GX = diff(log(F));
GY = diff(log(Y));
GF = diff(log(F));

GXX = abs( GX(:,C(:,1)) ./ GX(:,C(:,2)) );
GXY = abs( GX(:,C(:,1)) ./ GY(:,C(:,2)) );
GXF = abs( GX(:,C(:,1)) ./ GF(:,C(:,2)) );
GYY = abs( GY(:,C(:,1)) ./ GY(:,C(:,2)) );
GYF = abs( GY(:,C(:,1)) ./ GF(:,C(:,2)) );
GFF = abs( GF(:,C(:,1)) ./ GF(:,C(:,2)) );

%% As a function of distance, how much does a change in distance drive a change in df?

c = 6;
nFish = length(cave(c).fish);
nTimes = length(cave(c).t);

C = nchoosek(1:nFish,2);
nPairs = size(C,1);

[dist,df] = deal(zeros(nTimes,nPairs));
for k = 1:nPairs
    dist(:,k) = sqrt((cave(c).fish(C(k,1)).x - cave(c).fish(C(k,2)).x).^2 + (cave(c).fish(C(k,1)).y - cave(c).fish(C(k,2)).y).^2);
    df(:,k) = cave(c).fish(C(k,1)).freq(:,2) - cave(c).fish(C(k,2)).freq(:,2); 
end

ddist = smoothdata(diff(smoothdata(dist,'movmean',100)),'movmean',10);
ddf = smoothdata(diff(smoothdata(df,'movmean',100)),'movmean',10);

G = abs(ddf./ddist);
idx = ~isnan(G) & ~isinf(G) & G~=0;

G_dist = dist(1:end-1,:);

plot(G_dist(idx),G(idx),'.')
ylim([0,0.1])

distEdges = 0:1:250;
sensEdges = 0:0.0005:0.02;
N = hist3([G_dist(idx),G(idx)],'Edges',{distEdges,sensEdges});
imagesc(distEdges,sensEdges,rot90(flipud(N),-1))
set(gca,'YDir','normal');

%% Use MIDER for network inference of a sample trial

c = 5;
nFish = length(cave(c).fish);
f = [cave(c).fish.freq];
f = f(:,2:2:end);

df_norm = [cave(c).fish.x];
dist_norm = [cave(c).fish.y];
z = [cave(c).fish.z];
v = [zeros(1,nFish);sqrt(diff(df_norm).^2 + diff(dist_norm).^2)];
d = [zeros(1,nFish);diff(f)];


variablesF = mat2cell([repmat('F',nFish,1),num2str((1:nFish)')],ones(1,nFish));
variablesX = mat2cell([repmat('X',nFish,1),num2str((1:nFish)')],ones(1,nFish));
variablesY = mat2cell([repmat('Y',nFish,1),num2str((1:nFish)')],ones(1,nFish));
variablesZ = mat2cell([repmat('Z',nFish,1),num2str((1:nFish)')],ones(1,nFish));
variablesV = mat2cell([repmat('V',nFish,1),num2str((1:nFish)')],ones(1,nFish));
variablesD = mat2cell([repmat('D',nFish,1),num2str((1:nFish)')],ones(1,nFish));

variables = [variablesD];

% Run MIDER on all data
% options = [];
% miderOut = runMIDER(x,variables,options);

% Run MIDER on chunked data
nData = size(f,1);

nChunks =  2;
nWindow = ceil(nData/nChunks);

percOverlap = 0;
nOverlap = floor(nWindow*percOverlap/100);

[F,X,Y,Z,V,D] = deal([]);
for k = 1:size(f,2)
    F(:,:,k) = buffer(f(:,k),nWindow,nOverlap);
    X(:,:,k) = buffer(df_norm(:,k),nWindow,nOverlap);
    Y(:,:,k) = buffer(dist_norm(:,k),nWindow,nOverlap);
    Z(:,:,k) = buffer(z(:,k),nWindow,nOverlap);
    V(:,:,k) = buffer(v(:,k),nWindow,nOverlap);
    D(:,:,k) = buffer(d(:,k),nWindow,nOverlap);
end

nActualWindows = size(F,2);
options = [];
for k = 1:nActualWindows
    [miderOut,optionsOut] = runMIDER([squeeze(D(:,k,:))],variables,options);
    
    clf, hold on;
    
    xWin = mean(squeeze(X(:,k,:)));
    yWin = mean(squeeze(Y(:,k,:)));
    zWin = mean(squeeze(Z(:,k,:)));
    
    con_array = miderOut.con_array;
    crange = [min(con_array(con_array>0)),max(con_array(:))];
    
    plot(xWin,yWin,'o');
    for i=1:nFish
        for j=1:nFish
            if con_array(i,j) > 0
                plot([xWin(i),xWin(j)],[yWin(i),yWin(j)],'-','Color',vals2colormap(con_array(i,j),'jet',crange));
            end
        end
    end
    
    % Second, plot data points and names:
    text(xWin+0.05,yWin,variables,'FontSize',14)
                
   
%     for i=1:nFish
%         for j=1:nFish
%             if con_array (i,j) > 0 
%                 if miderOut.T(i,j) > 0 %Output.T(j,i)
%                     arrow([xWin(i),yWin(i)], [xWin(j),yWin(j)],...
%                         'Width', ,'Length',90,'BaseAngle',40,'TipAngle',30)  
%                 else
%                     if miderOut.T(j,i) > 0
%                         arrow([xWin(j),yWin(j)], [xWin(i),yWin(i)],...
%                         'Width', abs(50*(con_array(i,j))^1.5),'Length',90,'BaseAngle',40,'TipAngle',30)
%                     else
%                         arrow([xWin(j),yWin(j)], [xWin(i),yWin(i)],...
%                         'Width', abs(50*(con_array(i,j))^1.5),'Length',0)
%                     end
%                 end
%             end
%         end
%     end
    
%     plotResults(miderOut,nFish,variables,optionsOut)
    hold off;
    pause;
end

%% For each pair of fish, what is the mutual information between pairwise distance and df?
% For each part of the trial
% As a function of distance

c = 5;
nFish = length(cave(c).fish);
nTimes = length(cave(c).t);
col = jet(nFish);

% Pairwise states
C = nchoosek(1:nFish,2);
nPairs = size(C,1);

[dist,df] = deal(zeros(nTimes,nPairs));
for k = 1:nPairs
    dist(:,k) = sqrt((cave(c).fish(C(k,1)).x - cave(c).fish(C(k,2)).x).^2 + (cave(c).fish(C(k,1)).y - cave(c).fish(C(k,2)).y).^2);
    df(:,k) = cave(c).fish(C(k,1)).freq(:,2) - cave(c).fish(C(k,2)).freq(:,2); 
end

mmnorm = @(x) ( x - repmat(min(x),size(x,1),1) ) ./ ( repmat(max(x),size(x,1),1) - repmat(min(x),size(x,1),1) );
dist_smooth = smoothdata(dist,'movmean',100);
df_smooth = smoothdata(df,'movmean',100);

p_dist = histcounts(dist_smooth,0:5:300,'Normalization','probability');
p_df = histcounts(df_smooth,0:5:200,'Normalization','probability');
p_dist_df = histcounts2(dist_smooth,df_smooth,0:5:300,0:5:200,'Normalization','probability');

% Window data
nData = size(dist_smooth,1);
nChunks =  20;
nWindow = ceil(nData/nChunks);
percOverlap = 90;
nOverlap = floor(nWindow*percOverlap/100);

[DF,DST] = deal([]);
for k = 1:nPairs
%     clf
%     hold on
%     plot(mmnorm(dist_smooth(:,k)));
%     plot(mmnorm(df_smooth(:,k)));
%     hold off;
%     pause;
    
    DST(:,:,k) = buffer(dist_smooth(:,k),nWindow,nOverlap,'nodelay');
    DF(:,:,k) = buffer(df_smooth(:,k),nWindow,nOverlap,'nodelay');
end
nActualWindows = size(DF,2);


R = zeros(nActualWindows,nPairs);
for k = 1:nPairs
    for j = 1:nActualWindows
%         clf
%         hold on;
%         plot(mmnorm(DST(:,j,k)));
%         plot(mmnorm(DF(:,j,k)));
%         hold off;
%         pause;
        
        r = corrcoef(DST(:,j,k),DF(:,j,k),'Rows','complete');
        R(j,k) = r(1,2);
    end
end

%%

mean_dist = squeeze(mean(DST,1));
mean_df = squeeze(mean(DF,1));

clf, hold on;
col = jet(nFish);
for k = 1%nFish
    idx = find(C(:,1)==k | C(:,2)==k);
    subplot(2,2,1);
    Rfish = R(:,idx);
    
    plot(Rfish);%,'.','Color',col(k,:));
    
    subplot(2,2,2);
    plot(mean_dist(:,idx));
%     histogram(Rfish(:));

    subplot(2,2,3);
    plot(1./sqrt(mean_dist(:,idx)),R(:,idx),'.')
    
    subplot(2,2,4);
    plot(1./sqrt(mean_dist(:,idx)))
end

%%

crange = [-1,1];

for j = 1:nActualWindows-1
    clf, hold on;
    for k = 1:nPairs
        rectangle('Position',[C(k,1),C(k,2),1,1],'FaceColor',vals2colormap(R(j,k),'jet',crange))
    end
    hold off;
    pause;
end

%%
crange = [-1,1];

clf, hold on;
for j = 1:20%nActualWindows-1
    for k = 1:nPairs
        plotcube([1,1,1],[C(k,1),C(k,2),j],0.2,vals2colormap(R(j,k),'jet',crange))
    end
end
hold off;

%% Showing that correlations are not chance - use shuffled distribution

N = 10000;
R_shuffled = zeros(N,1);

pairSample1 = randsample(nPairs,N,true);
windowSample1 = randsample(nActualWindows,N,true);

pairSample2 = randsample(nPairs,N,true);
windowSample2 = randsample(nActualWindows,N,true);

for n = 1:N
    r = corrcoef(DST(:,windowSample1(n),pairSample1(n)),DF(:,windowSample2(n),pairSample2(n)),'Rows','complete');
    R_shuffled(n) = r(1,2);
end

clf, hold on;

histogram(R(:),30,'Normalization','probability')
histogram(R_shuffled,30,'Normalization','probability')
hold off;


%%

clf;
for k = 1:nPairs
    df_norm = mmnorm(mean_df(:,k));
    dist_norm = mmnorm(mean_dist(:,k));
    
    subplot(2,1,1), cla, hold on; 
    plot(df_norm)
    plot(dist_norm)
    hold off;
    
    r = corrcoef(df_norm,dist_norm);
    r(1,2)
    if r(1,2)<0
        dist_norm = mmnorm(-dist_norm);
    end
    
    [r,lag] = xcorr(df_norm,dist_norm);
    [~,maxidx] = max(abs(r));
    lag = lag(maxidx)
    
    subplot(2,1,2), cla, hold on; 
    if lag>0    
        plot(mmnorm(df_norm(floor(lag/2)+1:end)))
        plot(mmnorm(dist_norm(1:end-floor(lag/2))))
        plot(mmnorm(dist_norm));
    elseif lag<0    
        plot(mmnorm(df_norm(1:end-floor(-lag/2))))
        plot(mmnorm(dist_norm(floor(-lag/2)+1:end)))
    else
        plot(df_norm)
        plot(dist_norm)
    end
    hold off;
    
    
    pause;
end

%%

[warp,i_df,i_dist] = deal(cell(nPairs,1));
f = 3;


for k = 1:nPairs
    try
        clf;
    %     if C(k,1)==f || C(k,2)==f
            df_norm = mmnorm(mean_df(1:nActualWindows-1,k));
            dist_norm = mmnorm(mean_dist(1:nActualWindows-1,k));

            fprintf('\n\nPair %d',k);
            subplot(4,1,1), hold on;
            title(sprintf('Pair %d: Fish %d-Fish%d',k,C(k,1),C(k,2)));
            plot(df_norm);
            plot(dist_norm);
            hold off;

            r = corrcoef(df_norm,dist_norm);
            r = r(1,2);
            fprintf('\nCorrelation coeff (Raw): %.4f',r);
            if r<0
                dist_norm = mmnorm(-dist_norm);
            end

            [r,lag] = xcorr(df_norm,dist_norm);
            [~,maxidx] = max(abs(r));
            lag = lag(maxidx);

            if lag>0    
                df_norm = mmnorm(df_norm(floor(lag/2)+1:end));
                dist_norm = mmnorm(dist_norm(1:end-floor(lag/2)));
            elseif lag<0    
                df_norm = mmnorm(df_norm(1:end-floor(-lag/2)));
                dist_norm = mmnorm(dist_norm(floor(-lag/2)+1:end));
            end

            r = corrcoef(df_norm,dist_norm);
            r = r(1,2);
            fprintf('\nCorrelation coeff (Post-%d lag): %.4f',lag,r);

            subplot(4,1,2), hold on;
            plot(df_norm);
            plot(dist_norm);
            hold off;

            [~,i_df{k},i_dist{k}] = dtw(df_norm,dist_norm,20);
            warp{k} = i_df{k}-i_dist{k};

            subplot(4,1,3), hold on;
            plot(df_norm(i_df{k}));
            plot(dist_norm(i_dist{k}));
            hold off;

            r = corrcoef(df_norm(i_df{k}),dist_norm(i_dist{k}));
            r = r(1,2);
            fprintf('\nCorrelation coeff (Post-DTW): %.4f',r);

            subplot(4,1,4), hold on;
            plot(i_df{k},i_df{k}-i_dist{k});
            hold off;
    %     end
    catch
    end
pause;    
end

%%

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

mmnorm = @(x) ( x - repmat(min(x),size(x,1),1) ) ./ ( repmat(max(x),size(x,1),1) - repmat(min(x),size(x,1),1) );
dt = mean(diff(cave(1).t));
winInterval = 30; % seconds
winLength = round(winInterval/dt);

mean_dist = highpass(smoothdata(dist,'movmean',winLength,'omitnan'),0.001,1/dt,'ImpulseResponse','iir');
mean_df = highpass(smoothdata(df,'movmean',winLength,'omitnan'),0.001,1/dt,'ImpulseResponse','iir');


%%

for k = 1:nPairs    
    try
        clf;
    %     if C(k,1)==f || C(k,2)==f
            df_norm = mmnorm(mean_df(:,k));
            dist_norm = mmnorm(mean_dist(:,k));

            fprintf('\n\nPair %d',k);
            subplot(4,1,1), hold on;
            title(sprintf('Pair %d: Fish %d-Fish%d',k,C(k,1),C(k,2)));
            plot(df_norm);
            plot(dist_norm);
            hold off;

            r = corrcoef(df_norm,dist_norm);
            r = r(1,2);
            fprintf('\nCorrelation coeff (Raw): %.4f',r);
            if r<0
                dist_norm = mmnorm(-dist_norm);
            end

            [r,lag] = xcorr(df_norm,dist_norm);
            [~,maxidx] = max(abs(r));
            lag = lag(maxidx);

            if lag>0    
                df_norm = mmnorm(df_norm(floor(lag/2)+1:end));
                dist_norm = mmnorm(dist_norm(1:end-floor(lag/2)));
            elseif lag<0    
                df_norm = mmnorm(df_norm(1:end-floor(-lag/2)));
                dist_norm = mmnorm(dist_norm(floor(-lag/2)+1:end));
            end

            r = corrcoef(df_norm,dist_norm);
            r = r(1,2);
            fprintf('\nCorrelation coeff (Post-%d lag): %.4f',lag,r);

            subplot(4,1,2), hold on;
            plot(df_norm);
            plot(dist_norm);
            hold off;

            [~,i_df{k},i_dist{k}] = dtw(df_norm,dist_norm,round(60/dt));
            warp{k} = i_df{k}-i_dist{k};

            subplot(4,1,3), hold on;
            plot(df_norm(i_df{k}));
            plot(dist_norm(i_dist{k}));
            hold off;

            r = corrcoef(df_norm(i_df{k}),dist_norm(i_dist{k}));
            r = r(1,2);
            fprintf('\nCorrelation coeff (Post-DTW): %.4f',r);

            subplot(4,1,4), hold on;
            plot(i_df{k},i_df{k}-i_dist{k});
            hold off;
    %     end
    catch
    end
pause; 
end

%%

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

mmnorm = @(x) ( x - repmat(min(x),size(x,1),1) ) ./ ( repmat(max(x),size(x,1),1) - repmat(min(x),size(x,1),1) );
dt = mean(diff(cave(1).t));
winInterval = 15; % seconds
winLength = round(winInterval/dt);

mean_dist = highpass(smoothdata(dist,'movmean',winLength,'omitnan'),0.001,1/dt,'ImpulseResponse','iir');
mean_df = highpass(smoothdata(df,'movmean',winLength,'omitnan'),0.001,1/dt,'ImpulseResponse','iir');



%%

nWindow = round(round(60/dt)/2)*2;

pairLags = cell(nPairs,1);

for k = 1:nPairs    
%     try
        clf;
    %     if C(k,1)==f || C(k,2)==f
            df_norm = mmnorm(mean_df(:,k));
            dist_norm = mmnorm(mean_dist(:,k));
            
            r = corrcoef(df_norm,dist_norm);
            r = r(1,2);
            fprintf('\nCorrelation coeff (Raw): %.4f',r);
            if r<0
                dist_norm = mmnorm(-dist_norm);
            end
            
            fprintf('\n\nPair %d',k);
            subplot(3,1,1), hold on;
            title(sprintf('Pair %d: Fish %d-Fish%d',k,C(k,1),C(k,2)));
            plot(t,df_norm);
            plot(t,dist_norm);
            xlim([t(1),t(end)]);
            grid on;
            hold off;

            DST = mmnorm(buffer(dist_norm,nWindow,nWindow-1,'nodelay'));
            DF = mmnorm(buffer(df_norm,nWindow,nWindow-1,'nodelay'));
            T = buffer(t,nWindow,nWindow-1,'nodelay');
            
            [lag,maxr] = deal(zeros(size(DST,2),1));
            for j = 1:size(DST,2)
                [r,l] = xcorr(DST(:,j),DF(:,j),'unbiased',round(30/dt));
%                   [r,l,bounds] = crosscorr(DST(:,j),DF(:,j));
                
%                 clf;
%                 subplot(2,1,1), hold on;
%                 plot(T(:,j),DST(:,j),'b');
%                 plot(T(:,j),DF(:,j),'r');
%                 hold off;
%                 
%                 subplot(2,1,2), hold on;
%                 plot(l*dt,r,'k');
%                 hold off;
                
                [maxr(j),maxidx] = max(r);
                lag(j) = l(maxidx)*dt;
                
%                 lag(j)
%                 pause;
            end
            pairLags{k} = lag;
            
            subplot(3,1,2), hold on;
%             lag(abs(lag)>45) = NaN;
%             lag(maxr<0.4) = NaN;
            plot(mean(T),lag)
            xlim([t(1),t(end)]);
            grid on;
            hold off;
            
            subplot(3,1,3), hold on;
%             lag(abs(lag)>45) = NaN;
%             lag(maxr<0.4) = NaN;
            plot(mean(T),maxr)
            xlim([t(1),t(end)]);
            grid on;
            hold off;
            
%     end
pause;
end
pairLags = [horzcat(pairLags{:})];

%%

crange = [-1,1];
nWindows = size(pairLags,1);

for j = 1:nWindows
    clf, hold on;
    for k = 1:nPairs
        rectangle('Position',[C(k,1),C(k,2),1,1],'FaceColor',vals2colormap(R(j,k),'jet',crange))
    end
    hold off;
    pause;
end