function [foo, out] = analretent(in, rango)
% Usage [foo, out] = analretent(in)
% Where "in" is the reduced structure from catfish Happiness


for kk = length(in):-1:1
    tims(kk) = max(in(kk).tt);
end
mintim = max(tims);

if nargin == 1; rango = [0, mintim]; end;


%% Make Colors
    numcolors = 100;
     oclrs = hsv(numcolors); % Colors for plotting up to numcolors different fishes
%     shuff = [11 22 2 7 13 18 1 23 16 12 10 17 3 24 8 14 19 20 4 16 25 5 15 9 21];
     shuff = randperm(numcolors);
     clrs = zeros(numcolors, 3);
     for i=1:numcolors; clrs(i,:) = oclrs(shuff(i),:); end; 
     
    
figure(1); clf;
figure(2); clf; hold on;
figure(3); clf; hold on;

%% Make some plots
for j=1:min([numcolors length(in)])
    
    if length(find(in(j).tt > rango(1) & in(j).tt < rango(2))) > 1
    
        txy = find(in(j).tt > rango(1) & in(j).tt < rango(2));
        ttim = find(in(j).tim > rango(1) & in(j).tim < rango(2));
        
        figure(1); ax(j) = subplot(5,6,j); plot(in(j).tx(txy), in(j).ty(txy), '*', 'Color', clrs(j,:));
    
        figure(2); plot(in(j).tim(ttim), in(j).freq(ttim), '*', 'Color', clrs(j,:));
    
        figure(3); plot(in(j).tx(txy), in(j).ty(txy), '*', 'Color', clrs(j,:));
    
    end
    
end

linkaxes(ax, 'xy'); figure(1); subplot(5,6,1); axis([-200, 200, -200, 200]);

%% Let's do heat maps!

for jj = 1:length(in)
   
    stepsize = 10;
    hotness = zeros(40,40);

    for xx = -200:stepsize:200-stepsize
        for yy = -200:stepsize:200-stepsize
            
            hotness((210+xx)/stepsize,(210+yy)/stepsize) = hotness((210+xx)/stepsize,(210+yy)/stepsize) + length(find(in(jj).tx > xx & in(jj).tx < xx+stepsize & in(jj).ty > yy & in(jj).ty < yy+stepsize));        
        
        end
    end
    
    out(jj).hotness = hotness;
    
    figure(4); xa(jj) = subplot(5,6,jj); surf(1:40, 1:40, hotness); view(0,90); caxis([0 200]); colormap('HOT');
    
end
linkaxes(xa, 'xy'); 

%% Calculate distances to other fish

    pw = nchoosek(1:length(in), 2);
    for kk = 1:length(pw)
        
        X(:,1) = in(pw(kk,1)).x;
        X(:,2) = in(pw(kk,1)).y;
        Y(:,1) = in(pw(kk,2)).x;
        Y(:,2) = in(pw(kk,2)).y;

        asdf = pdist2(X,Y);
        out(kk).pair = pw(kk,:);
        for j=1:length(asdf)
            out(kk).pairdist(j) = asdf(j,j);
        end
    end
    
    pairlist = [out.pair];
    
    foo.pp(1,:) = 1:length(in)-1; 
    for kk = 2:length(in); foo.pp(kk,:) = sort([find(pairlist(:,1) == kk)' find(pairlist(:,2) == kk)']); end;

  %  >> figure(27); clf; hold on; for i=1:length(pp); plot(out(pp(i)).pairdist); end;                   

    
    

%% Get the frequency data and analyze!

load caves10data.mat	
load caves15data.mat	
load caves18data.mat	
load caves20data.mat	
load caves24data.mat	
load caves27data.mat	
load caves2data.mat	
load caves32data.mat	
load caves3data.mat	
load caves7data.mat
load caves12data.mat	
load caves16data.mat	
load caves19data.mat	
load caves21data.mat	
load caves25data.mat	
load caves28data.mat	
load caves30data.mat	
load caves33data.mat	
load caves4data.mat	
load caves8data.mat
load caves14data.mat	
load caves17data.mat	
load caves1data.mat	
load caves23data.mat	
load caves26data.mat	
load caves29data.mat	
load caves31data.mat	
load caves34data.mat	
load caves5data.mat	
vars = {'s5'; 's34'; 's31'; 's29'; 's26'; 's23'; 's1'; 's17'; 's14'; 's8'; ...
    's4'; 's33'; 's30'; 's28'; 's25'; 's21'; 's19'; 's16'; 's12'; 's7'; 's3'; 's32'; 's2'; 's27'; 's24'; 's20'; 's18'; 's15'; 's10'};


% Frequency versus std/var data

foo.fmean = [];
foo.fstd = [];
foo.fvar = [];
for kk = 1:length(vars)
    
        aa = eval(vars{kk,:});
    
    for jj = 1:length(aa)
        
        foo.fmean(end+1) = mean([aa(jj).freq]);
        foo.fstd(end+1) = std([aa(jj).freq]);
        foo.fvar(end+1) = var([aa(jj).freq]);
        
    end
    
end

foo.freqVstd = fitlm([foo.fmean], [foo.fstd]);
foo.freqVvar = fitlm([foo.fmean], [foo.fvar]);


