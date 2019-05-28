function [fish, sa, sq, spadat, tubular] = Miranda(cave, srf)



epochdur = 300;

%% cave 1, 1 fish, 1000 seconds

%% cave 2, 1 fish, 900 seconds

%% cave 3, 7 fish, good recordings of all fish through over 1200 seconds

currfishnum = 3; maxdur = 1200;

    [tmp, ~] = SpaceCorps(cave(currfishnum), 1, 1:length(cave(currfishnum).fish), [0 epochdur]);
    [~, cvQ] = overture(tmp);
    caveReal = cvQ.realoverlaps;
    caveRand = cvQ.randoverlaps;
    
    for j=2:floor(maxdur/epochdur)
       
        [tmp, ~] = SpaceCorps(cave(currfishnum), 1, 1:length(cave(currfishnum).fish), [(j-1)*epochdur, j*epochdur]);
        [~, cvQ] = overture(tmp);
        caveReal = caveReal + cvQ.realoverlaps;
        caveRand = caveRand + cvQ.randoverlaps;
        
    end



%% cave 4, 7 fish, good recordings of all fish through over 1200 seconds

currfishnum = 4; maxdur = 1200;

    for j=1:floor(maxdur/epochdur)
       
        [tmp, ~] = SpaceCorps(cave(currfishnum), 1, 1:length(cave(currfishnum).fish), [(j-1)*epochdur, j*epochdur]);
        [~, cvQ] = overture(tmp);
        caveReal = caveReal + cvQ.realoverlaps;
        caveRand = caveRand + cvQ.randoverlaps;
        
    end

%% cave 5, 7 fish, good recordings of all fish through over 1200 seconds

currfishnum = 5; maxdur = 1200;

    for j=1:floor(maxdur/epochdur)
       
        [tmp, ~] = SpaceCorps(cave(currfishnum), 1, 1:length(cave(currfishnum).fish), [(j-1)*epochdur, j*epochdur]);
        [~, cvQ] = overture(tmp);
        caveReal = caveReal + cvQ.realoverlaps;
        caveRand = caveRand + cvQ.randoverlaps;
        
    end

%% cave 6, 8 fish, good recordings of all fish through over 1200 seconds LOTS OF SOCIAL

currfishnum = 6; maxdur = 1200;

    for j=1:floor(maxdur/epochdur)
       
        [tmp, ~] = SpaceCorps(cave(currfishnum), 1, 1:length(cave(currfishnum).fish), [(j-1)*epochdur, j*epochdur]);
        [~, cvQ] = overture(tmp);
        caveReal = caveReal + cvQ.realoverlaps;
        caveRand = caveRand + cvQ.randoverlaps;
        
    end


%% cave 7, 7 fish, some gaps over 900 seconds

currfishnum = 7; maxdur = 900;

    for j=1:floor(maxdur/epochdur)
       
        [tmp, ~] = SpaceCorps(cave(currfishnum), 1, 1:length(cave(currfishnum).fish), [(j-1)*epochdur, j*epochdur]);
        [~, cvQ] = overture(tmp);
        caveReal = caveReal + cvQ.realoverlaps;
        caveRand = caveRand + cvQ.randoverlaps;
        
    end

%% cave 8, 10 fish, some gaps over 1300 seconds, lovely data

currfishnum = 8; maxdur = 1300;

    for j=1:floor(maxdur/epochdur)
       
        [tmp, ~] = SpaceCorps(cave(currfishnum), 1, 1:length(cave(currfishnum).fish), [(j-1)*epochdur, j*epochdur]);
        [~, cvQ] = overture(tmp);
        caveReal = caveReal + cvQ.realoverlaps;
        caveRand = caveRand + cvQ.randoverlaps;
        
    end


%% cave 9, 1 fish, over 1000 seconds

%% cave 10, 3 fish, over 900 seconds

currfishnum = 10; maxdur = 900;

    for j=1:floor(maxdur/epochdur)
       
        [tmp, ~] = SpaceCorps(cave(currfishnum), 1, 1:length(cave(currfishnum).fish), [(j-1)*epochdur, j*epochdur]);
        [~, cvQ] = overture(tmp);
        caveReal = caveReal + cvQ.realoverlaps;
        caveRand = caveRand + cvQ.randoverlaps;
        
    end

%% cave 11, 7 fish, over 900 seconds

currfishnum = 11; maxdur = 900;

    for j=1:floor(maxdur/epochdur)
       
        [tmp, ~] = SpaceCorps(cave(currfishnum), 1, 1:length(cave(currfishnum).fish), [(j-1)*epochdur, j*epochdur]);
        [~, cvQ] = overture(tmp);
        caveReal = caveReal + cvQ.realoverlaps;
        caveRand = caveRand + cvQ.randoverlaps;
        
    end

%% cave 12, 7 fish, over 1100 seconds

currfishnum = 12; maxdur = 1100;

    for j=1:floor(maxdur/epochdur)
       
        [tmp, ~] = SpaceCorps(cave(currfishnum), 1, 1:length(cave(currfishnum).fish), [(j-1)*epochdur, j*epochdur]);
        [~, cvQ] = overture(tmp);
        caveReal = caveReal + cvQ.realoverlaps;
        caveRand = caveRand + cvQ.randoverlaps;
        
    end

%% cave 13, 9 fish, over 1200 seconds, some complexity, some in and out

currfishnum = 13; maxdur = 1200;

    for j=1:floor(maxdur/epochdur)
       
        [tmp, ~] = SpaceCorps(cave(currfishnum), 1, 1:length(cave(currfishnum).fish), [(j-1)*epochdur, j*epochdur]);
        [~, cvQ] = overture(tmp);
        caveReal = caveReal + cvQ.realoverlaps;
        caveRand = caveRand + cvQ.randoverlaps;
        
    end


%% cave 14, 7 fish, over 1000 seconds

curlen = length(fish);

tmp = MajorTom(cave(14));
    for j=1:length(tmp); fish(j+curlen) = tmp(j); end
    
    [spadat(10).dat, ~] = SpaceCorps(cave(14), 1);
    [sa(10).sa, sq(10).sq] = overture(spadat(10).dat);


    
CaveIDX = length(fish)
CaveANA = 10;

%% surface 1, 12 fish, over 600 seconds

% TUBE FISH, index 7, 9, 11

curlen = length(fish);

tmp = MajorTom(srf(1), [0 srf(1).fish(1).freq(end,1)], [1 2 3 4 5 6 8 10 12]);
    for j=1:length(tmp); fish(j+curlen) = tmp(j); end

tmp = MajorTom(srf(1), [0 srf(1).fish(1).freq(end,1)], [7 9 11]);
    for j=1:length(tmp); tubular(j) = tmp(j); end    
    
    [spadat(11).dat, ~] = SpaceCorps(srf(1), 2, [1 2 3 4 5 6 8 10 12]);
    [sa(11).sa, sq(11).sq] = overture(spadat(11).dat);


%% surface 2, 21 fish, over 1000 seconds

% TUBE FISH, index 11, 18, 19

curlen = length(fish);
tmp = MajorTom(srf(2), [0 srf(2).fish(1).freq(end,1)], [1 2 3 4 5 6 7 8 10 12 13 14 15 16 17 20 21]);
    for j=1:length(tmp); fish(j+curlen) = tmp(j); end

curlen = length(tubular);
tmp = MajorTom(srf(2), [0 srf(2).fish(1).freq(end,1)], [11 18 19]);
    for j=1:length(tmp); tubular(j+curlen) = tmp(j); end    
    
    
     [spadat(12).dat, ~] = SpaceCorps(srf(2), 2, [1 2 3 4 5 6 7 8 10 12 13 14 15 16 17 20 21]);
     [sa(12).sa, sq(12).sq] = overture(spadat(12).dat);


%% surface 3, 31 fish, over 600 seconds

% TUBE FISH, index 22, 28, MISHMASH @ 385 Hz

curlen = length(fish);
tmp = MajorTom(srf(3), [0 srf(3).fish(1).freq(end,1)], [1 2 3 4 5 6 7 8 10 11 12 13 14 15 16 17 18 19 20 21 23 24 25 26 27 29 30 31]);
    for j=1:length(tmp); fish(j+curlen) = tmp(j); end

curlen = length(tubular);
tmp = MajorTom(srf(3), [0 srf(3).fish(1).freq(end,1)], [22 28]);
    for j=1:length(tmp); tubular(j+curlen) = tmp(j); end    
    
    
    [spadat(13).dat, ~] = SpaceCorps(srf(3), 2, [1 2 3 4 5 6 7 8 10 11 12 13 14 15 16 17 18 19 20 21 23 24 25 26 27 29 30 31]);
    [sa(13).sa, sq(13).sq] = overture(spadat(13).dat);


%% surface 4, 22 fish, over 600 seconds [1 2 3 4 5 6 7 8 10 12 13 14 15 17 18 19 21 22]

% TUBE FISH, index 11, 16, 20
curlen = length(fish);
tmp = MajorTom(srf(4), [0 srf(4).fish(1).freq(end,1)], [1 2 3 4 5 6 7 8 10 12 13 14 15 17 18 19 21 22]);
    for j=1:length(tmp); fish(j+curlen) = tmp(j); end

curlen = length(tubular);
tmp = MajorTom(srf(4), [0 srf(4).fish(1).freq(end,1)], [11 16 20]);
    for j=1:length(tmp); tubular(j+curlen) = tmp(j); end    
    
    
    [spadat(14).dat, ~] = SpaceCorps(srf(4), 2, [1 2 3 4 5 6 7 8 10 12 13 14 15 17 18 19 21 22]);
    [sa(14).sa, sq(14).sq] = overture(spadat(14).dat);


%% surface 5, 24 fish, over 1100 seconds
curlen = length(fish);
tmp = MajorTom(srf(5), [0 srf(5).fish(1).freq(end,1)]);
    for j=1:length(tmp); fish(j+curlen) = tmp(j); end

    [spadat(15).dat, ~] = SpaceCorps(srf(5), 2);
    [sa(15).sa, sq(15).sq] = overture(spadat(15).dat);



%% Stats for single fish

% figure(1); clf; 
%     subplot(121); plot([fish.meanfreq], '-*'); ylim([200 600]);
%     subplot(122); plot([fish.varfreq], '-*'); ylim([0 2]);

CaveRealOver = []; SurfRealOver = [];
CaveRandOver = []; SurfRandOver = [];

for j=1:CaveANA   
    CaveRealOver = [CaveRealOver sq(j).sq.realoverlaps];
    CaveRandOver = [CaveRandOver sq(j).sq.randoverlaps];    
end
for j=CaveANA+1:length(sq)
    SurfRealOver = [SurfRealOver sq(j).sq.realoverlaps];
    SurfRandOver = [SurfRandOver sq(j).sq.randoverlaps];        
end

lms = 0:0.05:1;
figure(27); clf;
    subplot(221); histogram(CaveRealOver, 'BinEdges', lms);
    subplot(222); histogram(CaveRandOver, 'BinEdges', lms);
    subplot(223); histogram(SurfRealOver, 'BinEdges', lms);
    subplot(224); histogram(SurfRandOver, 'BinEdges', lms);


