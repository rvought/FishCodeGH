function [srfVsrf, cvVcv, cvVsrf, data] = Miranda(cave, srf)



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
        caveReal = [caveReal, cvQ.realoverlaps];
        caveRand = [caveRand, cvQ.randoverlaps];
        
    end



%% cave 4, 7 fish, good recordings of all fish through over 1200 seconds

currfishnum = 4; maxdur = 1200;

    for j=1:floor(maxdur/epochdur)
       
        [tmp, ~] = SpaceCorps(cave(currfishnum), 1, 1:length(cave(currfishnum).fish), [(j-1)*epochdur, j*epochdur]);
        [~, cvQ] = overture(tmp);
        caveReal = [caveReal, cvQ.realoverlaps];
        caveRand = [caveRand, cvQ.randoverlaps];
        
    end

%% cave 5, 7 fish, good recordings of all fish through over 1200 seconds

currfishnum = 5; maxdur = 1200;

    for j=1:floor(maxdur/epochdur)
       
        [tmp, ~] = SpaceCorps(cave(currfishnum), 1, 1:length(cave(currfishnum).fish), [(j-1)*epochdur, j*epochdur]);
        [~, cvQ] = overture(tmp);
        caveReal = [caveReal, cvQ.realoverlaps];
        caveRand = [caveRand, cvQ.randoverlaps];
        
    end

%% cave 6, 8 fish, good recordings of all fish through over 1200 seconds LOTS OF SOCIAL

currfishnum = 6; maxdur = 1200;

    for j=1:floor(maxdur/epochdur)
       
        [tmp, ~] = SpaceCorps(cave(currfishnum), 1, 1:length(cave(currfishnum).fish), [(j-1)*epochdur, j*epochdur]);
        [~, cvQ] = overture(tmp);
            caveReal = [caveReal, cvQ.realoverlaps];
            caveRand = [caveRand, cvQ.randoverlaps];        
    end


%% cave 7, 7 fish, some gaps over 900 seconds

currfishnum = 7; maxdur = 900;

    for j=1:floor(maxdur/epochdur)
       
        [tmp, ~] = SpaceCorps(cave(currfishnum), 1, 1:length(cave(currfishnum).fish), [(j-1)*epochdur, j*epochdur]);
        [~, cvQ] = overture(tmp);
            caveReal = [caveReal, cvQ.realoverlaps];
            caveRand = [caveRand, cvQ.randoverlaps];
        
    end

%% cave 8, 10 fish, some gaps over 1300 seconds, lovely data

currfishnum = 8; maxdur = 1300;

    for j=1:floor(maxdur/epochdur)
       
        [tmp, ~] = SpaceCorps(cave(currfishnum), 1, 1:length(cave(currfishnum).fish), [(j-1)*epochdur, j*epochdur]);
        [~, cvQ] = overture(tmp);
            caveReal = [caveReal, cvQ.realoverlaps];
            caveRand = [caveRand, cvQ.randoverlaps];
        
    end


%% cave 9, 1 fish, over 1000 seconds

%% cave 10, 3 fish, over 900 seconds

currfishnum = 10; maxdur = 900;

    for j=1:floor(maxdur/epochdur)
       
        [tmp, ~] = SpaceCorps(cave(currfishnum), 1, 1:length(cave(currfishnum).fish), [(j-1)*epochdur, j*epochdur]);
        [~, cvQ] = overture(tmp);
            caveReal = [caveReal, cvQ.realoverlaps];
            caveRand = [caveRand, cvQ.randoverlaps];
        
    end

%% cave 11, 7 fish, over 900 seconds

currfishnum = 11; maxdur = 900;

    for j=1:floor(maxdur/epochdur)
       
        [tmp, ~] = SpaceCorps(cave(currfishnum), 1, 1:length(cave(currfishnum).fish), [(j-1)*epochdur, j*epochdur]);
        [~, cvQ] = overture(tmp);
            caveReal = [caveReal, cvQ.realoverlaps];
            caveRand = [caveRand, cvQ.randoverlaps];
        
    end

%% cave 12, 7 fish, over 1100 seconds

currfishnum = 12; maxdur = 1100;

    for j=1:floor(maxdur/epochdur)
       
        [tmp, ~] = SpaceCorps(cave(currfishnum), 1, 1:length(cave(currfishnum).fish), [(j-1)*epochdur, j*epochdur]);
        [~, cvQ] = overture(tmp);
            caveReal = [caveReal, cvQ.realoverlaps];
            caveRand = [caveRand, cvQ.randoverlaps];
        
    end

%% cave 13, 9 fish, over 1200 seconds, some complexity, some in and out

currfishnum = 13; maxdur = 1200;

    for j=1:floor(maxdur/epochdur)
       
        [tmp, ~] = SpaceCorps(cave(currfishnum), 1, 1:length(cave(currfishnum).fish), [(j-1)*epochdur, j*epochdur]);
        [~, cvQ] = overture(tmp);
            caveReal = [caveReal, cvQ.realoverlaps];
            caveRand = [caveRand, cvQ.randoverlaps];
        
    end


%% cave 14, 7 fish, over 1000 seconds

currfishnum = 14; maxdur = 1000;

    for j=1:floor(maxdur/epochdur)
       
        [tmp, ~] = SpaceCorps(cave(currfishnum), 1, 1:length(cave(currfishnum).fish), [(j-1)*epochdur, j*epochdur]);
        [~, cvQ] = overture(tmp);
            caveReal = [caveReal, cvQ.realoverlaps];
            caveRand = [caveRand, cvQ.randoverlaps];
        
    end

%% surface 1, 12 fish, over 600 seconds

% TUBE FISH, index 7, 9, 11
    
currfishnum = 1; maxdur = 600; goodfish = [1 2 3 4 5 6 8 10 12];

    [tmp, ~] = SpaceCorps(srf(currfishnum), 2, goodfish, [0 epochdur]);
    [~, sfQ] = overture(tmp);
    surfReal = sfQ.realoverlaps;
    surfRand = sfQ.randoverlaps;
    
    for j=2:floor(maxdur/epochdur)
       
        [tmp, ~] = SpaceCorps(srf(currfishnum), 2, goodfish, [(j-1)*epochdur, j*epochdur]);
        [~, sfQ] = overture(tmp);
            surfReal = [surfReal, sfQ.realoverlaps];
            surfRand = [surfRand, sfQ.randoverlaps];
        
    end
        

%% surface 2, 21 fish, over 1000 seconds

% TUBE FISH, index 11, 18, 19

currfishnum = 2; maxdur = 1000; goodfish = [1 2 3 4 5 6 7 8 10 12 13 14 15 16 17 20 21];

    for j=1:floor(maxdur/epochdur)
       
        [tmp, ~] = SpaceCorps(srf(currfishnum), 2, goodfish, [(j-1)*epochdur, j*epochdur]);
        [~, sfQ] = overture(tmp);
            surfReal = [surfReal, sfQ.realoverlaps];
            surfRand = [surfRand, sfQ.randoverlaps];
        
    end


%% surface 3, 31 fish, over 600 seconds

% TUBE FISH, index 22, 28, MISHMASH @ 385 Hz

currfishnum = 3; maxdur = 1000; goodfish = [1 2 3 4 5 6 7 8 10 11 12 13 14 15 16 17 18 19 20 21 23 24 25 26 27 29 30 31];

    for j=1:floor(maxdur/epochdur)
       
        [tmp, ~] = SpaceCorps(srf(currfishnum), 2, goodfish, [(j-1)*epochdur, j*epochdur]);
        [~, sfQ] = overture(tmp);
            surfReal = [surfReal, sfQ.realoverlaps];
            surfRand = [surfRand, sfQ.randoverlaps];
        
    end

%% surface 4, 22 fish, over 600 seconds [1 2 3 4 5 6 7 8 10 12 13 14 15 17 18 19 21 22]

% TUBE FISH, index 11, 16, 20

currfishnum = 4; maxdur = 600; goodfish = [1 2 3 4 5 6 7 8 10 12 13 14 15 17 18 19 21 22];

    for j=1:floor(maxdur/epochdur)
       
        [tmp, ~] = SpaceCorps(srf(currfishnum), 2, goodfish, [(j-1)*epochdur, j*epochdur]);
        [~, sfQ] = overture(tmp);
            surfReal = [surfReal, sfQ.realoverlaps];
            surfRand = [surfRand, sfQ.randoverlaps];
        
    end
   

%% surface 5, 24 fish, over 1100 seconds

currfishnum = 5; maxdur = 1100; 

    for j=1:floor(maxdur/epochdur)
       
        [tmp, ~] = SpaceCorps(srf(currfishnum), 2, 1:length(srf(currfishnum).fish), [(j-1)*epochdur, j*epochdur]);
        [~, sfQ] = overture(tmp);
            surfReal = [surfReal, sfQ.realoverlaps];
            surfRand = [surfRand, sfQ.randoverlaps];
        
    end

    data.surfReal = surfReal;
    data.surfRand = surfRand;
    data.caveRand = caveRand;
    data.caveReal = caveReal;
    
[srfVsrf.H, srfVsrf.P, srfVsrf.K] = kstest2(histcounts(surfReal,20), histcounts(surfRand,20))
[cvVcv.H, cvVcv.P, cvVcv.K] = kstest2(histcounts(caveReal,20), histcounts(caveRand,20))
[cvVsrf.H, cvVsrf.P, cvVsrf.K] = kstest2(histcounts(caveReal,20), histcounts(caveRand,20))


