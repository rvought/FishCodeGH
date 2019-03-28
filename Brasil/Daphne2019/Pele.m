function [fish, sa, sq, spadat] = Pele(cave, srf)



%% cave 1, 1 fish, 1000 seconds

fish(1) = MajorTom(cave(1));

% Let's do two sections of 500 seconds

%     fish(1) = MajorTom(cave(1), [0 500]);
%     fish(2) = MajorTom(cave(1), [500 1000]);

%% cave 2, 1 fish, 900 seconds

fish(2) = MajorTom(cave(2));

%     fish(3) = MajorTom(cave(2), [0 450]);
%     fish(4) = MajorTom(cave(2), [450 900]);

%% cave 3, 7 fish, good recordings of all fish through over 1200 seconds

curlen = length(fish);

tmp = MajorTom(cave(3));
    for j=1:length(tmp); fish(j+curlen) = tmp(j); end

%     fish(5) = MajorTom(cave(3), [0 600]);
%     fish(6) = MajorTom(cave(3), [600 1200]);

    [spadat(1).dat, ~] = SpaceCorps(cave(3), 1);
    [sa(1).sa, sq(1).sq] = overture(spadat(1).dat);

%     [tmpOUT, tmpCMBS] = SpaceCorps(cave(1), 1, [600 1200]);
%     [sa, sq] = overture(tmpOUT);



%% cave 4, 7 fish, good recordings of all fish through over 1200 seconds

curlen = length(fish);

tmp = MajorTom(cave(4));
    for j=1:length(tmp); fish(j+curlen) = tmp(j); end
    
    [spadat(2).dat, ~] = SpaceCorps(cave(4), 1);
    [sa(2).sa, sq(2).sq] = overture(spadat(2).dat);


%% cave 5, 7 fish, good recordings of all fish through over 1200 seconds
curlen = length(fish);

tmp = MajorTom(cave(5));
    for j=1:length(tmp); fish(j+curlen) = tmp(j); end
    
    [spadat(3).dat, ~] = SpaceCorps(cave(5), 1);
    [sa(3).sa, sq(3).sq] = overture(spadat(3).dat);


%% cave 6, 8 fish, good recordings of all fish through over 1200 seconds LOTS OF SOCIAL

curlen = length(fish);

tmp = MajorTom(cave(6));
    for j=1:length(tmp); fish(j+curlen) = tmp(j); end
    
    [spadat(4).dat, ~] = SpaceCorps(cave(6), 1);
    [sa(4).sa, sq(4).sq] = overture(spadat(4).dat);



%% cave 7, 7 fish, some gaps over 900 seconds

curlen = length(fish);

tmp = MajorTom(cave(7));
    for j=1:length(tmp); fish(j+curlen) = tmp(j); end
    
    [spadat(5).dat, ~] = SpaceCorps(cave(7), 1);
    [sa(5).sa, sq(5).sq] = overture(spadat(5).dat);


%% cave 8, 10 fish, some gaps over 1300 seconds, lovely data


%% cave 9, 1 fish, over 1000 seconds

% curlen = length(fish);
% 
% tmp = MajorTom(cave(9));
%     for j=1:length(tmp); fish(j+curlen) = tmp(j); end
%     
%     [spadat(6).dat, ~] = SpaceCorps(cave(9), 1);
%     [sa(6).sa, sq(6).sq] = overture(spadat(6).dat);


%% cave 10, 3 fish, over 900 seconds


%% cave 11, 7 fish, over 900 seconds


%% cave 12, 7 fish, over 1100 seconds



%% cave 13, 9 fish, over 1200 seconds, some complexity, some in and out



%% cave 14, 7 fish, over 1000 seconds



%% surface 1, 12 fish, over 600 seconds

% TUBE FISH, index 7, 9, 11

curlen = length(fish)

tmp = MajorTom(srf(1), [0 srf(1).fish(1).freq(end,1)], [1 2 3 4 5 6 8 10 12]);
    for j=1:length(tmp); fish(j+curlen) = tmp(j); end

    
    [spadat(6).dat, ~] = SpaceCorps(srf(1), 1, [1 2 3 4 5 6 8 10 12]);
    [sa(6).sa, sq(6).sq] = overture(spadat(6).dat);


%% surface 2, 21 fish, over 1000 seconds

% TUBE FISH, index 11, 18, 19



%% surface 3, 31 fish, over 600 seconds

% TUBE FISH, index 22, 28, MISHMASH @ 385 Hz


%% surface 4, 22 fish, over 600 seconds

% TUBE FISH, index 11, 16, 20
curlen = length(fish);
tmp = MajorTom(srf(4), [0 srf(4).fish(1).freq(end,1)], [1 2 3 4 5 6 7 8 10 12 13 14 15 17 18 19 21 22]);
    for j=1:length(tmp); fish(j+curlen) = tmp(j); end
    
    [spadat(7).dat, ~] = SpaceCorps(srf(4), 2, [1 2 3 4 5 6 7 8 10 12 13 14 15 17 18 19 21 22]);
    [sa(7).sa, sq(7).sq] = overture(spadat(7).dat);


%% surface 5, 24 fish, over 1100 seconds



%% Stats for single fish

% figure(1); clf; 
%     subplot(121); plot([fish.meanfreq], '-*'); ylim([200 600]);
%     subplot(122); plot([fish.varfreq], '-*'); ylim([0 2]);



