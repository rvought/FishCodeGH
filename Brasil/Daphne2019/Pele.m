function [fish, pairs, stts] = Pele(cave, srf)



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

tmp = MajorTom(cave(3));
    for j=1:length(tmp); fish(end+1) = tmp(j); end;

%     fish(5) = MajorTom(cave(3), [0 600]);
%     fish(6) = MajorTom(cave(3), [600 1200]);

    [spadat(1), cmbs(1)] = SpaceCorps(cave(1), 1);
    [sa(1), sq(1)] = overture(spadat(1));

%     [tmpOUT, tmpCMBS] = SpaceCorps(cave(1), 1, [600 1200]);
%     [sa, sq] = overture(tmpOUT);



%% cave 4, 7 fish, good recordings of all fish through over 1200 seconds

fish(4) = MajorTom(cave(4));

    [spadat(2), cmbs(2)] = SpaceCorps(cave(4), 1);
    [sa(2), sq(2)] = overture(spadat(2));


%% cave 5, 7 fish, good recordings of all fish through over 1200 seconds


%% cave 6, 8 fish, good recordings of all fish through over 1200 seconds LOTS OF SOCIAL


%% cave 7, 7 fish, some gaps over 900 seconds


%% cave 8, 10 fish, some gaps over 1300 seconds, lovely data


%% cave 9, 1 fish, over 1000 seconds

%     fish(7) = MajorTom(cave(9), [0 500]);
%     fish(8c) = MajorTom(cave(9), [500 1000]);


%% cave 10, 3 fish, over 900 seconds


%% cave 11, 7 fish, over 900 seconds


%% cave 12, 7 fish, over 1100 seconds



%% cave 13, 9 fish, over 1200 seconds, some complexity, some in and out



%% cave 14, 7 fish, over 1000 seconds



%% surface 1, 12 fish, over 600 seconds

% TUBE FISH, index 7, 9, 11
fish(5) = MajorTom(srf(1), [0 srf(1).fish(1).freq(end,1)], [1 2 3 4 5 6 8 10 12]);

    [spadat(3), cmbs(3)] = SpaceCorps(srf(1), 2);
    [sa(3), sq(3)] = overture(spadat(3));


%% surface 2, 21 fish, over 1000 seconds

% TUBE FISH, index 11, 18, 19



%% surface 3, 31 fish, over 600 seconds

% TUBE FISH, index 22, 28, MISHMASH @ 385 Hz


%% surface 4, 22 fish, over 600 seconds

% TUBE FISH, index 11, 16, 20
fish(6) = MajorTom(srf(4), [0 srf(4).fish(1).freq(end,1)], [1 2 3 4 5 6 8 10 12]);

    [spadat(4), cmbs(4)] = SpaceCorps(srf(4), 2);
    [sa(4), sq(4)] = overture(spadat(4));


%% surface 5, 24 fish, over 1100 seconds



%% Stats for single fish

% figure(1); clf; 
%     subplot(121); plot([fish.meanfreq], '-*'); ylim([200 600]);
%     subplot(122); plot([fish.varfreq], '-*'); ylim([0 2]);



