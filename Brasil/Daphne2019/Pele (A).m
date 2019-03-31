function [fish, w] = Pele(cave, srf)


CutJumpDist = 10; % We have bad XY coordinates, especially when fish bump into the
% electrodes.  If the fish jumps more that CutJumpDist cm between two samples, 
% then we consider that a bad measurement. Sample rate is between 4 and 5
% Hz, giving us a speed of 40-50cm/second.

solocave = [];
groupcave = [];
allsurf = [];
tubefish = [];

%% cave 1, 1 fish, 1000 seconds

% Let's do two sections of 500 seconds

    fish(1) = MajorTom(cave(1), [0 500], CutJumpDist); solocave = [solocave, length(fish)];
    fish(end+1) = MajorTom(cave(1), [500 1000], CutJumpDist);  solocave = [solocave, length(fish)];   
    
%% cave 2, 1 fish, 900 seconds

    fish(end+1) = MajorTom(cave(2), [0 450], CutJumpDist); solocave = [solocave, length(fish)];
    fish(end+1) = MajorTom(cave(2), [450 900], CutJumpDist); solocave = [solocave, length(fish)];

%% cave 3, 7 fish, good recordings of all fish through over 1200 seconds

    tmp = MajorTom(cave(3), [0 600], CutJumpDist);
        for z=1:length(tmp); fish(end+1) = tmp(z); groupcave = [groupcave, length(fish)]; end  
    tmp = MajorTom(cave(3), [600 1200], CutJumpDist);
        for z=1:length(tmp); fish(end+1) = tmp(z); groupcave = [groupcave, length(fish)]; end

        length(fish)

%     [tmpOUT, tmpCMBS] = SpaceCorps(cave(1), 1, [0 600]);
%     [sa, sq] = overture(tmpOUT);
% 
%     [tmpOUT, tmpCMBS] = SpaceCorps(cave(1), 1, [600 1200]);
%     [sa, sq] = overture(tmpOUT);



%% cave 4, 7 fish, good recordings of all fish through over 1200 seconds

    tmp = MajorTom(cave(4), [0 600], CutJumpDist);
        for z=1:length(tmp); fish(end+1) = tmp(z); groupcave = [groupcave, length(fish)]; end  
    tmp = MajorTom(cave(4), [600 1200], CutJumpDist);
        for z=1:length(tmp); fish(end+1) = tmp(z); groupcave = [groupcave, length(fish)]; end


%% cave 5, 7 fish, good recordings of all fish through over 1200 seconds

    tmp = MajorTom(cave(5), [0 600], CutJumpDist);
        for z=1:length(tmp); fish(end+1) = tmp(z); groupcave = [groupcave, length(fish)]; end  
    tmp = MajorTom(cave(5), [600 1200], CutJumpDist);
        for z=1:length(tmp); fish(end+1) = tmp(z); groupcave = [groupcave, length(fish)]; end


%% cave 6, 8 fish, good recordings of all fish through over 1200 seconds LOTS OF SOCIAL

    tmp = MajorTom(cave(6), [0 600], CutJumpDist);
        for z=1:length(tmp); fish(end+1) = tmp(z); groupcave = [groupcave, length(fish)]; end  
    tmp = MajorTom(cave(6), [600 1200], CutJumpDist);
        for z=1:length(tmp); fish(end+1) = tmp(z); groupcave = [groupcave, length(fish)]; end


%% cave 7, 7 fish, some gaps over 900 seconds


%% cave 8, 10 fish, some gaps over 1300 seconds, lovely data

    tmp = MajorTom(cave(8), [0 600], CutJumpDist);
        for z=1:length(tmp); fish(end+1) = tmp(z); groupcave = [groupcave, length(fish)]; end  
    tmp = MajorTom(cave(8), [600 1200], CutJumpDist);
        for z=1:length(tmp); fish(end+1) = tmp(z); groupcave = [groupcave, length(fish)]; end


%% cave 9, 1 fish, over 1000 seconds

    fish(end+1) = MajorTom(cave(9), [0 500], CutJumpDist); solocave = [solocave, length(fish)];
    fish(end+1) = MajorTom(cave(9), [500 1000], CutJumpDist); solocave = [solocave, length(fish)];

    length(fish)

%% cave 10, 3 fish, over 900 seconds


%% cave 11, 7 fish, over 900 seconds


%% cave 12, 7 fish, over 1100 seconds

    tmp = MajorTom(cave(12), [0 500], CutJumpDist);
        for z=1:length(tmp); fish(end+1) = tmp(z); groupcave = [groupcave, length(fish)]; end
    tmp = MajorTom(cave(12), [500 1000], CutJumpDist);
        for z=1:length(tmp); fish(end+1) = tmp(z); groupcave = [groupcave, length(fish)]; end

length(fish)
%% cave 13, 9 fish, over 1200 seconds, some complexity, some in and out




%% cave 14, 7 fish, over 1000 seconds



%% surface 1, 12 fish, over 600 seconds

% TUBE FISH, index 7, 9, 11

    tmp = MajorTom(srf(1), [0 300], CutJumpDist, [1:6 8 10 12]);
        for z=1:length(tmp); fish(end+1) = tmp(z); allsurf = [allsurf, length(fish)]; end
    tmp = MajorTom(srf(1), [0 300], CutJumpDist, [7 9 11]);
        for z=1:length(tmp); fish(end+1) = tmp(z); tubefish = [tubefish, length(fish)]; end
    tmp = MajorTom(srf(1), [300 600], CutJumpDist, [1:6 8 10 12]);
        for z=1:length(tmp); fish(end+1) = tmp(z); allsurf = [allsurf, length(fish)]; end
    tmp = MajorTom(srf(1), [300 600], CutJumpDist, [7 9 11]);
        for z=1:length(tmp); fish(end+1) = tmp(z); tubefish = [tubefish, length(fish)]; end



%% surface 2, 21 fish, over 1000 seconds

% TUBE FISH, index 11, 18, 19

    tmp = MajorTom(srf(2), [0 500], CutJumpDist, [1:10 12:17 20:21]);
        for z=1:length(tmp); fish(end+1) = tmp(z); allsurf = [allsurf, length(fish)]; end
%     tmp = MajorTom(srf(2), [0 500], CutJumpDist, [11 18 19]);
%         for z=1:length(tmp); fish(end+1) = tmp(z); tubefish = [tubefish, length(fish)]; end
    tmp = MajorTom(srf(2), [500 1000], CutJumpDist, [1:10 12:17 20:21]);
        for z=1:length(tmp); fish(end+1) = tmp(z); allsurf = [allsurf, length(fish)]; end
%     tmp = MajorTom(srf(2), [500 1000], CutJumpDist, [11 18 19]);
%         for z=1:length(tmp); fish(end+1) = tmp(z); tubefish = [tubefish, length(fish)]; end



%% surface 3, 31 fish, over 600 seconds

% TUBE FISH, index 22, 28, MISHMASH @ 385 Hz


%% surface 4, 22 fish, over 600 seconds

% TUBE FISH, index 11, 16, 20


%% surface 5, 24 fish, over 1100 seconds

    tmp = MajorTom(srf(5), [0 500], CutJumpDist);
        for z=1:length(tmp); fish(end+1) = tmp(z); allsurf = [allsurf, length(fish)]; end
    tmp = MajorTom(srf(5), [500 1000], CutJumpDist);
        for z=1:length(tmp); fish(end+1) = tmp(z); allsurf = [allsurf, length(fish)]; end
 length(fish)

%% Stats for single fish

figure(1); clf; hold on;
    plot([fish(groupcave).meanfreq], [fish(groupcave).meanvel], 'r*', 'MarkerSize', 6);
	plot([fish(solocave).meanfreq], [fish(solocave).meanvel], 'k*', 'MarkerSize', 6);
figure(2); clf; hold on;
	plot([fish(allsurf).meanfreq], [fish(allsurf).meanvel], 'r*', 'MarkerSize', 6);
	plot([fish(tubefish).meanfreq], [fish(tubefish).meanvel], 'k*', 'MarkerSize', 6);

w.tubefish = tubefish;
w.allsurf = allsurf;
w.solocave = solocave;
w.groupcave = groupcave;


