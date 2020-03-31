
   
StepSize = 300;

%% Immobilized Surface Fish in the grid

Sstds = [];

% Surface recording #1
for fishidx = [7 9 11] 
    for j = 1:2 % 0 to 600
        tt = srf(1).fish(fishidx).freq(:,1) > StepSize*(j-1) & srf(1).fish(fishidx).freq(:,1) < StepSize*j;
        Sstds(end+1) = nanstd(srf(1).fish(fishidx).freq(tt,2));
    end
end
% Surface recording #2
for fishidx = [11 18 19] 
    for j = 1:3 % 0 to 900
        tt = srf(2).fish(fishidx).freq(:,1) > StepSize*(j-1) & srf(2).fish(fishidx).freq(:,1) < StepSize*j;
        Sstds(end+1) = nanstd(srf(2).fish(fishidx).freq(tt,2));
    end
end
% Surface recording #3
for fishidx = [22 28] 
    for j = 1:2 % 0 to 600
        tt = srf(3).fish(fishidx).freq(:,1) > StepSize*(j-1) & srf(3).fish(fishidx).freq(:,1) < StepSize*j;
        Sstds(end+1) = nanstd(srf(3).fish(fishidx).freq(tt,2));
    end
end
% Surface recording #4
for fishidx = [11 16 20] 
    for j = 1:2 % 0 to 600
        tt = srf(4).fish(fishidx).freq(:,1) > StepSize*(j-1) & srf(4).fish(fishidx).freq(:,1) < StepSize*j;
        Sstds(end+1) = nanstd(srf(4).fish(fishidx).freq(tt,2));
    end
end
% Surface recording #5
for fishidx = [3 23] 
    for j = 1:3 % 0 to 900
        tt = srf(5).fish(fishidx).freq(:,1) > StepSize*(j-1) & srf(5).fish(fishidx).freq(:,1) < StepSize*j;
        Sstds(end+1) = nanstd(srf(5).fish(fishidx).freq(tt,2));
    end
end

Sstds = Sstds(~isnan(Sstds));
Sstds = Sstds(Sstds ~= 0);

fprintf('Mean Surface Immobile Var = %2.8f, Var %2.8f', mean(Sstds), var(Sstds));

%% Moving Surface Fish in the grid

Gstds = [];

idx = 1;
for fishidx = [1 2 3 4 5 6  8  10  12]
    for j = 1:2 % 0 to 600
        tt = srf(1).fish(fishidx).freq(:,1) > StepSize*(j-1) & srf(1).fish(fishidx).freq(:,1) < StepSize*j;
        Gstds(end+1) = nanstd(srf(1).fish(fishidx).freq(tt,2));
    end
end

% Surface recording #2
for fishidx = [1 2 3 4 5 6 7 8 9 10  12 13 14 15 16 17   20 21]
    for j = 1:3 % 0 to 900
        tt = srf(2).fish(fishidx).freq(:,1) > StepSize*(j-1) & srf(2).fish(fishidx).freq(:,1) < StepSize*j;
        Gstds(end+1) = nanstd(srf(2).fish(fishidx).freq(tt,2));
    end
end

% Surface recording #3
for fishidx = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21  23 24 25 26 27  29 30 31]
    for j = 1:2 % 0 to 600
        tt = srf(3).fish(fishidx).freq(:,1) > StepSize*(j-1) & srf(3).fish(fishidx).freq(:,1) < StepSize*j;
        Gstds(end+1) = nanstd(srf(3).fish(fishidx).freq(tt,2));
    end
end

% Surface recording #4
for fishidx = [1 2 3 4 5 6 7 8 9 10  12 13 14 15  17 18 19  21 22]
    for j = 1:2 % 0 to 600
        tt = srf(4).fish(fishidx).freq(:,1) > StepSize*(j-1) & srf(4).fish(fishidx).freq(:,1) < StepSize*j;
        Gstds(end+1) = nanstd(srf(4).fish(fishidx).freq(tt,2));
    end
end

% Surface recording #5
for fishidx = [1 2  4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22  24]
    for j = 1:3 % 0 to 900
        tt = srf(5).fish(fishidx).freq(:,1) > StepSize*(j-1) & srf(5).fish(fishidx).freq(:,1) < StepSize*j;
        Gstds(end+1) = nanstd(srf(5).fish(fishidx).freq(tt,2));
    end
end

Gstds = Gstds(~isnan(Gstds));
Gstds = Gstds(Gstds ~= 0);

%% Freely moving solitary fish in the cave

SCstds = [];

% Cave recording #1 - 1 fish
    for j = 1:3 % 0 to 900
        tt = cave(1).fish(1).freq(:,1) > StepSize*(j-1) & cave(1).fish(1).freq(:,1) < StepSize*j;
        SCstds(end+1) = nanstd(cave(1).fish(1).freq(tt,2));
    end
% Cave recording #2 - 1 fish
    for j = 1:3 % 0 to 900
        tt = cave(2).fish(1).freq(:,1) > StepSize*(j-1) & cave(2).fish(1).freq(:,1) < StepSize*j;
        SCstds(end+1) = nanstd(cave(2).fish(1).freq(tt,2));
    end
% Cave recording #9 - 1 fish
    for j = 1:3 % 0 to 900
        tt = cave(9).fish(1).freq(:,1) > StepSize*(j-1) & cave(9).fish(1).freq(:,1) < StepSize*j;
        SCstds(end+1) = nanstd(cave(9).fish(1).freq(tt,2));
    end

SCstds = SCstds(~isnan(SCstds));
SCstds = SCstds(SCstds ~=0);



