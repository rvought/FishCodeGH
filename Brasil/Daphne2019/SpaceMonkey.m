

StepSize = 300;
CfishArea = []; 
Coverlaps = [];
CpairAreas = [];

% [ strt, stp ] = brasilsamplelist(id, fish);

%% Cave analysis

idx = 3;
for j = 1:4 % 0 to 1200
    out = SpaceCorps(cave(idx), 1:length(cave(idx).fish), [StepSize*(j-1), StepSize*j]);
    foo(j) = SpaceAnal(out);
end

idx = 4;
for j = 1:4 % 0 to 1200
    out = SpaceCorps(cave(idx), 1:length(cave(idx).fish), [StepSize*(j-1), StepSize*j]);
    foo(end+1) = SpaceAnal(out);
end

idx = 5;
for j = 1:4 % 0 to 1200
    out = SpaceCorps(cave(idx), 1:length(cave(idx).fish), [StepSize*(j-1), StepSize*j]);
    foo(end+1) = SpaceAnal(out);
end

idx = 6;
for j = 1:3 % 0 to 900
    out = SpaceCorps(cave(idx), 1:length(cave(idx).fish), [StepSize*(j-1), StepSize*j]);
    foo(end+1) = SpaceAnal(out);
end

idx = 7;
for j = 1:3 % 0 to 900
    out = SpaceCorps(cave(idx), 1:length(cave(idx).fish), [StepSize*(j-1), StepSize*j]);
    foo(end+1) = SpaceAnal(out);
end

idx = 8;
for j = 1:4 % 0 to 1200
    out = SpaceCorps(cave(idx), 1:length(cave(idx).fish), [StepSize*(j-1), StepSize*j]);
    foo(end+1) = SpaceAnal(out);
end

idx = 10;
for j = 1:3 % 0 to 900
    out = SpaceCorps(cave(idx), 1:length(cave(idx).fish), [StepSize*(j-1), StepSize*j]);
    foo(end+1) = SpaceAnal(out);
end

idx = 11;
for j = 1:3 % 0 to 900
    out = SpaceCorps(cave(idx), 1:length(cave(idx).fish), [StepSize*(j-1), StepSize*j]);
    foo(end+1) = SpaceAnal(out);
end

idx = 12;
for j = 1:3 % 0 to 900
    out = SpaceCorps(cave(idx), 1:length(cave(idx).fish), [StepSize*(j-1), StepSize*j]);
    foo(end+1) = SpaceAnal(out);
end

idx = 13;
for j = 1:4 % 0 to 900
    out = SpaceCorps(cave(idx), 1:length(cave(idx).fish), [StepSize*(j-1), StepSize*j]);
    foo(end+1) = SpaceAnal(out);
end

idx = 14;
for j = 1:3 % 0 to 900
    out = SpaceCorps(cave(idx), 1:length(cave(idx).fish), [StepSize*(j-1), StepSize*j]);
    foo(end+1) = SpaceAnal(out);
end


for j=1:length(foo)
    
    CfishArea = [CfishArea, foo(j).AREAreal];
    Coverlaps = [Coverlaps, foo(j).overlaps];
    CpairAreas = [CpairAreas, foo(j).sumsareas];
    
end

%% Surface analysis

SfishArea = []; 
Soverlaps = [];
SpairAreas = [];

idx = 1;
fishidx = [1 2 3 4 5 6  8  10  12];
for j = 1:2 % 0 to 600
    out = SpaceCorps(srf(idx), fishidx, [StepSize*(j-1), StepSize*j]);
    oof(j) = SpaceAnal(out);
end

idx = 2;
fishidx = [1 2 3 4 5 6 7 8 9 10  12 13 14 15 16 17   20 21];
for j = 1:3 % 0 to 900
    out = SpaceCorps(srf(idx), fishidx, [StepSize*(j-1), StepSize*j]);
    oof(end+1) = SpaceAnal(out);
end

idx = 3;
fishidx = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21  23 24 25 26 27  29 30 31];
for j = 1:2 % 0 to 600
    out = SpaceCorps(srf(idx), fishidx, [StepSize*(j-1), StepSize*j]);
    oof(end+1) = SpaceAnal(out);
end

idx = 4;
fishidx = [1 2 3 4 5 6 7 8 9 10  12 13 14 15  17 18 19  21 22];
for j = 1:2 % 0 to 600
    out = SpaceCorps(srf(idx), fishidx, [StepSize*(j-1), StepSize*j]);
    oof(end+1) = SpaceAnal(out);
end

idx = 5;
fishidx = [1 2  4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22  24];
for j = 1:3 % 0 to 900
    out = SpaceCorps(srf(idx), fishidx, [StepSize*(j-1), StepSize*j]);
    oof(end+1) = SpaceAnal(out);
end



for j=1:length(oof)
    
    SfishArea = [SfishArea, oof(j).AREAreal];
    Soverlaps = [Soverlaps, oof(j).overlaps];
    SpairAreas = [SpairAreas, oof(j).sumsareas];
    
end

%% Do the stats for the paper

fprintf('Mean %2.4f and STD %2.4f and N %i for CAVE overlap percentages \n', mean(Coverlaps ./ CpairAreas), std(Coverlaps ./ CpairAreas), length(Coverlaps));
fprintf('Mean %2.4f and STD %2.4f and N %i for SURFACE overlap percentages \n', mean(Soverlaps ./ SpairAreas), std(Soverlaps ./ SpairAreas), length(Soverlaps));

% [aa,bb,cc,dd] = ttest2(Coverlaps, Soverlaps)
[aa,bb,cc,dd] = ttest2(Coverlaps ./ CpairAreas, Soverlaps ./ SpairAreas)
