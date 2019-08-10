function foo = SpaceAnal(out)


%% Get the box around the entire group of fishies and calculate hulls for each
    minX = [];
    minY = [];

for j=length(out):-1:1 % For every fish
    
    if length(out(j).realXY) > 1 % If there is data
        
    % Eliminate zeros, which occur when the fish was not tracked    
    realXY(1,:) = out(j).realXY(1,(find(out(j).realXY(1,:) ~= 0)));
    realXY(2,:) = out(j).realXY(2,(find(out(j).realXY(2,:) ~= 0)));
    
    samplelength(j) = length(out(j).realXY); % The length of each data set
    
    % Get the minX and minY so that we can move into an all-positive regime
    % for the poly2mask functiob below
    
    minX = min([minX, min(realXY(1,:))]);
    minY = min([minY, min(realXY(2,:))]);
    
    clear realXY;
    
    end
end

    minX = 1 + abs(minX);
    minY = 1 + abs(minY);
    minsamplength = round(0.75 * max(samplelength)); % Not fair to calculate an area for a fish with few data points

%% Generate the masks and calculate the areas and overlaps    

for j=length(out):-1:1 % For each fish
    
    if length(out(j).realXY) > minsamplength % Make sure that we have sufficient data to calculate an area
        
        % Again, eliminate zeros, which occur when the fish was not tracked (I'm a very lazy coder)    
        realXY(1,:) = out(j).realXY(1,(out(j).realXY(1,:) ~= 0));
        realXY(2,:) = out(j).realXY(2,(out(j).realXY(2,:) ~= 0));
        
        % Get the convex hull for each fish
        [foo.IDXreal{j}, foo.AREAreal(j)] = convhull(realXY(1,:)+minX, realXY(2,:)+minX);
        
        clear realXY;

    end

end

%% Calculate overlaps for all fish pairs

% Get each pairwise combination of fish in the sample
combos = combnk(1:length(out), 2);

    foo.overlaps = []; foo.sumsareas = []; 
    foo.A = []; foo.B = [];
    boxsize = 400;  % MAKE SURE THAT THIS IS LARGE ENOUGH SO FISH DON'T ESCAPE!!!

% For each pair...    
for j = 1:length(combos)

    a = combos(j,1); b = combos(j,2);
    
    if length(out(a).realXY) > minsamplength && length(out(b).realXY) > minsamplength % AGAIN, check to make sure we have data...
    
% figure(1);  
% subplot(131); imshow(poly2mask(out(a).realXY(1,IDXreal{combos(j,1)})+minX,out(a).realXY(2,IDXreal{combos(j,1)})+minY,boxsize,boxsize));
% subplot(132); imshow(poly2mask(out(b).realXY(1,IDXreal{combos(j,2)})+minX,out(b).realXY(2,IDXreal{combos(j,2)})+minY,boxsize,boxsize));
% subplot(133); imshow(poly2mask(out(a).realXY(1,IDXreal{combos(j,1)})+minX,out(a).realXY(2,IDXreal{combos(j,1)})+minY,boxsize,boxsize) & poly2mask(out(b).realXY(1,IDXreal{combos(j,2)})+minX,out(b).realXY(2,IDXreal{combos(j,2)})+minY,boxsize,boxsize));

% Overlapping region
foo.overlaps(end+1) = bwarea(poly2mask(out(a).realXY(1,:)+minX,out(a).realXY(2,:)+minY,boxsize,boxsize) & ...
    poly2mask(out(b).realXY(1,:)+minX,out(b).realXY(2,:)+minY,boxsize,boxsize));

% Combined region
foo.sumsareas(end+1) = bwarea(poly2mask(out(a).realXY(1,:)+minX,out(a).realXY(2,:)+minY,boxsize,boxsize) + ...
    poly2mask(out(b).realXY(1,:)+minX,out(b).realXY(2,:)+minY,boxsize,boxsize));

% Each by themselves

foo.A(end+1) = bwarea(poly2mask(out(a).realXY(1,:)+minX,out(a).realXY(2,:)+minY,boxsize,boxsize));
foo.B(end+1) = bwarea(poly2mask(out(b).realXY(1,:)+minX,out(b).realXY(2,:)+minY,boxsize,boxsize));

    end 

end




