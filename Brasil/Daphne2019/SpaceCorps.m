function [out, cmbs] = SpaceCorps(in, casu, tims)
% Usage out = SpaceCorps(in)
% Find spatial relations between fish
% Distance histograms, range overlap

%% Setup

if nargin < 3 % If the user did not specify the time frame, take the whole thing
    tims = [0 999999];
end

% tr is the data within our time range
    tr = in.fish(1).freq(:,1) > tims(1) & in.fish(1).freq(:,1) < tims(2);
    
% Constrains change in angle of random fish
    constrainer = pi/4; % pi/2 (90 degrees) works great
    
% OPTION 1: Loop to find bounding box fitted to the data
xminedge = []; xscale = [];
yminedge = []; yscale = [];

    for bb = 1:length(in.fish)
            xminedge = min([xminedge, min(in.fish(bb).x)]);
            yminedge = min([yminedge, min(in.fish(bb).y)]);
            xscale = max([xscale, max(in.fish(bb).x)]);
            yscale = max([yscale, max(in.fish(bb).y)]);        
    end    
    
    xscale = xscale + abs(xminedge);
    yscale = yscale + abs(yminedge);

% OPTION 2: Use preset bounding boxes    
% if casu == 1 % cave
%     xminedge = -150; xscale = 400; 
%     yminedge = -150; yscale = 275;
% end
% if casu == 2 % surface
%     xminedge = -150; xscale = 300;  
%     yminedge = -100; yscale = 275;
% end

% Spatial centers for histogram
    ctrs{1} = xminedge-20:10:xscale+xminedge+20;
    ctrs{2} = yminedge-20:10:yscale+yminedge+20;

% Distance centers for histogram
    dctrs = 1:10:500;

    
    
numfish = length(in.fish); % How many fish in this recording

%% Cycle through each fish

for j = numfish:-1:1 % For each fish
    
    % Record the size of the box for analysis
    out(j).xbounds = [xminedge xminedge+xscale];
    out(j).ybounds = [yminedge yminedge+yscale];
    
    % Get valid data (check frequency data)            
    out(j).valididx = find(~isnan(in.fish(j).freq(tr,2)));
    
    %% Generate corresponding random fish
        % Starts at a random spot within the grid
        %  rnd(j).xy(:,out(j).valididx(1)) = [rand(1,1)*xscale rand(1,1)*yscale];
        % Starts with random angle
        
        % Starts at the same spot as the real fish in the grid but with
        % jiggled angle
        out(j).rndXY(:,out(j).valididx(1)) = [in.fish(j).x(out(j).valididx(1)) in.fish(j).y(out(j).valididx(1))];
            firstheta = atan2(in.fish(j).y(out(j).valididx(2)) - in.fish(j).y(out(j).valididx(1)), in.fish(j).x(out(j).valididx(2)) - in.fish(j).x(out(j).valididx(1)));
            deltatheta = constrainer * (rand(1,1) - 0.5); % Set a random change in direction for our artificial fish
        out(j).randtheta(1) = firstheta + deltatheta;
        
    % Moves the same distance as the real fish, but with jiggled change in direction
    for rr = 2:length(out(j).valididx)
            tmpXY(1,:) = [in.fish(j).x(out(j).valididx(rr-1)), in.fish(j).y(out(j).valididx(rr-1))];
            tmpXY(2,:) = [in.fish(j).x(out(j).valididx(rr)), in.fish(j).y(out(j).valididx(rr))];
            out(j).realhowfar(rr) = pdist(tmpXY); % How far did the real fish travel?

            deltatheta = constrainer * (rand(1,1) - 0.5); % Set a random change in direction for our artificial fish
            
            if out(j).valididx(rr) - out(j).valididx(rr-1) == 1 % Contiguous data
                
                % OPTION 1: Use the previous random theta as seed for new angle
                % out(j).randtheta(rr) = out(j).randtheta(rr-1) + deltatheta;

                % OPTION 2: Jiggle the fish theta
                realtheta = atan2(tmpXY(2,2) - tmpXY(1,2), tmpXY(2,1) - tmpXY(1,1));            
                out(j).randtheta(rr) = realtheta + deltatheta;
                
                
                
            else % We had a gap in the data  
                
                % Calculate the last real fish angle as seed for new angle
                realtheta = atan2(tmpXY(2,2) - tmpXY(1,2), tmpXY(2,1) - tmpXY(1,1));            
                out(j).randtheta(rr) = realtheta + deltatheta;

            end

            % FINALLY - Movement with the same distance as real but jiggled angle        
            out(j).rndXY(:,out(j).valididx(rr)) = [out(j).rndXY(1,out(j).valididx(rr-1)) + out(j).realhowfar(rr)*cos(out(j).randtheta(rr)), out(j).rndXY(2,out(j).valididx(rr-1)) + out(j).realhowfar(rr)*sin(out(j).randtheta(rr))]; 
                  
            % BUT... we don't want to randomly escape from the grid, so if we
            % violate the boundaries, we 'reflect' back in
            if out(j).rndXY(1,end) < xminedge
                out(j).rndXY(1,end) = xminedge + xminedge - out(j).rndXY(1,end);
            end
            if out(j).rndXY(1,end) > xminedge + xscale
                out(j).rndXY(1,end) = (xminedge + xscale) - (out(j).rndXY(1,end) - (xminedge + xscale));
            end
            if out(j).rndXY(2,end) < yminedge
                out(j).rndXY(2,end) = yminedge + yminedge - out(j).rndXY(2,end);
            end
            if out(j).rndXY(2,end) > yminedge + yscale
                out(j).rndXY(2,end) = (yminedge + yscale) - (out(j).rndXY(2,end) - (yminedge + yscale));
            end
            
            % And a copy of the original data, just for fun
            out(j).realXY(:,out(j).valididx(rr)) = [in.fish(j).x(out(j).valididx(rr)), in.fish(j).y(out(j).valididx(rr))]; 

            % Completely random points
            % rnd(j).xy(:,out(j).valididx(rr)) = [rand(1,1)*xscale rand(1,1)*yscale]; 
    end
    
    %% Get statistics for each real and random fish
    
    tmp(j).xy(1,:) = in.fish(j).x(out(j).valididx);
    tmp(j).xy(2,:) = in.fish(j).y(out(j).valididx);  
    
    rtmp(j).xy(1,:) = out(j).rndXY(1,out(j).valididx);
    rtmp(j).xy(2,:) = out(j).rndXY(2,out(j).valididx);  

    % Spatial histogram
    out(j).realhist = hist3(tmp(j).xy', ctrs);
    out(j).randhist = hist3(rtmp(j).xy', ctrs);
    
    % Prepare for the allhist later
    out(j).allhist = zeros(1,length(dctrs));  
    out(j).allrandhist = zeros(1,length(dctrs));  
    
end

if numfish > 1 % We have more than one fish

    % Compare spatial histograms (fish against all other fish)
    %
    %
    %
    
    combos = combnk(1:numfish, 2); % All pairwise combinations of fish

    for p = length(combos):-1:1 % For each pair of fish

        cmbs(p).fishnums = combos(p,:); % Save the output combo
        
        % Get only entries for which we have data for both fish
        [~, idx1, idx2] = intersect(out(combos(p,1)).valididx, out(combos(p,2)).valididx);
               
        % Calculate the interfish-distance
        cmbs(p).dist = [];
        for jj = 1:length(idx1) % THERE HAS GOT TO BE A MORE EFFICIENT WAY OF DOING THIS...
            cmbs(p).dist(out(combos(p,1)).valididx(idx1(jj))) = pdist2(tmp(combos(p,1)).xy(:,idx1(jj))', tmp(combos(p,2)).xy(:,idx2(jj))');            
            cmbs(p).randdist(out(combos(p,1)).valididx(idx1(jj))) = pdist2(tmp(combos(p,1)).xy(:,idx1(jj))', tmp(combos(p,2)).xy(:,idx2(jj))');

        end
        
        % Distance histogram for each pair of fish
        if length(cmbs(p).dist) > 1
            cmbs(p).dhist = hist(cmbs(p).dist, dctrs);
        else
            cmbs(p).dhist = [];
        end
        
        % Assemble the histogram of distances of each fish to all others
        
        if sum(cmbs(p).dhist) > 0
            out(combos(p,1)).allhist = out(combos(p,1)).allhist + cmbs(p).dhist;        
            out(combos(p,2)).allhist = out(combos(p,2)).allhist + cmbs(p).dhist;        
        end
        
        % Compare spatial histograms (fish against each fish separately)
        
        
            
        

    end % Cycle through each pair of fish


%% Plots 

figure(casu); clf;

ax(1) = subplot(121); hold on;
    for z=1:length(in.fish)
        plot(in.fish(z).x(out(z).valididx), in.fish(z).y(out(z).valididx), '.', 'MarkerSize', 8);
%        plot(in.fish(z).x(out(z).valididx), in.fish(z).y(out(z).valididx), '*-');
    end
ax(2) = subplot(122); hold on;
    for z=1:length(out)
        plot(out(z).rndXY(1,out(z).valididx), out(z).rndXY(2,out(z).valididx), '.', 'MarkerSize', 8);
%        plot(out(z).rndXY(1,out(z).valididx), out(z).rndXY(2,out(z).valididx), '*-');
    end
linkaxes(ax, 'xy');    



end % If we have more than one fish


