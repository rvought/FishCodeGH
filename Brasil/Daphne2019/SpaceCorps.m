function [out, cmbs] = SpaceCorps(in, casu, feesh, tims)
% Usage out = SpaceCorps(in)
% Find spatial relations between fish and compare to 'jiggled' and 'randomized'
% distributions.
% Distance histograms, range overlap

%% Setup

% if nargin < 4 % If the user did not specify the time frame, take the whole thing
%     tims = [0 999999];
% end

% if nargin < 3 % If the user did not specify specific fish, use them all
%     feesh = 1:length(in.fish);
% end


    
% Constrains change in angle for jiggling and for the randomized fish
    constrainer = pi/2; % pi/2 (90 degrees) works great (converted to +/- below)
    
% OPTION 1: Loop to find bounding box fitted to the spatial extent of data
xminedge = []; xscale = [];
yminedge = []; yscale = [];

    for bb = 1:length(feesh)
            xminedge = min([xminedge, min(in.fish(feesh(bb)).x)]);
            yminedge = min([yminedge, min(in.fish(feesh(bb)).y)]);
            xscale = max([xscale, max(in.fish(feesh(bb)).x)]);
            yscale = max([yscale, max(in.fish(feesh(bb)).y)]);        
    end    
    
    xscale = xscale + abs(xminedge);
    yscale = yscale + abs(yminedge);

% OPTION 2: Use preset bounding boxes set by us   
% if casu == 1 % cave
%     xminedge = -150; xscale = 400; 
%     yminedge = -175; yscale = 275;
% end
% if casu == 2 % surface
%     xminedge = -150; xscale = 300;  
%     yminedge = -100; yscale = 275;
% end

% Spatial centers for 2D histogram (XXcm boxes that extend 20cm beyond edges)
% The boxes are identical size, although the overall region (defined by
% minedge and scale values) are not the same from recording to recording.
% howmanybins = 8;
% xbinsiz = (xscale+xminedge+20 + abs(xminedge-20)) / howmanybins;
% ybinsiz = (yscale+yminedge+20 + abs(yminedge-20)) / howmanybins;
%     ctrs{1} = xminedge-20:xbinsiz:xscale+xminedge+20;
%     ctrs{2} = yminedge-20:ybinsiz:yscale+yminedge+20;
% 
%     qwer = ctrs{1}
%     asdf = ctrs{2}

% THE CRITICAL PARAMETER (sort of...)
% Use standard centers for all of the data (extend well beyond limits of
% possible fish positions)

BoxLen = 50; % We've tried 50, 60, 75, 100. This is in cm.
    ctrs{1} = -300:BoxLen:300;
    ctrs{2} = -300:BoxLen:300;
    
% Distance bins between pairs of fish for histogram
    dctrs = 1:10:500;    
    
numfish = length(feesh); % How many fish in this recording

%% Cycle through each fish

for j = numfish:-1:1 % For each fish    
    
    % Record the size of the box for analysis (Calculated above)
    out(j).xbounds = [xminedge xminedge+xscale];
    out(j).ybounds = [yminedge yminedge+yscale];
    
    % Record the centers of the bins for analysis (Determined above)
    out(j).ctrs = ctrs;
    
    % tr is the data within our time range
        tr = in.fish(j).freq(:,1) > tims(1) & in.fish(j).freq(:,1) < tims(2);
    % Get valid data (check frequency data to omit NaNs)            
    out(j).valididx = tr(~isnan(in.fish(feesh(j)).freq(tr,2)));
    
    if ~isempty(out(j).valididx) % Make sure that we have data before proceding
        
    %% Generate corresponding jiggled and random fish
        
        % Both random and jiggled start at the same spot as the real fish in the grid 
        
        out(j).rndXY(:,out(j).valididx(1)) = [in.fish(feesh(j)).x(out(j).valididx(1)) in.fish(feesh(j)).y(out(j).valididx(1))]; % First XY for randome

        % Both random and jiggled fish have a jiggled initial theta
        
            firstheta = atan2(in.fish(feesh(j)).y(out(j).valididx(2)) - in.fish(feesh(j)).y(out(j).valididx(1)), in.fish(feesh(j)).x(out(j).valididx(2)) - in.fish(feesh(j)).x(out(j).valididx(1)));
            deltatheta = constrainer * (rand(1,1) - 0.5); % Set a random change in direction for our artificial fish

        out(j).randtheta(1) = firstheta + deltatheta; % First theta for Random
        
        % Copy into the jiggle data
        out(j).jigXY = out(j).rndXY; 
        out(j).jiggletheta(1) = out(j).randtheta(1);
        
    % Cycle for every valid data point in original recording
    % Both jiggled and random fish move the same distance as the real fish, but with jiggled change in direction
    for rr = 2:length(out(j).valididx)
        
            % Put the fish positions into a convenient format and calculate
            % the distance and angle of the real fish 
            tmpXY(1,:) = [in.fish(feesh(j)).x(out(j).valididx(rr-1)), in.fish(feesh(j)).y(out(j).valididx(rr-1))];
            tmpXY(2,:) = [in.fish(feesh(j)).x(out(j).valididx(rr)), in.fish(feesh(j)).y(out(j).valididx(rr))];
            out(j).realhowfar(rr) = pdist(tmpXY); % How far did the real fish travel?
            realtheta = atan2(tmpXY(2,2) - tmpXY(1,2), tmpXY(2,1) - tmpXY(1,1));            

            deltatheta = constrainer * (rand(1,1) - 0.5); % Set a random change in direction for our artificial fish
            
            if out(j).valididx(rr) - out(j).valididx(rr-1) == 1 % Contiguous data
                
                % Variant 1: Random - Use the previous random theta as seed for new angle
                out(j).randtheta(rr) = out(j).randtheta(rr-1) + deltatheta;

                % Variant 2: Jiggled - Use the fish theta as seed for new angle
                out(j).jiggletheta(rr) = realtheta + deltatheta;
                
            else % We had a gap in the data  
                
                % The "GAP" calculation isn't technically necessary but was
                % an issue I explored.
                % Calculate the last real fish angle as seed for new angle
                % We do this for both the random and jiggled data - This
                % has very little impact on the outcomes. 
                
                out(j).randtheta(rr) = realtheta + deltatheta;
                out(j).jiggletheta(rr) = realtheta + deltatheta;

            end

            % FINALLY - Movement with the same distance as real but with altered angles (next XY position)
            
            out(j).jigXY(:,out(j).valididx(rr)) = [out(j).jigXY(1,out(j).valididx(rr-1)) + out(j).realhowfar(rr)*cos(out(j).jiggletheta(rr)), out(j).jigXY(2,out(j).valididx(rr-1)) + out(j).realhowfar(rr)*sin(out(j).jiggletheta(rr))]; 
            out(j).rndXY(:,out(j).valididx(rr)) = [out(j).rndXY(1,out(j).valididx(rr-1)) + out(j).realhowfar(rr)*cos(out(j).randtheta(rr)), out(j).rndXY(2,out(j).valididx(rr-1)) + out(j).realhowfar(rr)*sin(out(j).randtheta(rr))]; 
                  
            % BUT... we don't want to wander out of the grid, so if we
            % violate the boundaries, we 'reflect' back in.  This is
            % extremely rudimentary - there are many more elegant
            % solutions.
            
            % Catch 'random' exits
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
            
            % Catch 'jiggle' exits
            if out(j).jigXY(1,end) < xminedge
                out(j).jigXY(1,end) = xminedge + xminedge - out(j).jigXY(1,end);
            end
            if out(j).jigXY(1,end) > xminedge + xscale
                out(j).jigXY(1,end) = (xminedge + xscale) - (out(j).jigXY(1,end) - (xminedge + xscale));
            end
            if out(j).jigXY(2,end) < yminedge
                out(j).jigXY(2,end) = yminedge + yminedge - out(j).jigXY(2,end);
            end
            if out(j).jigXY(2,end) > yminedge + yscale
                out(j).jigXY(2,end) = (yminedge + yscale) - (out(j).jigXY(2,end) - (yminedge + yscale));
            end
            
            % And a copy of the original data is put into the structure, just for plotting fun
            out(j).realXY(:,out(j).valididx(rr)) = [in.fish(feesh(j)).x(out(j).valididx(rr)), in.fish(feesh(j)).y(out(j).valididx(rr))]; 

    end

    
    %% Prepare data for analysis for each real and jiggled and random fish

    % Put the data into a convenient format for analysis
    tmp(j).xy(1,:) = in.fish(feesh(j)).x(out(j).valididx);
    tmp(j).xy(2,:) = in.fish(feesh(j)).y(out(j).valididx);  
    
%     rtmp(j).xy(1,:) = out(j).rndXY(1,out(j).valididx);
%     rtmp(j).xy(2,:) = out(j).rndXY(2,out(j).valididx);  
% 
%     jtmp(j).xy(1,:) = out(j).jigXY(1,out(j).valididx);
%     jtmp(j).xy(2,:) = out(j).jigXY(2,out(j).valididx);  

    rtmp(j).xy(1,:) = out(j).rndXY(1,:);
    rtmp(j).xy(2,:) = out(j).rndXY(2,:);  

    jtmp(j).xy(1,:) = out(j).jigXY(1,:);
    jtmp(j).xy(2,:) = out(j).jigXY(2,:);  

    % Spatial histogram (Heat map from above) for each fish
    out(j).realhist = hist3(tmp(j).xy', ctrs);
    out(j).randhist = hist3(rtmp(j).xy', ctrs);
    out(j).jighist = hist3(jtmp(j).xy', ctrs);
    
    % Prepare for the allhist later
    out(j).allhist = zeros(1,length(dctrs));  
    out(j).allrandhist = zeros(1,length(dctrs));  
    out(j).alljighist = zeros(1,length(dctrs));  
    
    end % If we have data for the current fish

end % Cycle through each fish



%% Cycle through pairs of fish
% if numfish > 1 % We have more than one fish
%     
%     combos = combnk(1:length(out), 2); % All pairwise combinations of fish
% 
%     for p = length(combos):-1:1 % For each pair of fish
% 
%         cmbs(p).fishnums = combos(p,:); % Save the identities of the fish to the output structure.
%         
%         % Get only entries for which we have data for both fish (shared indices)
%         [~, idx1, idx2] = intersect(out(combos(p,1)).valididx, out(combos(p,2)).valididx);
%         if ~isempty(idx1) && ~isempty(idx2)       
%         % Calculate the interfish-distance
%         
%         for jj = 1:length(idx1) % THERE HAS GOT TO BE A MORE EFFICIENT WAY OF DOING THIS...
%             cmbs(p).realdist(out(combos(p,1)).valididx(idx1(jj))) = pdist2(tmp(combos(p,1)).xy(:,idx1(jj))', tmp(combos(p,2)).xy(:,idx2(jj))');            
%             cmbs(p).randdist(out(combos(p,1)).valididx(idx1(jj))) = pdist2(rtmp(combos(p,1)).xy(:,idx1(jj))', rtmp(combos(p,2)).xy(:,idx2(jj))');
%             cmbs(p).jigdist(out(combos(p,1)).valididx(idx1(jj))) = pdist2(jtmp(combos(p,1)).xy(:,idx1(jj))', jtmp(combos(p,2)).xy(:,idx2(jj))');
%         end
%         
%         % Distance histogram for each pair of fish
%         if length(cmbs(p).realdist) > 1 % Then we have data for this fish
%             cmbs(p).realhist = hist(cmbs(p).realdist, dctrs);
%             cmbs(p).randhist = hist(cmbs(p).randdist, dctrs);
%             cmbs(p).jighist = hist(cmbs(p).jigdist, dctrs);
%         else % We don't have data for this particular fish
%             fprintf('This happens. \n ');
%         end
%         
%         % Assemble the histogram of distances of each fish to all others
%         
%         if sum(cmbs(p).realhist) > 0
%             out(combos(p,1)).allhist = out(combos(p,1)).allhist + cmbs(p).realhist;        
%             out(combos(p,2)).allhist = out(combos(p,2)).allhist + cmbs(p).realhist;        
%             
%             out(combos(p,1)).allrandhist = out(combos(p,1)).allrandhist + cmbs(p).randhist;        
%             out(combos(p,2)).allrandhist = out(combos(p,2)).allrandhist + cmbs(p).randhist;        
%             
%             out(combos(p,1)).alljighist = out(combos(p,1)).alljighist + cmbs(p).jighist;        
%             out(combos(p,2)).alljighist = out(combos(p,2)).alljighist + cmbs(p).jighist;        
%                         
%         end
%         
%         end
%         
%         % Compare spatial histograms (fish against each fish separately)
%         
%             % Simple sums of overlap
%         
%         
%         
%     end % Cycle through each pair of fish


%% Plots 

figure(casu); clf;

ax(1) = subplot(131); hold on;
    for z=1:length(feesh)
        if ~isempty(out(z).valididx)
        plot(in.fish(feesh(z)).x(out(z).valididx), in.fish(feesh(z)).y(out(z).valididx), '.', 'MarkerSize', 8);
%        plot(in.fish(z).x(out(z).valididx), in.fish(z).y(out(z).valididx), '*-');
        end
    end
ax(2) = subplot(132); hold on;
    for z=1:length(out)
        if ~isempty(out(z).valididx)
        plot(out(z).jigXY(1,:), out(z).jigXY(2,:), '.', 'MarkerSize', 8);
%        plot(out(z).jigXY(1,out(z).valididx), out(z).jigXY(2,out(z).valididx), '*-');
        end
    end
ax(3) = subplot(133); hold on;
    for z=1:length(out)
        if ~isempty(out(z).valididx)
        plot(out(z).rndXY(1,:), out(z).rndXY(2,:), '.', 'MarkerSize', 8);
%        plot(out(z).rndXY(1,out(z).valididx), out(z).rndXY(2,out(z).valididx), '*-');
        end
    end
linkaxes(ax, 'xy');    

% figure(casu+2); clf;
%     axx(1) = subplot(131); hold on;  
%     for z=1:length(cmbs)
%         plot(cmbs(z).realhist);
%     end
%     axx(2) = subplot(132); hold on;  
%     for z=1:length(cmbs)
%         plot(cmbs(z).jighist);
%     end
%     axx(3) = subplot(133); hold on;  
%     for z=1:length(cmbs)
%         plot(cmbs(z).randhist);
%     end
% linkaxes(axx, 'xy');    
    
%end % If we have more than one fish
cmbs = 1;

%% Analyses





