function out = caveposter(data, freqthresh, durthresh)
% Usage out = analretent(data, freqthresh, durthresh)
% data is the output of catfish (e.g. s27)
% freqthresh has to do with the framerate - 15 (Hz) is a good number for
% this as the usual sample rate is just under 20 Hz.
% durthresh is the number of consecutive samples to qualify as a good
% segment - 100 is a good number here.

numcolors = 30;
oclrs = hsv(numcolors); % Colors for plotting up to numcolors different fishes


%% Extract only the good position data

fprintf('Finding sections of cleanly sampled position data. \n'); 

for trk = length(data):-1:1;
    
tsteps = diff(data(trk).tt); % Get the intervals

    % Construct that gives us the starts and ends of useful data 
    z = zeros(1,length(tsteps)); % A list of zeros 
    z((1 ./ tsteps) > freqthresh) = 1; % Proper intervals are marked as 1.
    
    tmp(trk).strts = find(diff(z) == 1); % From 0 to 1
    tmp(trk).eds = find(diff(z) == -1); % From 1 back down to 0

    if ~isempty(tmp(trk).eds) && ~isempty(tmp(trk).strts)
    
    % Cleanup the ends - we always want pairs of starts/stops
    if (tmp(trk).eds(1) < tmp(trk).strts(1)); tmp(trk).strts = [1 tmp(trk).strts]; end; % If the first event is an end, then we have good data from the start.
    if (tmp(trk).strts(end) > tmp(trk).eds(end)); tmp(trk).eds(end+1) = length(z); end; % If the last event is a start, then we have good data to the end.

    end;
    
end

fprintf('Finding long-duration sections in the position data. \n'); 

for j = length(tmp):-1:1;
    
    durs = tmp(j).eds - tmp(j).strts;
    out(j).eds = tmp(j).eds(durs > durthresh);
    out(j).strts = tmp(j).strts(durs > durthresh);

end

clear tmp;

figure(1); clf; hold on;

% Plot the nice segments for each fish in Figure 1
for j=1:length(out); 
    if ~isempty(out(j).strts)
    for k=1:length(out(j).strts);
        plot(data(j).tx(out(j).strts(k):out(j).eds(k)), data(j).ty(out(j).strts(k):out(j).eds(k)), 'Color', oclrs(j,:));
    end
    end
end

%% Distance statistics

% Calculate the distances traveled

fprintf('Calculating velocities. \n'); 

for j=1:length(out); % For each fish

    if ~isempty(out(j).strts)
    
    for k = 1:length(out(j).strts) % For each clean section of data
       
        for p = out(j).strts(k):out(j).eds(k)-1 % Each time step
            % Calculate the distance for each time stamp
            out(j).trk(k).dist(p) = pdist([data(j).tx(p), data(j).ty(p); data(j).tx(p+1), data(j).ty(p+1)]);
            
        end
        
        out(j).trk(k).totaldist = sum(out(j).trk(k).dist); % Sum up the total distance
        out(j).trk(k).totaltime = data(j).tt(out(j).eds(k)) - data(j).tt(out(j).strts(k));
        out(j).trk(k).meanvel = sum(out(j).trk(k).dist)/out(j).trk(k).totaltime; % Calculate the mean velocity per sample
        
    end
    
    out(j).mmeanvel = mean([out(j).trk(:).meanvel]);

    end
    
end

%% Frequency statistics

fprintf('Finding sections of cleanly sampled frequency data. \n'); 

for trk = length(data):-1:1;

tsteps = diff(data(trk).tim); % Get the intervals

    % Construct that gives us the starts and ends of useful data 
    z = zeros(1,length(tsteps)); % A list of zeros 
    z((1 ./ tsteps) > freqthresh) = 1; % Proper intervals are marked as 1.
    
    tmp(trk).strts = find(diff(z) == 1); % From 0 to 1
    tmp(trk).eds = find(diff(z) == -1); % From 1 back down to 0

    % Cleanup the ends - we always want pairs of starts/stops
    
        if ~isempty(tmp(trk).eds) && ~isempty(tmp(trk).strts)

            if (tmp(trk).eds(1) < tmp(trk).strts(1)); tmp(trk).strts = [1 tmp(trk).strts]; end; % If the first event is an end, then we have good data from the start.
            if (tmp(trk).strts(end) > tmp(trk).eds(end)); tmp(trk).eds(end+1) = length(z); end; % If the last event is a start, then we have good data to the end.

        end
end

for j = length(tmp):-1:1;
    
    durs = tmp(j).eds - tmp(j).strts;
    out(j).Feds = tmp(j).eds(durs > durthresh);
    out(j).Fstrts = tmp(j).strts(durs > durthresh);

end

for j=1:length(out); % For each fish

    tmpfreq = [];

    if ~isempty(out(j).Fstrts)
        
    for k = 1:length(out(j).Fstrts) % For each clean section of data
        
        % Calculate the mean frequency
        out(j).trk(k).meanfreq = mean(data(j).freq(out(j).Fstrts(k):out(j).Feds(k)));
        % Concatenate the frequency 
        tmpfreq = [tmpfreq data(j).freq(out(j).Fstrts(k):out(j).Feds(k))];
        
    end
    
    out(j).mmeanfreq = mean([out(j).trk(:).meanfreq]);
    out(j).stdfreq = std(tmpfreq);
    
    end
    
end




