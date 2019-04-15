function QuickPlotGrid(in, idx, freqs, tims)
% Usage: QuickPlotGrid(in, freqs, tims)
% Generates F1 with frequency traces and F2 with grid positions
% freqs and tims are optional. 
% Freqs in Hz, tims in seconds.


%% Setup
if exist idx == 0
    idx = 1:length(in.fish); % Default is show all fish
end    
if nargin < 3 
    freqs = [250 500]; % Default frequency range for EODs
end    

    freqs = sort(freqs);

    figure(1); clf; hold on; % Frequency plot
    figure(2); clf; hold on; % Grid position plot
    % Make the Grid position plot square
        wp=get(gcf,'Position');
        set(gcf, 'Position', [wp(1),wp(2),560,560])


%% For each fish in the sample, plot frequency and xy plots
for j = 1:length(idx)

        tt = find(~isnan(in.fish(idx(j)).freq(:,2))); % all valid data
        
        if nargin == 4 % use only data in the range specified by the user
            tims = sort(tims);
            uu = find(in.fish(idx(j)).freq(:,1) > tims(1) & in.fish(idx(j)).freq(:,1) < tims(2));
            tt = intersect(tt, uu);
        end
   
   figure(1); plot(in.fish(idx(j)).freq(tt,1), in.fish(idx(j)).freq(tt,2), '.', 'MarkerSize', 8); ylim(freqs);
   figure(2); plot(in.fish(idx(j)).x(tt), in.fish(idx(j)).y(tt), '.', 'MarkerSize', 8); 
        axis([-200, 250, -200, 250]); % May need adjustment
    
end
