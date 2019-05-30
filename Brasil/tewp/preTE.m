function [dFdist, orig] = preTE(data, fishidx)

if nargin < 2
    fishidx = 1:length(data.fish);
end

fishnum = length(fishidx); % How many fish
combos = combnk(fishidx, 2); % Each pair of fish that we plan to analyze

Fs = 1/mean(diff(data.fish(fishidx(1)).freq(:,1)));
[b,a] = butter(3, 0.05/(Fs/2), 'low');

%% Copy original data into structure
for k = fishnum:-1:1
    
         orig(fishidx(k)).EOD = data.fish(fishidx(k)).freq(:,2);
         orig(fishidx(k)).xy(:,1) = data.fish(fishidx(k)).x;
         orig(fishidx(k)).xy(:,2) = data.fish(fishidx(k)).y;
         orig(fishidx(k)).tim = data.fish(fishidx(k)).freq(:,1);
         
end

%% Calculate dF and distance between each pair
for p = length(combos):-1:1 % For each pair of fish

         dFdist(p).fishnums = combos(p,:); % Save the identities of the fish to the output structure
 
         
         % Use embedded function to calcuate distance and dF
         [currdist, currdF] = distdfcalc(data.fish(combos(p,1)), data.fish(combos(p,2)));
         
         dFdist(p).distance = filtfilt(b,a,currdist);
         dFdist(p).dF = filtfilt(b,a,currdF);
         dFdist(p).tim = data.fish(combos(p,1)).freq(:,1);

end


%% Embedded function to calculate dF and distance
function [dist, dF] = distdfcalc(A, B)

for j=length(A(1).freq(:,1)):-1:1 % For every time step in the sample
   
    if ~isnan(A.freq(j,2)) && ~isnan(B.freq(j,2)) % If both have real data at the moment         
        
        % Calculate distance
        dist(j) = sqrt((A.x(j) - B.x(j)).^2 + (A.y(j) - B.y(j)).^2);
        
        % Calculate dF
        dF(j) = abs(A.freq(j,2) - B.freq(j,2));
        
    end
    
end

% Fill in missing data.  This is dangerous - need reality check somewhere!!

        dist(dist == 0) = NaN;
        dist = fillmissing(dist, 'pchip');
        dF(dF == 0) = NaN;
        dF = fillmissing(dF, 'pchip');

end % End of distdfcalc embedded function

end % End of preTE function