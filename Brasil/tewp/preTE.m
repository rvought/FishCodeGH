function out = preTE(data, fishidx)

if nargin < 2
    fishidx = 1:length(data.fish);
end

fishum = length(fishidx);
combos = combnk(fishidx, 2);

for p = length(combos):-1:1 % For each pair of fish

         out(p).fishnums = combos(p,:); % Save the identities of the fish to the output structure
 
         % THIS SECTION IS WASTEFUL: Copies original data into structure
         out(p).EODA = data.fish(combos(p,1)).freq(:,2);
         out(p).EODB = data.fish(combos(p,2)).freq(:,2);
         out(p).tim = data.fish(combos(p,2)).freq(:,1);
         
         % Use embedded function to calcuate distance and dF
         [currdist, currdF] = distdfcalc(data.fish(combos(p,1)), data.fish(combos(p,2)));
         
         out(p).distance = currdist;
         out(p).dF = currdF;

end


%% Embedded function
function [dist, dF] = distdfcalc(A, B)

for j=length(A(1).freq(:,1)):-1:1 % For every time step in the sample
   
    if ~isnan(A.freq(j,2)) && ~isnan(B.freq(j,2)) % If both have real data at the moment         
        
        % Calculate distance
        dist(j) = sqrt((A.x(j) - B.x(j)).^2 + (A.y(j) - B.y(j)).^2);
        
        % Calculate dF
        dF(j) = abs(A.freq(j,2) - B.freq(j,2));
        
    end
    
end

% Fill in missing data.  This is dangerous - need reality check.

        dist(dist == 0) = NaN;
        dist = fillmissing(dist, 'pchip');
        dF(dF == 0) = NaN;
        dF = fillmissing(dF, 'pchip');

end % End of distdfcalc embedded function

end % End of preTE function