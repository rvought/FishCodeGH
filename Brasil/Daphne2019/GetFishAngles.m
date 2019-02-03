function out = GetFishAngles(in, tims)
% out = GetFishAngles(in, tims) where 'in' is the Ravi/Manu stucture and 
% 'tims' is time range in seconds (e.g. [100 500]).
% This takes Grid data and extracts the angles of Movement of the fish
% Uses only data with sequential samples - skips the skips.
%

% If user fails to give a range, use the entire sample
if nargin < 2
    tims = [0 999999];
end

% tr is the data within our time range
    tr = in.fish(1).freq(:,1) > tims(1) & in.fish(1).freq(:,1) < tims(2);
    

numfish = length(in.fish); % How many fish in this recording

for j = numfish:-1:1 % For each fish
    
    % Get valid data (check frequency data)            
    tmp(j).valididx = find(~isnan(in.fish(j).freq(tr,2)));
    
    % Check for sequential data only
    seqdat = find(diff(tmp(j).valididx) == 1);
    
    % Cycle through sequential data
    
    for k = 1:length(seqdat)
                        
        %  Change in X and Y
        dX = in.fish(j).x(seqdat(k)) - in.fish(j).x(seqdat(k)+1); 
        dY = in.fish(j).y(seqdat(k)) - in.fish(j).y(seqdat(k)+1);

        % Angle of movement
        out(j).realang(seqdat(k)) =  atan2(dY,dX);
        
    end
end