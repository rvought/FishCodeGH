function [out, alldFs] = dFanalysis(data)
% January 2020 
% load SurfaceDataRev2018a.mat and/or CaveDataRev2018a.mat
% This cycles through each entry in the data and extracts pairwise dFs

minimuminteraction = 50; % The minimum number of shared sample times to be considered for dF analysis
alldFs = [];
goodidxs = [];

for kk=1:length(data)
    if length(data(kk).fish) > 1
        goodidxs = [goodidxs kk];
    end
end

%for kk = 3:3 % for testing - I happen to know cave(3) is a good sample
for kk = goodidxs(end:-1:1)    
    
    numfish = length(data(kk).fish); % How many fish in this recording
        
    if numfish > 1 % If we have more than one fish
    
        combos = combnk(1:numfish, 2); % All pairwise combinations of fish
        
        for p = length(combos):-1:1 % For each pair of fish
            
            pair(p).fishnums = combos(p,:); % Save the output combo
            
            % DO BOTH FISH OVERLAP IN TIME?
                          
            clear FD sharedidxs Descartes dF ia ib aXY bXY idxa idxb;                

            % Get good data indices for both fish (remove NaNs for fish EOD frequency)
            
            idxa = find(~isnan(data(kk).fish(combos(p,1)).freq(:,2)));
            idxb = find(~isnan(data(kk).fish(combos(p,2)).freq(:,2)));
                        
            sharedidxs = intersect(idxa, idxb);
            
            if ~isempty(sharedidxs) % FISH MIGHT BE INTERACTING!!!!!!
                
                if length(sharedidxs) > minimuminteraction
               
                    % Calculate the dF
                                            
                        dF = abs(data(kk).fish(combos(p,1)).freq(sharedidxs,2) - data(kk).fish(combos(p,2)).freq(sharedidxs,2));
                        % dF = abs(data(kk).fish(combos(p,1)).freq(ia,2) - data(kk).fish(combos(p,2)).freq(ib,2));

                    % How far apart are they?
                    
                        aXY(1,:) = data(kk).fish(combos(p,1)).x(sharedidxs);
                        aXY(2,:) = data(kk).fish(combos(p,1)).y(sharedidxs);
                        bXY(1,:) = data(kk).fish(combos(p,2)).x(sharedidxs);
                        bXY(2,:) = data(kk).fish(combos(p,2)).y(sharedidxs);
                        
                        Descartes = pdist2(aXY', bXY'); % Euclidian distances
                        Descartes = Descartes(eye(length(Descartes)) == 1); % Extract data from matrix (identity only)
                                        
                    pair(p).dFmean = nanmean(dF);
                    pair(p).dFvar = nanvar(dF);
                    pair(p).meanDist = nanmean(Descartes);
                    pair(p).varDist = nanvar(Descartes);

                    pair(p).sharedidxs = sharedidxs;
                    pair(p).sharedtims = data(kk).fish(combos(p,1)).freq(sharedidxs,1);
                    pair(p).descartes = Descartes;
                    pair(p).dF = dF;

                            % maxlen = min([length(dF), length(Descartes)]);
                            if abs(length(dF) - length(Descartes)) > 0; fprintf('Yowza!'); end
                    FD(:,1) = dF;
                    FD(:,2) = Descartes;
                    [pair(p).covDistdF, pair(p).covDistdFpval] = corrcoef(FD);
                    
                end % If we have sufficient samples of interaction
                
            end % Fish might be interacting
           
                       
        end % For this pair of fish
        
        
        
    end % If we have more than one fish
    
    % Now we concatenate for plotting
    
    for j = 1:length(pair)
        alldFs = [alldFs pair(j).dF']; 
    end
    
        out(kk).pair = pair;
    
end % For every data entry


