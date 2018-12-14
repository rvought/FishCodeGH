function out = bCorr(data, rango)

startim = rango(1); endtim = rango(2);

if length(data) == 1
    figure(1); clf; 
    subplot(211); hold on; xlabel('Samples'); ylabel('Velocity');
    subplot(212); hold on; xlabel('Time, s'); ylabel('EOD Frequency Hz'); ylim([200 500]);
end
figure(2); clf;
subplot(121); hold on; xlabel('Mean EOD Frequency, Hz'); ylabel('EOD variance, Hz');
subplot(122); hold on; xlabel('Mean EOD Frequency, Hz'); ylabel('Mean velocity');

for j=1:length(data) % For each recording session
    
    numfish = length(data(j).fish); % How many fish in this recording
    
    %% DATUMS FOR EACH FISH
    
    for ff = length(data(j).fish):-1:1 

        % Current epoch.  freq(:,1) is timestamps, freq(:,2) is EOD freq
        cts = find(data(j).fish(ff).freq(:,1) > startim & data(j).fish(ff).freq(:,1) < endtim); 
        cts = cts(~isnan(data(j).fish(ff).freq(cts,2))); % only take non NaN data (ARGH!)
        
        
        % Mean and variance of frequency of this fish
        
        out(j).fish(ff).mFreq = mean(data(j).fish(ff).freq(cts,2));
        out(j).fish(ff).vFreq = var(data(j).fish(ff).freq(cts,2));
        out(j).fish(ff).nFreq = length(cts);

        % Instantaneous distance and velocity of the fish
        
        for vv = 2:length(cts)
            if cts(vv)-cts(vv-1) == 1 % ONLY TAKE DATA WHERE WE HAVE CONSECUTIVE SAMPLES!!!!
                X(1,1) = data(j).fish(ff).x(cts(vv-1)); % First X
                X(2,1) = data(j).fish(ff).y(cts(vv-1)); % First Y
                X(1,2) = data(j).fish(ff).x(cts(vv));   % Second X
                X(2,2) = data(j).fish(ff).y(cts(vv));   % Second Y
                
                dist(vv) = pdist(X);
                vel(vv) = dist(vv) * (data(j).fish(ff).freq(cts(vv),1) - data(j).fish(ff).freq(cts(vv-1),1));
                
            end
        end
                
        out(j).fish(ff).mVel = mean(vel);
        out(j).fish(ff).mfiltVel = mean(medfilt1(vel,5));
        out(j).fish(ff).totalDist = sum(dist);
        out(j).fish(ff).totalfiltDist = sum(medfilt1(dist,5));

        if length(data) == 1
        figure(1);
        subplot(211); plot(medfilt1(vel,5), 'LineWidth', 2);
        subplot(212); plot(data(j).fish(ff).freq(cts,1), data(j).fish(ff).freq(cts,2), 'LineWidth', 2)
        end
        
        figure(2);
        % EOD variability against mean EODf
        subplot(121); plot(out(j).fish(ff).mFreq, out(j).fish(ff).vFreq, '.', 'MarkerSize', 20)
        % Mean velocity against mean EODf
        subplot(122); plot(out(j).fish(ff).mFreq, out(j).fish(ff).mfiltVel, '.', 'MarkerSize', 20)

    end
    
    
    %% PAIRWISE INTERACTIONS BETWEEN FISH
    
    if numfish > 1 % If we have more than one fish
    
        combos = combnk(1:numfish, 2); % All pairwise combinations of fish
        
        for p = 1:length(combos) % For each pair of fish
            
            out(j).pair(p).fishnums = combos(p,:); % Save the output combo
            
            % DO BOTH FISH APPEAR DURING THIS EPOCH AND OVERLAP IN TIME?
            
              firstfishtim =  find(data(j).fish(combos(p,1)).freq(:,1) > startim & data(j).fish(combos(p,1)).freq(:,1) < endtim); 
            if ~isempty(firstfishtim)  
              secondfishtim = find(data(j).fish(combos(p,2)).freq(:,1) > startim & data(j).fish(combos(p,2)).freq(:,1) < endtim); 
            if ~isempty(secondfishtim)
                
            clear FD Descartes dF ;                

            sharedtims = intersect(data(j).fish(combos(p,1)).freq(firstfishtim,1), data(j).fish(combos(p,2)).freq(secondfishtim,1));
            
            if ~isempty(sharedtims) % FISH MIGHT BE INTERACTING!!!!!!
                
                if length(sharedtims) > 50 % MEETS MINIMUM NUMBER OF SAMPLES (THIS NEEDS TO BE EDITABLE)
                
                for yy = length(sharedtims):-1:1
                    
                    % What is the dF?                    
                        dF(yy) = abs(data(j).fish(combos(p,1)).freq(data(j).fish(combos(p,1)).freq(:,1) == sharedtims(yy),2) - ...
                            data(j).fish(combos(p,2)).freq(data(j).fish(combos(p,2)).freq(:,1) == sharedtims(yy),2) );                    
                    
                    % How far apart are they?
                        XY(1,1) = data(j).fish(combos(p,1)).x(data(j).fish(combos(p,1)).freq(:,1) == sharedtims(yy));
                        XY(2,1) = data(j).fish(combos(p,1)).y(data(j).fish(combos(p,1)).freq(:,1) == sharedtims(yy));
                        XY(1,2) = data(j).fish(combos(p,2)).x(data(j).fish(combos(p,2)).freq(:,1) == sharedtims(yy));
                        XY(2,2) = data(j).fish(combos(p,2)).y(data(j).fish(combos(p,2)).freq(:,1) == sharedtims(yy));
                        
                        Descartes(yy) = pdist(XY);
                                        
                end
                                
                out(j).pair(p).dFmean = mean(dF);
                out(j).pair(p).dFvar = var(dF);
                out(j).pair(p).meanDist = mean(Descartes);
                out(j).pair(p).varDist = var(Descartes);
                
                out(j).pair(p).sharedtims = sharedtims;
                out(j).pair(p).descartes = Descartes;
                out(j).pair(p).dF = dF;
                
                        % maxlen = min([length(dF), length(Descartes)]);
                        if abs(length(dF) - length(Descartes)) > 0; fprintf('Yowza!'); end
                FD(:,1) = dF;
                FD(:,2) = Descartes;
                [out(j).pair(p).covDistdF, out(j).pair(p).covDistdFpval] = corrcoef(FD);
                
                
                % Calculate Correlations
                
                
                
                
                
                
                
                
                
                end
            end % Fish might be interacting
           
           
            end % Did we share time in this epoch?
            end % Did we share time in this epoch?
            
        end % For this pair of fish
        
        
        
    end % If we have more than one fish
    
    
    
    
    
    
end % For every entry