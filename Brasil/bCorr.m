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

figure(3); clf;
subplot(311); hold on; xlabel('Samples'); ylabel('Correlation R');
subplot(312); hold on; xlabel('Samples'); ylabel('Mutual Information');
subplot(313); hold on; xlabel('Samples'); ylabel('Cross Correlation peak');


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
                
            sharedidx = []; 
            FD = []; 
            Descartes =[];
            dF = []; 

            sharedtims = intersect(data(j).fish(combos(p,1)).freq(firstfishtim,1), data(j).fish(combos(p,2)).freq(secondfishtim,1));
            
            if ~isempty(sharedtims) % FISH MIGHT BE INTERACTING!!!!!!
                
                if length(sharedtims) > 50 % MEETS MINIMUM NUMBER OF SAMPLES (THIS NEEDS TO BE EDITABLE)
                length(sharedtims)
                for yy = length(sharedtims):-1:1
                                        
                    sharedidx(yy) = find(data(j).fish(combos(p,1)).freq(:,1) == sharedtims(yy));
                    
                    % What is the dF?                    
                        dF(yy) = abs(data(j).fish(combos(p,1)).freq(sharedidx(yy),2) - ...
                            data(j).fish(combos(p,2)).freq(sharedidx(yy),2) );                    
                    
                    % How far apart are they?
                        XY(1,1) = data(j).fish(combos(p,1)).x(sharedidx(yy));
                        XY(2,1) = data(j).fish(combos(p,1)).y(sharedidx(yy));
                        XY(1,2) = data(j).fish(combos(p,2)).x(sharedidx(yy));
                        XY(2,2) = data(j).fish(combos(p,2)).y(sharedidx(yy));
                        
                        Descartes(yy) = pdist(XY);
                                        
                end
                                                
                out(j).pair(p).dFmean = mean(dF);
                out(j).pair(p).dFvar = var(dF);
                out(j).pair(p).meanDist = mean(Descartes);
                out(j).pair(p).varDist = var(Descartes);
                
                out(j).pair(p).sharedtims = sharedtims;
                out(j).pair(p).descartes = Descartes;
                out(j).pair(p).dF = dF;
                
                        if abs(length(dF) - length(Descartes)) > 0; fprintf('Yowza!'); end
                        
                FD(:,1) = dF;
                FD(:,2) = Descartes;
                [out(j).pair(p).covDistdF, out(j).pair(p).covDistdFpval] = corrcoef(FD);
                
                
                % Calculate Correlations
                
                stepsize = 10; % How many seconds to move forward
                analtime = 20; % Window for correlation analysis in seconds
                stepz = (max(sharedtims) - analtime) / stepsize;
                out(j).corr(p).Fs = 1/stepsize; 
                startimothy = min(sharedtims);
                
                % Get rid of NaNs from the data (fillmissing linear) and
                % subtract the means for clean crosscorrelation analyses
                
                curdistrack = fillmissing(Descartes, 'linear') - mean(fillmissing(Descartes, 'linear'));
                curdFs = fillmissing(dF, 'linear') - mean(fillmissing(dF, 'linear'));                
               
                for kk = 1:stepz

                    curridx = sharedidx(sharedtims > startimothy & sharedtims < startimothy+analtime);
 
                    if length(curridx) > 25 % Half of the cutoff above
                    [r, pVal] = corrcoef(curdistrack(curridx), curdFs(curridx)); 
                       out(j).corr(p).r(kk) = r(2);
                       out(j).corr(p).p(kk) = pVal(2);
                       
                    out(j).corr(p).MIs(kk) = mi(curdistrack(curridx), curdFs(curridx)) / analtime;
                    
                    aa = xcorr(curdistrack(curridx), curdFs(curridx));
                    
                    [~, idx] = max(abs(aa));
                    out(j).corr(p).xcorr(kk) = aa(idx);
                    out(j).corr(p).xcorrtime(kk) = 2*idx/length(aa);
                    end
                    
                    startimothy = startimothy + stepsize;
                    
                end
                                
                end % Had enough samples to be worth analysis
            end % Fish might be interacting
           
           
            end % Did we share time in this epoch?
            end % Did we share time in this epoch?
            
            figure(3); subplot(311); plot(out(j).corr(p).r, '*-')
            figure(3); subplot(312); plot(out(j).corr(p).MIs, '*-')
            figure(3); subplot(313); plot(abs(out(j).corr(p).r), '*-')
            
        end % For this pair of fish
                
        
    end % If we have more than one fish
    
    
    
    
    
    
end % For every entry