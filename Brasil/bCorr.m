function out = bCorr(data, rango)

startim = rango(1); endtim = rango(2);

if length(data) == 1
    figure(1); clf; 
    subplot(211); hold on; xlabel('Time, s'); ylabel('Frequency, Hz');
    subplot(223); hold on; xlabel('cm'); ylabel('cm'); 
    subplot(224); hold on; xlabel('Samples'); ylabel('Velocity'); 
end

figure(2); clf;
subplot(121); hold on; xlabel('Mean EOD Frequency, Hz'); ylabel('EOD variance, Hz');
subplot(122); hold on; xlabel('Mean EOD Frequency, Hz'); ylabel('Mean velocity');

figure(3); clf;
subplot(311); hold on; xlabel('Samples'); ylabel('Correlation R Dist dF');
subplot(312); hold on; xlabel('Samples'); ylabel('Mutual Information Dist dF');
subplot(313); hold on; xlabel('Samples'); ylabel('Cross Correlation peak Dist dF');

figure(4); clf;
subplot(311); hold on; xlabel('Samples'); ylabel('Correlation R EODs');
subplot(312); hold on; xlabel('Samples'); ylabel('Mutual Information EODs');
subplot(313); hold on; xlabel('Samples'); ylabel('Cross Correlation peak EODs');


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
        
        % Some basic (and probably not very useful) measurements
        
            out(j).fish(ff).mVel = mean(vel);
            out(j).fish(ff).mfiltVel = mean(medfilt1(vel,5));
            out(j).fish(ff).totalDist = sum(dist);
            out(j).fish(ff).totalfiltDist = sum(medfilt1(dist,5));

        % If we only did one location, make a plot of the raw data    
        if length(data) == 1
            figure(1);            
            subplot(211); plot(data(j).fish(ff).freq(cts,1), data(j).fish(ff).freq(cts,2), 'LineWidth', 2)
            subplot(223); plot(data(j).fish(ff).x(cts), data(j).fish(ff).y(cts), '*');
            subplot(224); plot(medfilt1(vel,5), 'LineWidth', 1.5);
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
                
                
                % Calculate Correlations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                stepsize = 5; % How many seconds to move forward
                analtime = 100; % Window for correlation analysis in seconds
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
                        
                    % PEARSON CORRELATION COEFICIENT, dF versus distance    
                    [r, pVal] = corrcoef(curdistrack(curridx), curdFs(curridx)); 
                       out(j).corr(p).ddr(kk) = r(2);
                       out(j).corr(p).ddp(kk) = pVal(2);
                    
                    % MUTUAL INFORMATION, dF versus distance
                    out(j).corr(p).ddMIs(kk) = mi(curdistrack(curridx), curdFs(curridx)) / analtime;
                    
                    % CROSS CORRELATION, dF versus distance
                    aa = xcorr(curdistrack(curridx), curdFs(curridx));
                      [~, idx] = max(abs(aa));
                      out(j).corr(p).ddxcorr(kk) = aa(idx);
                      out(j).corr(p).ddxcorrtime(kk) = 2*idx/length(aa);
                    
                    % PEARSON CORRELATION COEFICIENT, EOD vs EOD
                    
                    firstEOD = data(j).fish(combos(p,1)).freq(curridx,2);
                    secondEOD = data(j).fish(combos(p,2)).freq(curridx,2);
                    
                    [r, pVal] = corrcoef(firstEOD - mean(firstEOD), secondEOD - mean(secondEOD) );
                        out(j).corr(p).ffr(kk) = r(2);
                        out(j).corr(p).ffp(kk) = pVal(2);

                    % MUTUAL INFORMATION, EOD versus EOD
                    out(j).corr(p).ffMIs(kk) = mi(data(j).fish(combos(p,1)).freq(curridx,2), data(j).fish(combos(p,2)).freq(curridx,2)) / analtime;
                    
                    % CROSS CORRELATION, EOD vs EOD
                    aa = xcorr(data(j).fish(combos(p,1)).freq(curridx,2) - mean(combos(p,1)).freq(curridx,2), ...
                        data(j).fish(combos(p,2)).freq(curridx,2) - mean(data(j).fish(combos(p,2)).freq(curridx,2)) );
                      [~, idx] = max(abs(aa));
                      out(j).corr(p).ffxcorr(kk) = aa(idx);
                      out(j).corr(p).ffxcorrtime(kk) = 2*idx/length(aa);

                    end
                    
                    startimothy = startimothy + stepsize;
                    
                end
                                
                end % Had enough samples to be worth analysis
            end % Fish might be interacting
           
           
            end % Did we share time in this epoch?
            end % Did we share time in this epoch?
            
            figure(3); subplot(311); plot(out(j).corr(p).ddr, '*-')
            figure(3); subplot(312); plot(out(j).corr(p).ddMIs, '*-')
            figure(3); subplot(313); plot(abs(out(j).corr(p).ddxcorr), '*-')
            
            figure(4); subplot(311); plot(out(j).corr(p).ffr, '*-')
            figure(4); subplot(312); plot(out(j).corr(p).ffMIs, '*-')
            figure(4); subplot(313); plot(abs(out(j).corr(p).ffxcorr), '*-')
            
        end % For this pair of fish
                
        
    end % If we have more than one fish
    
    
    
    
    
    
end % For every entry