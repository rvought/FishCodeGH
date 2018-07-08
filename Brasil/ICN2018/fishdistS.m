function out = fishdistS(in, dur)
% out = fishdist(in, dur)
% God help me code!

startim = dur(1); endtim = dur(2);

for j=1:length(in) % For each entry
    
    numfish = length(in(j).s); % How many fish in this recording
    
    % SINGLE FISH DATA
    
    for ff = 1:length(in(j).s) % for each fish
        
        cts = find(in(j).s(ff).tim > startim & in(j).s(ff).tim < endtim); % Current epoch
        
        % Mean and variance of frequency of this fish
        out(j).fish(ff).mfreq = mean(in(j).s(ff).freq(cts));
        out(j).fish(ff).vfreq = var(in(j).s(ff).freq(cts));
        out(j).fish(ff).Nfreq = length(cts);

        % Velocity and distance of the fish
        
        for vv = 2:length(cts)
            X(1,1) = in(j).s(ff).tx(cts(vv-1));
            X(2,1) = in(j).s(ff).ty(cts(vv-1));
            X(1,2) = in(j).s(ff).tx(cts(vv));
            X(2,2) = in(j).s(ff).ty(cts(vv));
            dist(vv) = pdist(X);
            vel(vv) = dist(vv) / (in(j).s(ff).tt(cts(vv)) - in(j).s(ff).tt(cts(vv-1)));
        end
        
        %vel = medfilt1(vel,5);
        
        out(j).fish(ff).mvel = mean(vel);
        out(j).fish(ff).totaldist = sum(dist);

    end
        
    % PAIRWISE INTERACTIONS BETWEEN FISH
    
    if numfish > 1 % If we have more than one fish
    
        combos = combnk(1:numfish, 2); % All pairwise combinations of fish
        
        for p = 1:length(combos) % For each pair of fish
            
            out(j).pair(p).fishnums = combos(p,:); % Save the output combo
            
            % DO BOTH FISH APPEAR DURING THIS EPOCH AND OVERLAP IN TIME?
            
                firstfishtim =  find(in(j).s(combos(p,1)).tim > startim & in(j).s(combos(p,1)).tim < endtim); 
            if ~isempty(firstfishtim)  
                secondfishtim = find(in(j).s(combos(p,2)).tim > startim & in(j).s(combos(p,2)).tim < endtim); 
            if ~isempty(secondfishtim)
                
            clear FD sharedtims Descartes dF ;                

            sharedtims = intersect(in(j).s(combos(p,1)).tim(firstfishtim), in(j).s(combos(p,2)).tim(secondfishtim));
            
            if ~isempty(sharedtims) % FISH MIGHT BE INTERACTING!!!!!!
                if length(sharedtims) > 50

                
                for yy = 1:length(sharedtims)
                    
                    % What is the dF?                    
                        dF(yy) = abs(in(j).s(combos(p,1)).freq(in(j).s(combos(p,1)).tim == sharedtims(yy)) - ...
                            in(j).s(combos(p,2)).freq(in(j).s(combos(p,2)).tim == sharedtims(yy)) );                    
                    
                    % How far apart are they?
                        XY(1,1) = in(j).s(combos(p,1)).tx(in(j).s(combos(p,1)).tim == sharedtims(yy));
                        XY(2,1) = in(j).s(combos(p,1)).ty(in(j).s(combos(p,1)).tim == sharedtims(yy));
                        XY(1,2) = in(j).s(combos(p,2)).tx(in(j).s(combos(p,2)).tim == sharedtims(yy));
                        XY(2,2) = in(j).s(combos(p,2)).ty(in(j).s(combos(p,2)).tim == sharedtims(yy));
                        
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
                
                end % minimum length of interaction
            end % Fish might be interacting
           
           
            end % Did we share time in this epoch?
            end % Did we share time in this epoch?
            
        end % For this pair of fish
        
        
        
    end % If we have more than one fish
    
    
    
    
    
    
end % For every entry
        
figure(1); hold on; 
for k = 18:22
    for j=1:length(out(k).pair)
        plot(out(k).pair(j).meanDist, out(k).pair(j).dFmean, 'm*'); 
    end
end

        
        
        
