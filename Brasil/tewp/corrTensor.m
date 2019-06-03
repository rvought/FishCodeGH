function out = corrTensor(dFdist, orig, in, rango)
% This is the tensor plot based on dF correlations
figure(1); clf; hold on;

thresh = 0.5;

if nargin > 3
    startim = rango(1);
    endtim = rango(2);
end
if nargin == 3
    startim = 0;
    endtim = dFdist(1).tim(end);
end

numcolors = 40;
oclrs = lines(numcolors); % 40 colors for plotting up to 40 different fishes
clrs = zeros(numcolors, 3);
% figure(3); clf; hold on; for kk=1:20; plot(kk,kk, '*', 'MarkerEdgeColor', oclrs(kk,:)); end;
%     shuff = randperm(numcolors);
%     for i=1:numcolors; clrs(i,:) = oclrs(shuff(i),:); end


for k=1:length(in(1).Corr)

%% dF vs Distance    
    figure(1); clf; 
    subplot(121); hold on; axis([-100, 200, -200, 100]);
    subplot(122); hold on; axis([-100, 200, -200, 100]);
    
    for j=1:length(dFdist)
    
    aa = dFdist(j).fishnums(1);
    bb = dFdist(j).fishnums(2);

    subplot(121);
    % If the correlation is below threshold, plot a thin grey line
    if abs(in(j).Corr(k)) < thresh   
        plot([orig(aa).xy(k,1), orig(bb).xy(k,1)],  [orig(aa).xy(k,2), orig(bb).xy(k,2)], '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
    end
    
    % If the correlation is above threshold, plot with color and thickness
    if abs(in(j).Corr(k)) > thresh   
        
        if in(j).Corr(k) < 0
            plot([orig(aa).xy(k,1), orig(bb).xy(k,1)],  [orig(aa).xy(k,2), orig(bb).xy(k,2)], 'm-', 'LineWidth', 5*abs(in(j).Corr(k)));
        end
        if in(j).Corr(k) >= 0
            plot([orig(aa).xy(k,1), orig(bb).xy(k,1)],  [orig(aa).xy(k,2), orig(bb).xy(k,2)], 'g-', 'LineWidth', 5*in(j).Corr(k));
        end
        
    end
    
    subplot(122);
    % If the correlation is below threshold, plot a thin grey line
    if abs(in(j).eodCorr(k)) < thresh   
        plot([orig(aa).xy(k,1), orig(bb).xy(k,1)],  [orig(aa).xy(k,2), orig(bb).xy(k,2)], '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
    end
    
    % If the correlation is above threshold, plot with color and thickness
    if abs(in(j).eodCorr(k)) > thresh   
        
        if in(j).eodCorr(k) < 0
            plot([orig(aa).xy(k,1), orig(bb).xy(k,1)],  [orig(aa).xy(k,2), orig(bb).xy(k,2)], 'm-', 'LineWidth', 5*abs(in(j).eodCorr(k)));
        end
        if in(j).eodCorr(k) >= 0
            plot([orig(aa).xy(k,1), orig(bb).xy(k,1)],  [orig(aa).xy(k,2), orig(bb).xy(k,2)], 'g-', 'LineWidth', 5*in(j).eodCorr(k));
        end
        
    end
    
    
    end
        
    % Plot fish position dots
    subplot(121)
    for j=1:length(orig)
        plot(orig(j).xy(k,1), orig(j).xy(k,2), '.', 'MarkerSize', 32, 'Color', oclrs(j,:));
    end
    
    subplot(122);
    for j=1:length(orig)
        plot(orig(j).xy(k,1), orig(j).xy(k,2), '.', 'MarkerSize', 32, 'Color', oclrs(j,:));
    end
    
    
    pause(0.01);
    
end

out = 1;
% 
%         % Get the mean location for each fish 
%             for j=1:numfish
%                fishX(j) =  mean(dat(j).tx(dat(j).tim > startim & dat(j).tim < endtim));
%                fishY(j) =  mean(dat(j).ty(dat(j).tim > startim & dat(j).tim < endtim));
%             end
%     
%         combos = combnk(1:numfish, 2); % All pairwise combinations of fish
%         
%         for p = 1:length(combos) % For each pair of fish
%             
%             out.pair(p).fishnums = combos(p,:); % Save the output combo
%             
%             % DO BOTH FISH APPEAR DURING THIS EPOCH AND OVERLAP IN TIME?
%                 % Get first fish tims
%                 firstfishidx =  find(dat(combos(p,1)).tim > startim & dat(combos(p,1)).tim < endtim); 
%             if ~isempty(firstfishidx)
%                 % Get second fish tims
%                 secondfishidx = find(dat(combos(p,2)).tim > startim & dat(combos(p,2)).tim < endtim); 
%             if ~isempty(secondfishidx)
%                 
%             clear FD sharedtims Descartes dF ;                
% 
%             sharedtims = intersect(dat(combos(p,1)).tim(firstfishidx), dat(combos(p,2)).tim(secondfishidx));
%             
%             if ~isempty(sharedtims) % FISH MIGHT BE INTERACTING!!!!!! (That is, they share tims in this epoch
%                 
%                 if length(sharedtims) > minsharedtim % Set a minimum limit for the amount of time interacting.
%                 
%                 for yy = 1:length(sharedtims)
%                     
%                     % What is the dF?                    
%                         dF(yy) = abs(dat(combos(p,1)).freq(dat(combos(p,1)).tim == sharedtims(yy)) - ...
%                             dat(combos(p,2)).freq(dat(combos(p,2)).tim == sharedtims(yy)) );                    
%                     
%                     % How far apart are they?
%                         XY(1,1) = dat(combos(p,1)).tx(dat(combos(p,1)).tim == sharedtims(yy));
%                         XY(2,1) = dat(combos(p,1)).ty(dat(combos(p,1)).tim == sharedtims(yy));
%                         XY(1,2) = dat(combos(p,2)).tx(dat(combos(p,2)).tim == sharedtims(yy));
%                         XY(2,2) = dat(combos(p,2)).ty(dat(combos(p,2)).tim == sharedtims(yy));
%                         Descartes(yy) = pdist(XY);
%                                         
%                 end
%                                 
%                 out.pair(p).dFmean = mean(dF);
%                 out.pair(p).dFvar = var(dF);
%                 out.pair(p).meanDist = mean(Descartes);
%                 out.pair(p).varDist = var(Descartes);
%                 
%                 out.pair(p).sharedtims = sharedtims;
%                 out.pair(p).descartes = Descartes;
%                 out.pair(p).dF = dF;
%                 
%                         % maxlen = min([length(dF), length(Descartes)]);
%                         if abs(length(dF) - length(Descartes)) > 0; fprintf('Yowza!'); end
%                 FD(:,1) = dF;
%                 FD(:,2) = Descartes;
%                 [out.pair(p).covDistdF, out.pair(p).covDistdFpval] = corrcoef(FD);
%                 
%                 % TENSOR PLOT!
%                                 
%                 if out.pair(p).covDistdFpval(1,2) < pvalthresh % If it is a significant interaction, plot it!
%                     if out.pair(p).covDistdF(1,2) < 0 % Negative correlation
% plot([fishX(combos(p,1)), fishX(combos(p,2))], [fishY(combos(p,1)), fishY(combos(p,2))], 'm-', 'LineWidth', 0.1*abs(floor(log10(out.pair(p).covDistdFpval(1,2)))));
%                     end               
%                     if out.pair(p).covDistdF(1,2) > 0 % Positive correlation
% plot([fishX(combos(p,1)), fishX(combos(p,2))], [fishY(combos(p,1)), fishY(combos(p,2))], 'g-', 'LineWidth', 0.1*abs(floor(log10(out.pair(p).covDistdFpval(1,2)))));
%                     end               
%                 else
%                     
% plot([fishX(combos(p,1)), fishX(combos(p,2))], [fishY(combos(p,1)), fishY(combos(p,2))], '-', 'LineWidth', 1, 'Color', [0.6 0.6 0.6]);
% 
%                 end                
%                 
%                 end % Fish interacted for sufficient time
%             end % Fish might be interacting
%            
%            
%             end % Did we share time in this epoch?
%             end % Did we share time in this epoch?
%             
%             
%             
%         end % For this pair of fish
%         
%         % At the last minute, plot the positions of each of the fish.
%             for j=1:numfish
%                plot(fishX(j), fishY(j), '.k', 'MarkerSize', 40) %, 'Color', clrs(j,:)); 
%                axis off; set(gcf,'color','w');
%                % axis([-50, 200, -100, 75]); % For f3
%                axis([-150, 100, -50, 100]);
%             end
%         
%         
%     end % If we have more than one fish
