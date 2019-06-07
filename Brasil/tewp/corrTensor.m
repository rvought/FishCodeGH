function out = corrTensor(dFdist, orig, in, makemovie, rango)
% This is the tensor plot based on dF correlations
figure(1); clf; hold on;

thresh = 0.5;
opacitylevel = 0.5;

if makemovie == 1
    writerObj = VideoWriter('delmenowplease.avi');
    writerObj.FrameRate = 15;
    open(writerObj);
end

if nargin > 4
    startim = rango(1);
    endtim = rango(2);
end
if nargin == 4
    startim = 0;
    for j=1:length(dFdist); atmp(j)=length(dFdist(j).dF); end
    endtim = max(atmp);
end

numcolors = 40;
oclrs = lines(numcolors); % 40 colors for plotting up to 40 different fishes
clrs = zeros(numcolors, 3);
% figure(3); clf; hold on; for kk=1:20; plot(kk,kk, '*', 'MarkerEdgeColor', oclrs(kk,:)); end;
%     shuff = randperm(numcolors);
%     for i=1:numcolors; clrs(i,:) = oclrs(shuff(i),:); end

fprintf('Starting Animation');

% Get the total length (many entries will be empty)
for tmp=length(in):-1:1; ourmaxK(tmp) = length(in(tmp).alltim); end

for k = 1:max(ourmaxK)

    fishinthisepoch = [];
    
%% dF vs Distance    
    figure(1); clf; set(gcf,'color','w');
    set(gcf,'Position', [600, 600, 1200, 500]);
    
%     subplot(121); hold on; axis([-100, 200, -150, 150]); % cave 3
%     text(-75, 120, 'dF vs. distance', 'FontSize', 24);
%     subplot(122); hold on; axis([-100, 200, -150, 150]);
%     text(-75, 120, 'EODf vs. EODf', 'FontSize', 24);

%     subplot(121); hold on; axis([-75, 225, -150, 150]); % cave 4,5,6
%     text(-50, 120, 'dF vs. distance', 'FontSize', 24);
%     subplot(122); hold on; axis([-75, 225, -150, 150]);
%     text(-50, 120, 'EODf vs. EODf', 'FontSize', 24);
    
    subplot(121); hold on; axis([-75, 225, -250, 100]); % cave  7
    text(-50, 70, 'dF vs. distance', 'FontSize', 24);
    subplot(122); hold on; axis([-75, 225, -250, 100]);
    text(-50, 70, 'EODf vs. EODf', 'FontSize', 24);

    for j=1:length(in)
        
    aa = dFdist(j).fishnums(1);
    bb = dFdist(j).fishnums(2);
    
    % Do we have any data for these fish?
    if ~isempty(in(j).alltim)
    % Do we have data for both fish?
    if sum(ismember(in(j).tt, in(j).alltim(k))) == 1 
        
        fishinthisepoch = [fishinthisepoch, aa, bb];


        % Get the current XY values for each fish
        posIDX = find(orig(aa).tim >= in(j).alltim(k), 1);
        fishAA = [orig(aa).xy(posIDX,1), orig(aa).xy(posIDX,2)];
        fishBB = [orig(bb).xy(posIDX,1), orig(bb).xy(posIDX,2)];
        
    subplot(121);
    % If the correlation is below threshold, plot a thin grey line
    if abs(in(j).Corr(k)) < thresh   
%        plot([orig(aa).xy(k,1), orig(bb).xy(k,1)],  [orig(aa).xy(k,2), orig(bb).xy(k,2)], '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
        plot([fishAA(1), fishBB(1)],  [fishAA(2), fishBB(2)], '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
    end
    
    % If the correlation is above threshold, plot with color and thickness
    if abs(in(j).Corr(k)) > thresh   
        
        if in(j).Corr(k) < 0
%            h = plot([orig(aa).xy(k,1), orig(bb).xy(k,1)],  [orig(aa).xy(k,2), orig(bb).xy(k,2)], 'm-', 'LineWidth', 8*abs(in(j).Corr(k)));
            h = plot([fishAA(1), fishBB(1)],  [fishAA(2), fishBB(2)], 'm-', 'LineWidth', 8*abs(in(j).Corr(k)));
            h.Color(4)=opacitylevel;
        end
        if in(j).Corr(k) >= 0
%            h = plot([orig(aa).xy(k,1), orig(bb).xy(k,1)],  [orig(aa).xy(k,2), orig(bb).xy(k,2)], 'g-', 'LineWidth', 8*in(j).Corr(k));
            h = plot([fishAA(1), fishBB(1)],  [fishAA(2), fishBB(2)], 'g-', 'LineWidth', 8*in(j).Corr(k));
            h.Color(4)=opacitylevel;
        end
        
    end
    box on;
    set(gca,'Yticklabel', [], 'Xticklabel', []);
    
    subplot(122);
    % If the correlation is below threshold, plot a thin grey line
    if abs(in(j).eodCorr(k)) < thresh   
%        plot([orig(aa).xy(k,1), orig(bb).xy(k,1)],  [orig(aa).xy(k,2), orig(bb).xy(k,2)], '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
        plot([fishAA(1), fishBB(1)],  [fishAA(2), fishBB(2)], '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.5);
    end
    
    % If the correlation is above threshold, plot with color and thickness
    if abs(in(j).eodCorr(k)) > thresh   
        
        if in(j).eodCorr(k) < 0
%            h = plot([orig(aa).xy(k,1), orig(bb).xy(k,1)],  [orig(aa).xy(k,2), orig(bb).xy(k,2)], 'm-', 'LineWidth', 8*abs(in(j).eodCorr(k)));
            h = plot([fishAA(1), fishBB(1)],  [fishAA(2), fishBB(2)], 'm-', 'LineWidth', 8*abs(in(j).eodCorr(k)));
            h.Color(4)=opacitylevel;
        end
        if in(j).eodCorr(k) >= 0
%            h = plot([orig(aa).xy(k,1), orig(bb).xy(k,1)],  [orig(aa).xy(k,2), orig(bb).xy(k,2)], 'g-', 'LineWidth', 8*in(j).eodCorr(k));
            h = plot([fishAA(1), fishBB(1)],  [fishAA(2), fishBB(2)], 'g-', 'LineWidth', 8*in(j).eodCorr(k));
            h.Color(4)=opacitylevel;
        end
        
    end
    box on;
    set(gca,'Yticklabel', [], 'Xticklabel', []);    
    
    end % The if
    end % Testing if there is any data
    end % For each dFdist

    
    % Plot fish position dots
    
    for j=1:length(orig)
        % Plot only the fish that have correlations   
        if sum(ismember(fishinthisepoch, j)) > 1
        subplot(121);
        plot(orig(j).xy(posIDX,1), orig(j).xy(posIDX,2), '.', 'MarkerSize', 36, 'Color', oclrs(j,:));
        subplot(122);
        plot(orig(j).xy(posIDX,1), orig(j).xy(posIDX,2), '.', 'MarkerSize', 36, 'Color', oclrs(j,:));    
        end
    end
    
if makemovie == 1
    mframe = getframe(gcf);
    writeVideo(writerObj,mframe);
end
    
    pause(0.01);
        
end




if makemovie == 1
    close(writerObj);
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
