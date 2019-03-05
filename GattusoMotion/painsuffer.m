function stims = painsuffer(st, intra, volt, tub, pos, keyb)


fprintf('There are %i events in this file.\n', length(keyb.times));

smp = input('Which events do you want to examine? e.g. [4 5 6 7]');

maxdur = 120;

curdx = 1;

% IF STATEMENT TO MANAGE VOLTAGES
if isempty(volt)
    voltidx = input('What is the size index of this object?');
    volttxt = input('What is the size of this object?');
    sizeIDX = voltidx;
    sizeTXT = volttxt;
    voltint = 1;
    vlength = 1;
else
    voltint = volt.interval;
    vlength = volt.length;
    
end

% IF STATEMENT TO MANAGE INTRA
if isempty(intra)
    intra=pos;
    
else
    intra=intra;
end
    
    
% The various sample frequencies
vFs = 1/voltint;
tFs = 1/tub.interval;
pFs = 1/pos.interval;
iFs = 1/intra.interval;

% Make time sequences
vtim = 1/vFs:1/vFs:vlength/vFs;
ttim = 1/tFs:1/tFs:tub.length/tFs;
ptim = 1/pFs:1/pFs:pos.length/pFs;
itim = 1/iFs:1/iFs:intra.length/iFs;

    

%% This assigns the indexes for the stimulus identity (idx(x,1)) and stimulus size (idx(x,2))

for j = smp % For every stimulus that our user asked to examine
        if (isempty(volt) == 0)               
        curV = volt.values(find(vtim > keyb.times(j), 1)); % Voltage at start of stimulus
        
        if curV > 1.10 && curV < 1.40; sizeIDX = 1; sizeTXT = 'S1'; end % S1
        if curV > 0.80 && curV < 1.00; sizeIDX = 2; sizeTXT = 'S2'; end % S2
        if curV > 0.60 && curV < 0.75; sizeIDX = 3; sizeTXT = 'S3'; end % S3
        if curV > 0.30 && curV < 0.60; sizeIDX = 4; sizeTXT = 'S4'; end % S4
        if curV > 4.95; sizeIDX = 5; sizeTXT = 'M1'; end   %              % M1
        if curV > 4.70 && curV < 4.95; sizeIDX = 6; sizeTXT = 'M2'; end % M2
        if curV > 3.40 && curV < 3.60; sizeIDX = 7; sizeTXT = 'M3'; end % M3
        if curV > 2.40 && curV < 2.80; sizeIDX = 8; sizeTXT = 'L'; end % L  
        end

        
        fprintf('This is stim %i with size %s. \n', keyb.codes(j), sizeTXT);
        
        startim = keyb.times(j);
        endtim = keyb.times(j)+maxdur;
        if j ~= length(keyb.times)
            endtim = min([keyb.times(j)+maxdur keyb.times(j+1)]);
        end
        
        % Plot the stimulus trace with spikes
        close(figure(1)); figure(1); clf; 
    
            %plottuberous
        tt = find(ttim > startim & ttim < endtim);    
            ax(3) = subplot(313); plot(ttim(tt), tub.values(tt), 'g', 'LineWidth', 2);
            hold on
        ii = find(itim > startim & itim < endtim);
            ax(2) = subplot(311); plot(itim(ii), intra.values(ii), 'b', 'LineWidth', 2);
            hold on;
            maxY = ax(2).YLim(2) - 0.20*(ax(2).YLim(2) - ax(2).YLim(1)) ;
            spikz = st.times(st.times > startim & st.times < endtim);
            Ys = ones(length(spikz))*maxY;
            plot(spikz, Ys, 'r*');
        pp = find(ptim > startim & ptim < endtim);
        ax(1) = subplot(312); plot(ptim(pp), pos.values(pp), 'k', 'LineWidth', 2);
            hold on;

            
          
            
            linkaxes(ax, 'x');
            
            
            numstims = input('How many epochs of this stimulus?  ');
            
            for k = 1:numstims
               
                fprintf('Click on start and end of stimulus.\n');
%                 [Xs, ~] = ginput(2);
                a = impoint;
                b = impoint; 
                Pstart = getPosition(a);
                Pend = getPosition(b);
                Xs(1) = Pstart(:,1);
                Xs(2) = Pend(:,1);
                
               
                
                subplot(312)
                hold on
                plot([Xs(1), Xs(1)], [-1 1], 'g');
                hold on
                plot([Xs(2), Xs(2)], [-1 1], 'r');
               
                % Extract stimulus
                stims(curdx).pos = pos.values(ptim > Xs(1) & ptim < Xs(2));
                stims(curdx).pFs = pFs;
                % Size of stimulus
                stims(curdx).sizeDX = sizeIDX; stims(curdx).size = sizeTXT;
                % Extract spikes
                stims(curdx).st = spikz(spikz > Xs(1) & spikz < Xs(2)) - Xs(1);
                
                % Extract frequencies of the stimulus from the user!
                freqs = input('List frequencies or 99 for CV stimulus: \n');
                
                if freqs == 99; cv = mode(abs(diff(stims(curdx).pos))); freqs = 0; end
                if freqs ~= 99; cv = 0; end
                
                stims(curdx).freqs = freqs; stims(curdx).cv = cv;
                
                stims(curdx).dur = Xs(2) - Xs(1);
                
                                
                curdx = curdx + 1;
               
                
                
                
            end
            
            
            stims(curdx).tub = input('Hit 1 for constant tuberous and 2 for oscillating.');
            stims(curdx).date = input('Input date');
            stims(curdx).filename = input('File Number');
            stims(curdx).USID = input('USID');

            
        
end



%% Take the data to make our spiketimes structure[typeNs
% What sizes do we have?

% sizeIDX = unique(typeNsize(:,2)); 
% stimIDX = unique(typeNsize(:,1)); 
% 
% 
% for j = 1:length(sizeIDX) % Loop for each size
%     
%     for k = 1:length(stimIDX) % Loop for each stimulus 
%               
%        ouridxs = find(typeNsize(:,2) == sizeIDX(j) & typeNsize(:,1) == stimIDX(k));
% %  %  length(ouridxs) 
% %        for p = 1:length(ouridxs)
% %            sz(j).spiketimes{p} = st.times(st.times > keyb.times(ouridxs(p)) & st.times < keyb.times(ouridxs(p)) + dur(stimIDX(k)-47)) - keyb.times(ouridxs(p));
% %        end
%         sz(j).spiketimes = st.times(st.times > keyb.times(ouridxs(p)) & st.times < keyb.times(ouridxs(p)) + dur(stimIDX(k)-47)) - keyb.times(ouridxs(p));
%         
%     end
% end


% for kk = 1:length(keyb.times);
%     
%     stimstartim = keyb.times(kk);
%     
%     
%     sz(kk).pos = pos.values(find(ptim > stimstartim & ptim < stimstartim + dur(keyb.codes(kk)-47)));
%     sz(kk).ptim = ptim(ptim > stimstartim & ptim < stimstartim+dur(keyb.codes(kk)-47)) - stimstartim;
%     sz(kk).tub = tub.values(ttim > stimstartim & ttim < stimstartim+dur(keyb.codes(kk)-47));
%     sz(kk).ttim = ttim(ttim > stimstartim & ttim < stimstartim+dur(keyb.codes(kk)-47)) - stimstartim;
%     sz(kk).vFs=vFs;
%     sz(kk).tFs=tFs;
%     sz(kk).pFs=pFs;
%     sz(kk).dur = dur(keyb.codes(kk)-47);
%     sz(kk).spiketimes = st.times(st.times > keyb.times(kk) & st.times < keyb.times(kk) + dur(keyb.codes(kk)-47)) - keyb.times(kk);
% end


%% USEFUL CODE
% 
% % PICK A SIZE/STIM COMBO
% sizeIwant = 1; % Pick a size (1-8)
% stimIwant = 48; % Pick a stimulus (48-57)
% 
% figure(1); clf;
%     ax(1) = subplot(212); plot(dat(sizeIwant).stim(stimIwant).ptim, dat(sizeIwant).stim(stimIwant).pos);
% 
%     ax(2) = subplot(211); hold on;
% for pp = 1:length(dat(sizeIwant).stim(stimIwant).spiketimes) % Usually once, but sometimes more repetitions of same stimulus
%     for k=1:length(dat(1).stim(48).spiketimes{1}); 
%         plot(dat(1).stim(48).spiketimes{1}(k), pp-1, 'r*'); 
%     end;
% end;
% 
% 





