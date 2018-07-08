function [typeNsize, sz] = voltmeter(st, volt, tub, pos, keyb)


vFs = 1/volt.interval;
tFs = 1/tub.interval;
pFs = 1/pos.interval;

vtim = 1/vFs:1/vFs:volt.length/vFs;
ttim = 1/tFs:1/tFs:tub.length/tFs;
ptim = 1/pFs:1/pFs:pos.length/pFs;

max(vtim)
vtim(end)

%% This assigns the indexes for the stimulus identity (idx(x,1)) and stimulus size (idx(x,2))
for j = 1:length(keyb.times)
    
        typeNsize(j,1) = keyb.codes(j); % This is the stimulus code % WILL WANT TO REASSIGN
   
%         if keyb.codes(j)==48; idx(j,1) = 1; end; %code 48
%         if keyb.codes(j)==49; idx(j,1) = 2; end; %code 49
%         if keyb.codes(j)==50; idx(j,1) = 3; end; %code 50
%         if keyb.codes(j)==51; idx(j,1) = 4; end; %code 51
%         if keyb.codes(j)==52; idx(j,1) = 5; end; %code 52
%         if keyb.codes(j)==53; idx(j,1) = 6; end; %code 53
%         if keyb.codes(j)==54; idx(j,1) = 7; end; %code 54
%         if keyb.codes(j)==55; idx(j,1) = 8; end; %code 55
%         if keyb.codes(j)==56; idx(j,1) = 9; end; %code 56
%         if keyb.codes(j)==57; idx(j,1) = 10; end; %code 57
        
        curV = volt.values(find(vtim > keyb.times(j), 1)); % Voltage at start of stimulus
        
        if curV > 1.10 && curV < 1.40; typeNsize(j,2) = 1; end; % S1
        if curV > 0.80 && curV < 1.00; typeNsize(j,2) = 2; end; % S2
        if curV > 0.60 && curV < 0.75; typeNsize(j,2) = 3; end; % S3
        if curV > 0.30 && curV < 0.60; typeNsize(j,2) = 4; end; % S4
        if curV > 4.95; typeNsize(j,2) = 5; end; %              % M1
        if curV > 4.70 && curV < 4.95; typeNsize(j,2) = 6; end; % M2
        if curV > 3.40 && curV < 3.60; typeNsize(j,2) = 7; end; % M3
        if curV > 2.40 && curV < 2.80; typeNsize(j,2) = 8; end; % L
        
end

dur = [45 45 26 26 24 24 24 24 24 40 40];



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


for kk = 1:length(keyb.times);
    
    stimstartim = keyb.times(kk);
    
    
    sz(kk).pos = pos.values(find(ptim > stimstartim & ptim < stimstartim + dur(keyb.codes(kk)-47)));
    sz(kk).ptim = ptim(ptim > stimstartim & ptim < stimstartim+dur(keyb.codes(kk)-47)) - stimstartim;
    sz(kk).tub = tub.values(ttim > stimstartim & ttim < stimstartim+dur(keyb.codes(kk)-47));
    sz(kk).ttim = ttim(ttim > stimstartim & ttim < stimstartim+dur(keyb.codes(kk)-47)) - stimstartim;
    sz(kk).vFs=vFs;
    sz(kk).tFs=tFs;
    sz(kk).pFs=pFs;
    sz(kk).dur = dur(keyb.codes(kk)-47);
    sz(kk).spiketimes = st.times(st.times > keyb.times(kk) & st.times < keyb.times(kk) + dur(keyb.codes(kk)-47)) - keyb.times(kk);
end


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





