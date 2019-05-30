% figure(22); clf; hold on; clear TE;
% for k = 1000:1020
%     for j=1:30 
%         [TE(j),~,~] = transferEntropyPartition(data.trial{k}(1,:), data.trial{k}(2,:), 1, j*5); 
%     end
% plot(TE)
% end
% 
% %% Plot raw signals
% figure(23); clf; 
% 
% for k = 1000:1:1005
% 
% dF = data.trial{k}(2,:) - min(data.trial{k}(2,:));
% 
% subplot(211); hold on; plot(data.trial{k}(1,:)/max(data.trial{k}(1,:)), 'b'); 
% subplot(212); hold on; plot(dF/max(dF), 'm')
% 
% end

f = cave(3).fish(1).freq(:,2);
ff = cave(3).fish(3).freq(:,2);

for j=length(cave(3).fish(1).freq(:,1)):-1:1
   
    if ~isnan(cave(3).fish(1).freq(j,2)) && ~isnan(cave(3).fish(5).freq(j,2))
        % Calculate distance
        dist(j) = sqrt((cave(3).fish(1).x(j) - cave(3).fish(3).x(j)).^2 + (cave(3).fish(1).y(j) - cave(3).fish(3).y(j)).^2);
        % Calculate dF
        dF(j) = abs(cave(3).fish(1).freq(j,2) - cave(3).fish(3).freq(j,2));
    end
    
end

figure(1); clf; subplot(311); plot(f); hold on; plot(ff);
figure(1); subplot(312); plot(dist);
figure(1); subplot(313); plot(dF);



%%

dist(:,k) = sqrt((dat.fish(C(k,1)).x - dat.fish(C(k,2)).x).^2 + (dat.fish(C(k,1)).y - dat.fish(C(k,2)).y).^2);
df(:,k) = abs(dat.fish(C(k,1)).freq(:,2) - dat.fish(C(k,2)).freq(:,2));
