function fish = bByFish(in)

allfish = [in.pair.fishnums];

maxfish = max(allfish);

for j=1:maxfish
   
    idx = find(allfish == j);
    
    corsum(j).ffr = zeros(1,length(in.corr(1).ffr));
%    corsum(j).ddr = zeros(1,length(in.corr(1).ddr));
        
    for k=1:length(idx)
        
        if idx(k) == 1 
            curdx = 1; 
        else
            if mod(idx(k),2) == 1
                curdx = (idx(k)-1)/2;
            else
                curdx = idx(k)/2;
            end
        end
       
        figure(j); 
            subplot(311); hold on; plot(in.corr(curdx).ffr)
            subplot(312); hold on; plot(in.corr(curdx).ffMIs)
            subplot(313); hold on; plot(in.corr(curdx).ffxcorr)
%             subplot(311); hold on; plot(in.corr(curdx).ddr)
%             subplot(312); hold on; plot(in.corr(curdx).ddMIs)
%             subplot(313); hold on; plot(in.corr(curdx).ddxcorr)
            
            corsum(j).ffr = corsum(j).ffr + abs(in.corr(curdx).ddr);
%            corsum(j).ddr = corsum(j).ddr + abs(in.corr(curdx).ddr);
        
    end
    
    figure(j); subplot(311); plot(corsum(j).ffr/max(corsum(j).ffr), 'k-', 'LineWidth', 2);
%     figure(j); subplot(311); plot(corsum(j).ddr/max(corsum(j).ddr), 'k-', 'LineWidth', 2);
    
end

fish = 1;