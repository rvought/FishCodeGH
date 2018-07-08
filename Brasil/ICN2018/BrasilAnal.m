function [out asdf] = BrasilAnal(in)


nM = 5; % Medfilt value for cleaning trajectories
ep = 60; % Epoch duration in seconds

figure; subplot(121); hold on; subplot(122); hold on;

ctrials = [3]; % 5
        

strials = []; % 21

for j = 1:length(ctrials)
    
    for k = 1:length(in(ctrials(j)).c)
        tt = find(in(ctrials(j)).c(k).tim > ep & in(ctrials(j)).c(k).tim < ep*2);
        x = medfilt1(in(ctrials(j)).c(k).tx(tt), nM);
        y = medfilt1(in(ctrials(j)).c(k).ty(tt), nM);
            asdf(k) = polyshape(x, y);
        %subplot(121); plot(in(ctrials(j)).c(k).tx(tt), in(ctrials(j)).c(k).ty(tt), '*');
        subplot(121); plot(x, y, '-');
        subplot(122); plot(in(ctrials(j)).c(k).tim(tt), in(ctrials(j)).c(k).freq(tt), '*');
        
        khull{k} = convhull(x, y);
        subplot(121); plot(x(khull{k}), y(khull{k}), 'LineWidth', 3);
        
    end
end

% for j = 1:length(strials)
%     
%     for k = 1:length(in(strials(j)).s)
%         tt = find(in(strials(j)).s(k).tim > ep & in(strials(j)).s(k).tim < ep*2);
%         x = medfilt1(in(strials(j)).s(k).tx(tt), nM);
%         y = medfilt1(in(strials(j)).s(k).ty(tt), nM);
%         %subplot(121); plot(in(ctrials(j)).c(k).tx(tt), in(ctrials(j)).c(k).ty(tt), '*');
%         subplot(121); plot(x, y, '-');
%         subplot(122); plot(in(strials(j)).s(k).tim(tt), in(strials(j)).s(k).freq(tt), '*');
%         
%         khull{k} = convhull(x, y);
%         subplot(121); plot(x(khull{k}), y(khull{k}), 'LineWidth', 3);
%         
%     end
%     
% end




out =1;