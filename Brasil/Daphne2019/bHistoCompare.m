function out = bHistoCompare(in)


Xedges = -150:10:250;
Yedges = -200:10:200;


for j=length(in.fish):-1:1
    
    out(j).foo = zeros(length(Xedges)-1,length(Yedges)-1);
    out(j).Xs = Xedges(1:end-1);
    out(j).Ys = Yedges(1:end-1);
    
    % For each X
        for k = 1:length(Xedges)-1
            
            % For each Y
                for m = 1:length(Yedges)-1
                    
                    out(j).foo(k,m) = out(j).foo(k,m) + length(find(in.fish(j).x > Xedges(k) ...
                        & in.fish(j).x < Xedges(k+1) & in.fish(j).y > Yedges(m) ...
                        & in.fish(j).y < Yedges(m+1)));
                    
                end
        end
                      
    
    
end

%% Plot
cmp(:,1) = 0.999*ones(1,128);
cmp(:,2) = (256:-1:129)/256;
cmp(:,3) = 0.999*ones(1,128);

figure(1); clf; set(gcf,'renderer','Painters')
for zz = 1:min([length(out), 9])
    subplot(3,3,zz); surf(out(zz).Xs, out(zz).Ys, out(zz).foo, 'EdgeColor', 'none'); view(0,90);
end
colormap(cmp);

if length(out) > 9
   howmany = ceil(length(out)/9);
   
   for kk = 2:howmany
    figure(kk); clf; set(gcf,'renderer','Painters')
    for zz = 1:min([9, length(out)-9*(kk-1)])
        subplot(3,3,zz); surf(out(9*(kk-1)+zz).Xs, out(zz).Ys, out(9*(kk-1)+zz).foo, 'EdgeColor', 'none'); view(0,90);
    end
    colormap(cmp);
   end
    
end
%linkaxes(ax, 'xy'); 
