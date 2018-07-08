zz = 1;

edges = {-250:5:250, -250:5:250};

endtim = 60;

for j=1:length(f)
    for k = 1:length(f(j).c)
        tt = find(f(j).c(k).tx < endtim);
        if length(tt) > 50
        tmp = zeros(length(tt),2);
        Kidx = convhull(f(j).c(k).tx(tt), f(j).c(k).ty(tt));
        
        cavefish(zz).convhullx = f(j).c(k).tx(Kidx);
        cavefish(zz).convhully = f(j).c(k).ty(Kidx);
        
        for p=length(tmp):-1:1
            tmp(p,1) = f(j).c(k).tx(p); tmp(p,2) = f(j).c(k).ty(p);
        end
       cavefish(zz).histo = hist3(tmp, 'Edges', edges);
       cavefish(zz).site = j;
       zz=zz+1;
        end
    end   
    
    for k = 1:length(f(j).s)
        tt = find(f(j).s(k).tx < endtim);
        if length(tt) > 50
        tmp = zeros(length(tt),2);
        Kidx = convhull(f(j).s(k).tx(tt), f(j).s(k).ty(tt));
        surfish(zz).convhullx = f(j).s(k).x(Kidx);
        surfish(zz).convhully = f(j).s(k).y(Kidx);
        for p=length(tmp):-1:1
            tmp(p,1) = f(j).s(k).x(p); tmp(p,2) = f(j).s(k).y(p);
        end
       surfish(zz).histo = hist3(tmp, 'Edges', edges);
       surfish(zz).site = j;
       zz=zz+1;
        end
    end

    
    
end


% surf(cavefish(10).histo, 'FaceColor', 'interp', 'EdgeColor', 'none'); view(0,90);
% caxis([0 100]);


figure(27); clf;
currfish = [32 33 34 35 36 37 38 39 40 41];
for j=1:length(currfish)
    figure(j); clf;
    surf(cavefish(currfish(j)).histo, 'FaceColor', 'interp', 'EdgeColor', 'none'); 
    view(0,90); caxis([0 50]); 
        hold on; plot((250+cavefish(currfish(j)).convhully)/5, (250+cavefish(currfish(j)).convhullx)/5, 'y');
    figure(27); hold on;
    plot(-cavefish(currfish(j)).convhullx, -cavefish(currfish(j)).convhully);
    
end



figure(87); clf;
currfish = [13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33];
for j=1:length(currfish)
    figure(j+50); clf;
    surf(cavefish(currfish(j)).histo, 'FaceColor', 'interp', 'EdgeColor', 'none'); 
    view(0,90); caxis([0 50]); 
        hold on; plot((250+cavefish(currfish(j)).convhully)/5, (250+cavefish(currfish(j)).convhullx)/5, 'y');
    figure(87); hold on;
    plot(-cavefish(currfish(j)).convhullx, -cavefish(currfish(j)).convhully);
    
end