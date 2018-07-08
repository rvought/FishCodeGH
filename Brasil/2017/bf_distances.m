function out = bf_distances(in)


numFish = length(in);

pairFish = nchoosek(1:numFish, 2);

for j = 1:length(pairFish)
    
    [~, f1, f2] = intersect(in(pairFish(j,1)).tim, in(pairFish(j,2)).tim);
    
    if length(f1) > 10

            A(:,1) = in(pairFish(j,1)).tx(f1);
            A(:,2) = in(pairFish(j,1)).ty(f1);
            B(:,1) = in(pairFish(j,2)).tx(f2);
            B(:,2) = in(pairFish(j,2)).ty(f2);
        
            for kk = 1:length(A)
               X(1,:) = A(kk,:);
               X(2,:) = B(kk,:);
                out(j).dist(kk) = pdist(X);                
            end
        
            % out(j).dist = pdist2(A,B);
            out(j).pair = pairFish(j,:);
            out(j).f1 = f1;
            out(j).f2 = f2;
            out(j).tim = in(pairFish(j,1)).tim(f1);
            clear A B
    end
    
    
end
    
    
%% Plotting

length(pairFish)

if length(pairFish) < 70
for k = 1:length(pairFish)

    figure;
    subplot(121); 
        plot(in(pairFish(k,1)).tx, in(pairFish(k,1)).ty, '-c', 'LineWidth', 3);
        hold on;
        plot(in(pairFish(k,2)).tx, in(pairFish(k,2)).ty, '-g', 'LineWidth', 3);

        plot(in(pairFish(k,1)).tx(out(k).f1), in(pairFish(k,1)).ty(out(k).f1), '-*b');
        plot(in(pairFish(k,2)).tx(out(k).f2), in(pairFish(k,2)).ty(out(k).f2), '-*r');
        
        xlim([-100 250]); ylim([-150 100]); 
    subplot(122);
        plot(in(pairFish(k,1)).tim, in(pairFish(k,1)).freq, '-c', 'LineWidth', 3);
        hold on;
        plot(in(pairFish(k,2)).tim, in(pairFish(k,2)).freq, '-g', 'LineWidth', 3);

        plot(in(pairFish(k,1)).tim(out(k).f1), in(pairFish(k,1)).freq(out(k).f1), '-*b');
        plot(in(pairFish(k,2)).tim(out(k).f2), in(pairFish(k,2)).freq(out(k).f2), '-*r');
        ylim([200 600]);
end
end


if length(pairFish) >= 70
for kk = 1:16:length(pairFish)
    figure;
    for jj = 1:2:15
    subplot(4,4,jj); 
        plot(in(pairFish(kk+jj-1,1)).tx, in(pairFish(kk+jj-1,1)).ty, '-*');
        hold on;
        plot(in(pairFish(kk+jj-1,2)).tx, in(pairFish(kk+jj-1,2)).ty, '-*');
        xlim([-200 250]); ylim([-150 100]); 
    subplot(4,4,jj+1);
        plot(in(pairFish(kk+jj-1,1)).tim, in(pairFish(kk+jj-1,1)).freq, '-*');
        hold on;
        plot(in(pairFish(kk+jj-1,2)).tim, in(pairFish(kk+jj-1,2)).freq, '-*');
        ylim([250 600]);
    end
    
end
end
