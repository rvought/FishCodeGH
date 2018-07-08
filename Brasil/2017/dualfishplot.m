function out = dualfishplot(data, evt, fishnum, tim)
% Usage dualfishplot(data, in)

idx{1} = find(data(fishnum(1)).tt > tim(1) & data(fishnum(1)).tt < tim(end));
numcols(1) = floor(length(idx{1})/2);

idx{2} = find(data(fishnum(2)).tt > tim(1) & data(fishnum(2)).tt < tim(end));
numcols(2) = floor(length(idx{2})/2);

afish = autumn(numcols(1));
bfish = winter(numcols(2));

wwq{1} = 0.2:2.3/numcols(1):2.5;
wwq{2} = 0.2:2.3/numcols(2):2.5;
wwr{1} = 2.5:2.3/numcols(1):4.8;
wwr{2} = 2.5:2.3/numcols(2):4.8;

%% Plot the spatial data in Figure 1
figure(1); clf; 
hold on;

    for j=1:numcols(1);
        if ~isempty(data(fishnum(1)).tx(idx{1}(j))) && ~isempty(data(fishnum(1)).tx(idx{1}(j)+1))
            plot([data(fishnum(1)).tx(idx{1}(j)), data(fishnum(1)).tx(idx{1}(j)+1)], [data(fishnum(1)).ty(idx{1}(j)), data(fishnum(1)).ty(idx{1}(j)+1)], 'Color', afish(j,:), 'LineWidth', wwq{1}(j));
        end
    end
    for j=1:numcols(2);
        if ~isempty(data(fishnum(2)).tx(idx{2}(j))) && ~isempty(data(fishnum(2)).tx(idx{2}(j)+1))
        plot([data(fishnum(2)).tx(idx{2}(j)), data(fishnum(2)).tx(idx{2}(j)+1)], [data(fishnum(2)).ty(idx{2}(j)), data(fishnum(2)).ty(idx{2}(j)+1)], 'Color', bfish(j,:), 'LineWidth', wwq{2}(j));
        end
    end
    
    for j=numcols(1):-1:1;
        if ~isempty(data(fishnum(1)).tx(idx{1}(j+numcols(1)))) && ~isempty(data(fishnum(1)).tx(idx{1}(j+numcols(1)-1)))
        plot([data(fishnum(1)).tx(idx{1}(j+numcols(1))), data(fishnum(1)).tx(idx{1}(j+numcols(1)-1))], [data(fishnum(1)).ty(idx{1}(j+numcols(1))), data(fishnum(1)).ty(idx{1}(j+numcols(1)-1))], 'Color', afish(1+numcols(1)-j,:), 'LineWidth', wwr{1}(j));
        end
    end
    for j=numcols(2):-1:1;
        if ~isempty(data(fishnum(2)).tx(idx{2}(j+numcols(2)))) && ~isempty(data(fishnum(2)).tx(idx{2}(j+numcols(2)-1)))
        plot([data(fishnum(2)).tx(idx{2}(j+numcols(2))), data(fishnum(2)).tx(idx{2}(j+numcols(2)-1))], [data(fishnum(2)).ty(idx{2}(j+numcols(2))), data(fishnum(2)).ty(idx{2}(j+numcols(2)-1))], 'Color', bfish(1+numcols(2)-j,:), 'LineWidth', wwr{2}(j));
        end
    end

    %% Distance between fish plot

[~, f1idx, f2idx] = intersect(data(fishnum(1)).tt(idx{1}), data(fishnum(2)).tt(idx{2}));    
length(f1idx)
    
    for j=1:length(f1idx);
        dist(j) = pdist([data(fishnum(1)).tx(idx{1}(f1idx(j))), data(fishnum(1)).ty(idx{1}(f1idx(j))); data(fishnum(2)).tx(idx{2}(f2idx(j))), data(fishnum(2)).ty(idx{2}(f2idx(j)))]);
    end
    length(dist)
    figure(2); clf; ax(2) = subplot(212); 
    plot(data(fishnum(1)).tt(idx{1}(f1idx)), dist, 'k');
    
clear idx numcols afish bfish wwq wwr
    
%% Plot the frequency data in Figure 2
figure(2); ax(1) = subplot(211);

hold on;

idx{1} = find(data(fishnum(1)).tim > tim(1) & data(fishnum(1)).tim < tim(end));
numcols(1) = floor(length(idx{1})/2);

idx{2} = find(data(fishnum(2)).tim > tim(1) & data(fishnum(2)).tim < tim(end));
numcols(2) = floor(length(idx{2})/2);

afish = autumn(numcols(1));
bfish = winter(numcols(2));

wwq{1} = 0.2:2.3/numcols(1):2.5;
wwq{2} = 0.2:2.3/numcols(2):2.5;
wwr{1} = 2.5:2.3/numcols(1):4.8;
wwr{2} = 2.5:2.3/numcols(2):4.8;


    for j=1:numcols(1);
        if ~isempty(data(fishnum(1)).tim(idx{1}(j))) && ~isempty(data(fishnum(1)).tim(idx{1}(j)+1))
            plot([data(fishnum(1)).tim(idx{1}(j)), data(fishnum(1)).tim(idx{1}(j)+1)], [data(fishnum(1)).freq(idx{1}(j)), data(fishnum(1)).freq(idx{1}(j)+1)], 'Color', afish(j,:), 'LineWidth', 2);
        end
    end
    for j=1:numcols(2);
        if ~isempty(data(fishnum(2)).tim(idx{2}(j))) && ~isempty(data(fishnum(2)).tim(idx{2}(j)+1))
            plot([data(fishnum(2)).tim(idx{2}(j)), data(fishnum(2)).tim(idx{2}(j)+1)], [data(fishnum(2)).freq(idx{2}(j)), data(fishnum(2)).freq(idx{2}(j)+1)], 'Color', bfish(j,:), 'LineWidth', 2);
        end
    end
    
    for j=numcols(1):-1:1;
        if ~isempty(data(fishnum(1)).tim(idx{1}(j+numcols(1)))) && ~isempty(data(fishnum(1)).tim(idx{1}(j+numcols(1)-1)))
        plot([data(fishnum(1)).tim(idx{1}(j+numcols(1))), data(fishnum(1)).tim(idx{1}(j+numcols(1)-1))], [data(fishnum(1)).freq(idx{1}(j+numcols(1))), data(fishnum(1)).freq(idx{1}(j+numcols(1)-1))], 'Color', afish(1+numcols(1)-j,:), 'LineWidth', 2);
        end
    end
    for j=numcols(2):-1:1;
        if ~isempty(data(fishnum(2)).tim(idx{2}(j+numcols(2)))) && ~isempty(data(fishnum(2)).tim(idx{2}(j+numcols(2)-1)))
        plot([data(fishnum(2)).tim(idx{2}(j+numcols(2))), data(fishnum(2)).tim(idx{2}(j+numcols(2)-1))], [data(fishnum(2)).freq(idx{2}(j+numcols(2))), data(fishnum(2)).freq(idx{2}(j+numcols(2)-1))], 'Color', bfish(1+numcols(2)-j,:), 'LineWidth', 2);
        end
    end

ylim([250 550]);


    
 linkaxes(ax, 'x');



%% Add Clicks

[xs, ~] = ginput(2);


    fclickidx(1) = find(data(fishnum(1)).tim < xs(1), 1, 'last');
    fclickidx(2) = find(data(fishnum(1)).tim < xs(2), 1, 'last');
    figure(2); subplot(211); plot(data(fishnum(1)).tim(fclickidx(1)), data(fishnum(1)).freq(fclickidx(1)), 'm>');
    figure(2); subplot(211); plot(data(fishnum(1)).tim(fclickidx(2)), data(fishnum(1)).freq(fclickidx(2)), 'mo');

    fclickidx(1) = find(data(fishnum(2)).tim < xs(1), 1, 'last');
    fclickidx(2) = find(data(fishnum(2)).tim < xs(2), 1, 'last');
    figure(2); subplot(211); plot(data(fishnum(2)).tim(fclickidx(1)), data(fishnum(2)).freq(fclickidx(1)), 'g>');
    figure(2); subplot(211); plot(data(fishnum(2)).tim(fclickidx(2)), data(fishnum(2)).freq(fclickidx(2)), 'go');
    
    
    sclickidx(1) = find(data(fishnum(1)).tt < xs(1), 1, 'last');
    sclickidx(2) = find(data(fishnum(1)).tt < xs(2), 1, 'last');
    figure(1); plot(data(fishnum(1)).tx(sclickidx(1)), data(fishnum(1)).ty(sclickidx(1)), 'm>', 'MarkerSize', 5);
    figure(1); plot(data(fishnum(1)).tx(sclickidx(2)), data(fishnum(1)).ty(sclickidx(2)), 'mo');

    sclickidx(1) = find(data(fishnum(2)).tt < xs(1), 1, 'last');
    sclickidx(2) = find(data(fishnum(2)).tt < xs(2), 1, 'last');
    figure(1); plot(data(fishnum(2)).tx(sclickidx(1)), data(fishnum(2)).ty(sclickidx(1)), 'g>', 'MarkerSize', 5);
    figure(1); plot(data(fishnum(2)).tx(sclickidx(2)), data(fishnum(2)).ty(sclickidx(2)), 'go');
    
    
    
    
    out = 1;
    