function out = SpaceCorps(in)

ttf2 = find(~isnan(in.fish(2).freq(:,2)));
ttf3 = find(~isnan(in.fish(3).freq(:,2)));

F1(1,:) = in.fish(2).x(ttf2);
F1(2,:) = in.fish(2).y(ttf2);

F2(1,:) = in.fish(3).x(ttf3);
F2(2,:) = in.fish(3).y(ttf3);

ctrs{1} = -150:5:250;
ctrs{2} = -150:5:250;

f1h = hist3(F1', ctrs);
f2h = hist3(F2', ctrs);

figure(1); clf; 
subplot(121); surf(f1h); view(0,90);
subplot(122); surf(f2h); view(0,90);

out = 1;

