function out = SpaceCorps(in)


F1(1,:) = in.fish(2).x;
F1(2,:) = in.fish(2).y;

F2(1,:) = in.fish(3).x;
F2(2,:) = in.fish(3).y;

ctrs(:,1) = -250:10:250;
ctrs(:,2) = -250:10:250;

f1h = hist3(F1, ctrs);
f2h = hist3(F2, ctrs);

figure(1); clf; 
subplot(121); surf(f1h); view(0,90);
subplot(122); surf(f1h); view(0,90);

