function gridplotter(in, tracks, tim)
% Usage: gridpotter(in, tracks, tim)

figure(1); clf; subplot(121); hold on; subplot(122); hold on;

if nargin == 1
    tracknum = length(in);
    tracks = 1:length(in);
end
if nargin > 1
    tracknum = length(tracks);
end
if nargin < 3
    tim = [0 10000];
end


for i=1:tracknum
    
    tt = find(in(tracks(i)).tim > tim(1) & in(tracks(i)).tim < tim(2));
    txy = find(in(tracks(i)).tt > tim(1) & in(tracks(i)).tt < tim(2));
    
    subplot(121); plot(in(tracks(i)).tim(tt), in(tracks(i)).freq(tt), '*');
        text(in(tracks(i)).tim(end), in(tracks(i)).freq(1)+5, num2str(tracks(i)));
    subplot(122); plot(in(tracks(i)).tx(txy), in(tracks(i)).ty(txy), '*');
    
end;

subplot(121); ylim([250, 500]);
subplot(122); axis([-150, 200, -200, 150]);