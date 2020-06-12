
%% Initialize variables

c.fftCh1amp = [];
c.fftCh2amp = [];
c.fftCh1freq = [];
c.fftCh2freq = [];
c.rmsCh1 = [];
c.rmsCh2 = [];
c.lightlevel = [];
c.temper = [];
c.idx = [];

c.timmy = [];

%% Load data

% load 7May20b66.mat

%% Concatonate

c.fftCh1amp = [c.fftCh1amp [out.fftCh1peakamp]];
c.fftCh2amp = [c.fftCh2amp [out.fftCh2peakamp]];
    
c.fftCh1freq = [c.fftCh1freq [out.fftCh1peakfreq]];
c.fftCh2freq = [c.fftCh2freq [out.fftCh2peakfreq]];

maxidx = max(c.idx);
    if ~isempty(maxidx)
        c.idx = [c.idx ([out.idx]+ maxidx)];
    else
        c.idx = [out.idx];
    end

c.rmsCh1 = [c.rmsCh1 [out.rmsCh1]];
c.rmsCh2 = [c.rmsCh2 [out.rmsCh2]];

c.lightlevel = [c.lightlevel [out.light]];
c.temper = [c.temper [out.temp]];

% Extract the times

%prelen = length(c.timmy);

%for j=1:length(out)
    %mon = str2num(out(j).fileinfo(16:17));
    %day = str2num(out(j).fileinfo(19:20));
    %yr = str2num(out(j).fileinfo(22:25));
    %hr = str2num(out(j).fileinfo(27:28));
    %mi = str2num(out(j).fileinfo(30:31));
    %ss = str2num(out(j).fileinfo(33:34));
    %c.timmy(j+prelen) = datenum(yr,mon,day,hr,mi,ss);
%

%% Plot 

figure(4); clf;

ax(1)=subplot(311); 
    plot(c.idx, c.fftCh1amp, 'LineWidth', 4); 
    hold on; plot(c.idx, c.fftCh2amp, 'LineWidth', 4);

ax(2)=subplot(312); 
    plot(c.idx, c.rmsCh1, 'LineWidth', 4); 
    hold on; plot(c.idx, c.rmsCh2, 'LineWidth', 4);

ax(3)=subplot(313); plot(1:length(c.lightlevel), c.lightlevel, 'LineWidth', 4);
linkaxes(ax, 'x');



%% Overlay plot

figure(27); clf; hold on;
thresh = 220; % Select based on light values

preidx = 20;
postidx = 50;

lv = zeros(1, max(c.idx));
lv(c.lightlevel>thresh) = 1;
dls = diff(lv);

darkidxs = find(dls == -1);

meany.fft = zeros(length(darkidxs), preidx+postidx);
meany.rms = zeros(length(darkidxs), preidx+postidx);

for j=1:length(darkidxs)
    
        tt = find(c.idx > darkidxs(j)-preidx & c.idx < darkidxs(j)+postidx);

    subplot(211); hold on;
        plot(c.idx(tt)-darkidxs(j), c.fftCh1amp(tt), 'k*');
    subplot(212); hold on;
        plot(c.idx(tt)-darkidxs(j), c.rmsCh1(tt), 'k*');
        
        meany.fft(j,c.idx(tt)-darkidxs(j)+preidx) = c.fftCh1amp(tt);        
        meany.rms(j,c.idx(tt)-darkidxs(j)+preidx) = c.rmsCh1(tt);
        
end

meany.fft(meany.fft == 0) = NaN;
meany.rms(meany.rms == 0) = NaN;

subplot(211); plot(1-preidx:postidx, nanmean(meany.fft), 'r-', 'LineWidth', 2);
subplot(212); plot(1-preidx:postidx, nanmean(meany.rms), 'b-', 'LineWidth', 2);

% save 6789May2020.mat c meany
