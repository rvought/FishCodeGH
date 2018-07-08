function [out, dat] = sta(spikes, stim, Fs, window)
% out = sta(spikes, stim, Fs)
% spikes is a list of spike times (in seconds)
% stim is the stimulus samples
% Fs is the sample rate for the stimulus

%% Setup

if nargin == 4; pre = window(1); post = window(2); end;
if nargin == 3; pre = 0.050; post = 0.010; end;

% This is adding a pad for the end...
tailor = zeros(1,Fs);
stim = [stim tailor];


 %% Do a random permutation of the spikes - same number with the same distribution of ISIs
  spikeintervals = diff(spikes);
  randspikeintervals = spikeintervals(randperm(length(spikeintervals)));
  randspiketimes(1) = spikes(1);
  for i=1:length(randspikeintervals); randspiketimes(end+1) = randspiketimes(end) + randspikeintervals(i); end;

% A Time series
tim = 1/Fs:1/Fs:length(stim)/Fs;

vel = diff(stim(1:10:end)); vel(end+1) = vel(end);
vtim = tim(1:10:end);
acc = diff(vel(1:2:end)); acc(end+1) = acc(end);
atim = vtim(1:2:end);



% Time series for the output
out.tim = -pre:1/Fs:post;
 out.tim = out.tim(1:end-1);
     len = length(out.tim);
%     out.avg = zeros(1,len);
%     dat(:,1) = zeros(1,len);

%% Extract chunks of stimuli around each spike
% Get the first datum

% What is the first useful spike?

firstspikeidx = find(spikes > pre, 1);

    out.avg = stim(tim >= spikes(firstspikeidx) - pre & tim <= spikes(firstspikeidx) + post);
    out.rand = out.avg;
    out.vel = vel(vtim >= spikes(firstspikeidx) - pre & vtim <= spikes(firstspikeidx) + post);
    out.rvel = out.vel;
    out.acc = acc(atim >= spikes(firstspikeidx) - pre & atim <= spikes(firstspikeidx) + post);
    out.racc = out.acc;
    
% Cycle through the rest of the data

for j = length(spikes):-1:firstspikeidx+1;
    newone = stim(tim >= spikes(j) -pre & tim <= spikes(j) +post);
        newrandone = stim(tim >= randspiketimes(j) -pre & tim <= randspiketimes(j) +post);
    newone = newone(1:len);
        newrandone = newrandone(1:len);
    out.avg = out.avg + newone;
        out.rand = out.rand + newrandone;
    dat(:,j) = newone;
    rdat(:,j) = newrandone;
    
    newvel = vel(vtim >= spikes(j) -pre & vtim <= spikes(j) +post);
        newrvel = vel(vtim >= randspiketimes(j) -pre & vtim <= randspiketimes(j) +post);
    out.vel = out.vel + newvel;
        out.rvel = out.rvel + newrvel;
        
    newacc = acc(atim >= spikes(j) -pre & atim <= spikes(j) +post);
        newracc = acc(atim >= randspiketimes(j) -pre & atim <= randspiketimes(j) +post);
    out.acc = out.acc + newacc;
        out.racc = out.racc + newracc;
           
end;

%% Calculate

    % Divide by the number of spikes that we have
    out.avg = out.avg / length(spikes);
    out.rand = out.rand / length(spikes);
    out.vel = out.vel / length(spikes);
    out.rvel = out.rvel / length(spikes);
    out.acc = out.acc / length(spikes);
    out.racc = out.racc / length(spikes);
    out.vtim = out.tim(1:10:end);
    out.atim = out.vtim(1:2:end);
    
    
figure(1); clf; 
ax(1) = subplot(311); plot(out.tim, out.avg, out.tim, out.rand);
ax(2) = subplot(312); plot(out.vtim, out.vel, out.vtim, out.rvel);
ax(3) = subplot(313); plot(out.atim, out.acc, out.atim, out.racc);
linkaxes(ax, 'x');
    
    
    
    
    
    % Get the STD
    for k = length(out.avg):-1:1;
        out.std(k) = std(dat(k,:));
    end;

    % Get the peak index in the average STA
    [~, idx] = max(abs(out.avg));

    pos(1) = -1; neg(1) = -1;

for p = 1:length(spikes);
    % figure(4); hold on; plot(dat(:,p));
    if dat(idx,p) > 0; pos = [pos p]; end;
    if dat(idx,p) < 0; neg = [neg p]; end;
end;

 pos = pos(2:end);
 %neg = neg(2:end);
 
 out.pos = dat(:,pos(1)); 
 %out.neg = dat(:,neg(1));
 
 for l = 2:length(pos);
     out.pos = out.pos + dat(:,pos(l));
 end;
 %for l = 2:length(neg);
 %    out.neg = out.neg + dat(:,neg(l));
 %end;
 
 out.pos = out.pos / length(pos);
 
 
 
 
 
 
 
%  %% Do a random permutation %%% CODE SHOULD GO SOMEWHERE ELSE
%   spikeintervals = diff(spikes);
%  randspikeintervals = spikeintervals(randperm(length(spikeintervals)));
%  
%  figure(3);
%  subplot(311); hold on; % Plot original spikes
%  for i=1:length(spikes); plot([spikes(i) spikes(i)], [0 1]); end;
%  
%  % Reconstruct spikes from spikeintervals
%  
%  newspiketimes(1) = spikes(1);
%  for i=1:length(spikeintervals); newspiketimes(end+1) = newspiketimes(end) + spikeintervals(i); end;
% subplot(312); hold on; % Plot reconstructed original spikes
% for i=1:length(newspiketimes); plot([newspiketimes(i) newspiketimes(i)], [0 1]); end;
%  
%  randspiketimes(1) = spikes(1);
%  for i=1:length(randspikeintervals); randspiketimes(end+1) = randspiketimes(end) + randspikeintervals(i); end;
% subplot(313); hold on; % Plot randomized spikes
%  for i=1:length(randspiketimes); plot([randspiketimes(i) randspiketimes(i)], [0 1]); end;
%   
 
% out.neg = out.neg / length(neg);
% out.pos = mean(foo(pos).raw);
% out.neg = mean(foo(neg).raw);
% 
% 
% 
% figure(1);
% subplot(211); plot(out.tim, out.avg);
% hold on; 
% plot(out.tim, out.pos, 'k');
% plot(out.tim, out.neg, 'r');
% subplot(212); plot(out.tim, out.std);
% figure(2);
% hold on;
% for k=1:length(foo);
%     plot(out.tim, foo(k).raw);
% end;

%average
