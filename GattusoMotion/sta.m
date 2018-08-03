function out = sta(spikes, stim, Fs, window)
% out = sta(spikes, stim, Fs)
% spikes is a list of spike times (in seconds)
% stim is the stimulus samples
% Fs is the sample rate for the stimulus

%% Setup

if nargin == 3; pre = 1; post = 0.2; end % Default pre is 200 msec pre spike and 20 msec post
if nargin == 4; pre = window(1); post = window(2); end % Usually don't use this - usually all of the data

% This is adding a pad for the end so that the STA go to the last spike...
    tailor = zeros(1,Fs);
    stim = [stim tailor];


%% Do a random permutation of the spikes - same number with the same distribution of ISIs
  spikeintervals = diff(spikes);
  randspikeintervals = spikeintervals(randperm(length(spikeintervals)));
  randspiketimes(1) = spikes(1);
  for i=1:length(randspikeintervals) 
      randspiketimes(end+1) = randspiketimes(end) + randspikeintervals(i); 
  end

%% Handle the stimulus

    tim = 1/Fs:1/Fs:length(stim)/Fs; % A time series for the original stimulus

% Time series for the output

    out.tim = -pre:1/Fs:post; % Time for the STA
    out.tim = out.tim(1:end-1); % For some reason shorten it by one sample... 
    
    len = length(out.tim);
    out.avg = zeros(1,len); % Initialize the STA
    out.ravg = out.avg;
    out.std = out.avg;
    out.rstd = out.avg;

%% Extract chunks of stimuli around each spike

    spikes = spikes(spikes > pre); % Make sure the first spike happens with complete pre.
    randspiketimes = randspiketimes(randspiketimes > pre); % Make sure the first spike happens with complete pre.
    

for j = length(spikes):-1:1

% Position    
    % Spikes
    newone = stim(tim >= spikes(j) -pre & tim <= spikes(j) +post); % Get stimulus chunk
    newone = newone(1:len); % Trim to correct numer of samples (our 'find' command above will have some jitter)
    out.avg = out.avg + newone;

%     % Shuffled spikes
%     newrandone = stim(tim >= randspiketimes(j) -pre & tim <= randspiketimes(j) +post);
%     newrandone = newrandone(1:len);
%     out.ravg = out.ravg + newrandone;
%     
    dat(:,j) = newone;
%     rdat(:,j) = newrandone;
%            
end

%% Calculate

    % Divide by the number of spikes that we have
    out.avg = out.avg / length(spikes);
%     out.ravg = out.ravg / length(spikes);
        
figure; clf; 
    ax(1) = subplot(211); plot(out.tim, out.avg, 'LineWidth', 2); 
%     hold on; plot(out.tim, out.ravg);
    
        
    % Get the STD
    for k = length(out.avg):-1:1
        out.std(k) = std(dat(k,:));
    end
%     % Get the randSTD
%     for k = length(out.ravg):-1:1
%         out.rstd(k) = std(rdat(k,:));
%     end

    ax(2) = subplot(212); plot(out.tim, out.std, 'LineWidth', 2); 
    linkaxes(ax, 'x');
    xlim([-pre, post]);
    
    
%     % Get the peak index in the average STA
%     [~, idx] = max(abs(out.avg));
% 
%     pos(1) = -1; neg(1) = -1;
% 
% for p = 1:length(spikes)
%     % figure(4); hold on; plot(dat(:,p));
%     if dat(idx,p) > 0; pos = [pos p]; end
%     if dat(idx,p) < 0; neg = [neg p]; end
% end
% 
%  pos = pos(2:end);
%  %neg = neg(2:end);
%  
%  out.pos = dat(:,pos(1)); 
%  %out.neg = dat(:,neg(1));
%  
%  for l = 2:length(pos)
%      out.pos = out.pos + dat(:,pos(l));
%  end
%  %for l = 2:length(neg);
%  %    out.neg = out.neg + dat(:,neg(l));
%  %end;
%  
%  out.pos = out.pos / length(pos);
% end
 
 
 
 
 
 
 
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
