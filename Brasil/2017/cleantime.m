function out = cleantime(tracks, particle, plotornot)
% Usage out = cleantime(tracks, particle)


% Colors for plotting
    clrs(1,:)='k*-'; clrs(2,:)='r-*'; clrs(3,:)='c-*'; 
    clrs(4,:)='g-*'; clrs(5,:)='b-*'; clrs(6,:)='m-*';
    clrs(7,:)='y-*';
% Need to find a bigger color palate to use - this sucks. **************

%% How many fish do we have?
    uniqfish = unique([tracks.id]);
    fprintf('Manu found %i fish. \n', length(uniqfish));

% Extract the indices for each fish from Manu's painful structure
    foo = [tracks.t]; goo = [tracks.f1]; % Just to extract times and freqs

    for i= 1:length(uniqfish)
        fishidxs{i} = find([tracks.id] == uniqfish(i));
            tim{i} = foo(fishidxs{i}); % Get the relevant times
            freq{i} = goo(fishidxs{i}); % Get the relevant frequencies
    end;
    clear foo goo;
    
%% Separate the frequency data for each fish
for k = 1:length(uniqfish) % For each fish

    % Reset the temporary variables.
    maxtim = 0; % Keep track of the lastest time for the current fish
    keeptims = []; % We are only keeping time points for which a frequency was registered

    for i = 1:length(fishidxs{k}) % For each data point for the current fish

        if tim{k}(i) > maxtim % Make sure we are advancing in time
            maxtim = tim{k}(i);
                if isnan(freq{k}(i)) == 0 % Only take data where there is a real frequency measurement
                    keeptims = [keeptims i];
                end
        end
    end

    out(k).tim = tim{k}(keeptims);
    out(k).freq = freq{k}(keeptims);

end

% Plot all of the frequency data.
    if nargin == 3
        figure(1); clf; hold on; for i = 1:length(out); plot(out(i).tim, out(i).freq, '*'); end;
    end
%% Align the frequency and position data - Extract position data

for ii = 1:length(uniqfish) % For each fish...
    
   for jj = 1:length(out(ii).tim)
       out(ii).posidx(jj) = find(particle.t == out(ii).tim(jj));
   end;
   
        % These are the Xs and Ys for all times - sample and hold
       out(ii).x = particle.fish(ii).x; 
       out(ii).y = particle.fish(ii).y;

        % These are the positions ONLY when there was a detected EOD
        % frequency.  These are the data that should be used.
       out(ii).tx = particle.fish(ii).x(out(ii).posidx)';
       out(ii).ty = particle.fish(ii).y(out(ii).posidx)';
       out(ii).tt = particle.t(out(ii).posidx);

       if nargin == 3
        figure(2); clf; hold on; for i=1:length(out); plot(out(i).tx, out(i).ty, '*'); end;
       end
       
       
end



%% Divide each fish into contingious temporal epochs.  This is important for finding errors and calculating speeds/distances 
       % THIS SHOULD BE DONE AFTER CONCATENATION OF FILES
% % SET THE MINIMUM DURATION FOR AN EPOCH
% 
% mindur = 5; % Number of contiguous samples to be considered an epoch
% 
% for kk = 1:length(uniqfish) % For each fish...
%     
%     samplen = mode(diff(out(kk).tt)); % Get the mode intersample interval - presumably at the sample rate
%   
%     % Get the starts and ends of the epochs
%     z = zeros(1,length(out(kk).tt)); 
%     z(out(kk).tt == samplen) = 1;
%     z = diff(z);
%         zups = find(z == 1); % Start indices
%         zdns = find(z == -1); % End indicess
%         
%         if zdns(1) < zups(1) % We need to figure out how to handle an epoch that started before the start of the sample.
%             zdns(1:end-1) = zdns(2:end); zdns = zdns(1:end-1);
%         end
%         
%         for pp = 1:length(zups);
%         end
%         
% end

       
       
%% Get the instantaneous velocity 
    % THIS TOO SHOULD HAPPEN AFTER THE CONCATENATION
%        timstep = out(1).tt(5) - out(1).tt(4);
% %       out(ii).idist = -0.1 * ones(1, length(out(ii).tx));
%        out(ii).idist(1) = 0;
%        for jj = 2:length(out(ii).tt) 
% %           if (out(ii).tt(jj) - out(ii).tt(jj-1)) == timstep;
%                out(ii).idist(jj) = sqrt((out(ii).tx(jj-1)-out(ii).tx(jj))^2 + (out(ii).ty(jj-1)-out(ii).ty(jj))^2);
% %           end
%        end
%        
% end


% figure(3); clf; hold on; for i=1:length(out); plot(out(i).tt, out(i).idist, '*'); end;


% out.tim = tim(keeptims);
% out.freq = freq(keeptims);

%>> %for j = 1:7; out(j) = cleantime(tim(fish{j}), freq(fish{j})); end;
%>> for j = 1:7; out(j) = cleantime(tim(fish{j}), freq(fish{j})); end; 
