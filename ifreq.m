function z = ifreq(y, Fs)
% ifreq computes zero crossings accurately by interpolation between
% samples. Also reports instanteous frequency. To keep things simple, the
% threshold is assumed to be zero.
% z = ifreq(signal, samplerate);
% z.tim are zero crossing times
% z.freq are the instanteous freqs
% z.diff are the intervals

%% We could make this user defined, but it is easy enough to work with
% this set to zero by default.

    threshold = 0;

%% This gives us postive going zerocrossings
% We take the samples 1 to end-1 and compare that to 2 to end - offset by
% one sample and find cases where value is below threshold and next above.

    pres = y(1:end-1); posts = y(2:end);
    zcs_idx = find((pres <= threshold) & (posts > threshold));

%% Calculate the actual zero crossing between the samples

% How large is the step from below to above the threshold
    amplitude = y(zcs_idx+1)-y(zcs_idx);   
% Assume a line between these points and find what percentage 
    pct = (threshold - y(zcs_idx)) ./ amplitude;    
% So now we add the fraction to the index
    fract = zcs_idx + pct;                

%% Output the data into our structure z   
    
% Inter-Zercrossing-Interval    
    z.diff = diff(fract)/Fs;
% The inverse, which is the instanteous frequency
    z.freq = 1./(diff(fract)/Fs);     
% And the times at which these events occur
    z.tims = cumsum(diff(fract)/Fs)+fract(1)/Fs; 
