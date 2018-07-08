function signal = physfilt( ndata, Fs, hilord );
% Signal = physfilt( signal, Fs, [high low order] );
% signal is the data
% Fs is the samplerate in Hz
% Optional... 
%     "low" freq cutoff for lowpass
%     "high" freq cutoff of highpass
%     order of filter (odd number, usual 3-7)
% Default is  low=4500, high=750, order=7

signal = ndata - mean(ndata);

if nargin == 2 
	order = 7; low = 1500; high = 500;
%	order = 7; low = 3600; high = 1000;
%	order = 7; low = 1400; high = 400;
end

if nargin == 3
	high = hilord(1); low = hilord(2); order = hilord(3);
end

        Wnlow = low*2/Fs;
        Wnhigh = high*2/Fs;
        [bb,aa] = butter(order,Wnlow,'low');
        [dd,cc] = butter(order,Wnhigh,'high');
        signal = filtfilt(bb,aa,signal);
        signal =filtfilt(dd,cc,signal);
