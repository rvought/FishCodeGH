function out = makestim(in, Fs, offset, tailval, filename)
% makestim(in, Fs, offset, filename)
% "in" is position signal in cm

z = zeros(1,ceil(Fs*0.1));

sclr = 0.00628; % Measured at 1Hz, 2V pk-pk (velocity out)

% First derivative

out = diff(in)*0.8 / sclr;

out = out + (offset / 10000);

out(end+1) = offset/10000;

z = z + tailval/10000;

wavwrite([out z], Fs, 16, filename);



