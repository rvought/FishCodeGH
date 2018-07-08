function out = cereal(param)

param = param+2;

s = serial('COM5');
%set(s,'BaudRate',115200);
set(s,'BaudRate',19200);

fopen(s); 
fprintf(s,'t 0');
fprintf(s,'s r0xc8 256');
fprintf(s,'s r0xca 100000');
fprintf(s,'s r0xcb 200000');
fprintf(s,'s r0xcc 200000');
fprintf(s,'s r0xcd 200000');
fprintf(s,'s r0x24 21');
fprintf(s,'t 1');
fprintf(s,'r0x24 0');

fprintf(s,'s r0x19 20000'); % ok Set analog scaling to 4000 counts per 10V.
fprintf(s,'s r0xcb 700000'); % ok Set velocity to 7000 counts/second
fprintf(s,'s r0xcc 200000'); % ok Set acceleration to 200000 counts/second2
fprintf(s,'s r0xcd 200000'); % ok Set deceleration to 200000 counts/second2
fprintf(s,'s r0x24 22'); % ok Amplifier set in Analog Position Mode
fprintf(s,'t 1'); % ok This command will guarantee all new move parameters are in effect.





fprintf(s, 'g r0x32');
out = fscanf(s);
fclose(s);
delete(s);
clear s
