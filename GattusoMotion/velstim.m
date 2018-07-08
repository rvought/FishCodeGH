clear all;

Fs = 1000;


tim001 = 1/Fs:1/Fs:0.1;
tim002 = 1/Fs:1/Fs:0.2;
tim005 = 1/Fs:1/Fs:0.5;
tim010 = 1/Fs:1/Fs:1.0;
tim020 = 1/Fs:1/Fs:2.0;
tim040 = 1/Fs:1/Fs:4.0;
tim050 = 1/Fs:1/Fs:5.0;
tim100 = 1/Fs:1/Fs:10.0;

% 500 msec of zeros (for beginning and end of stimuli)
zipper = zeros(1, length(tim005));

% Initial offset and reset from/to zero
        cvinddn = 0:-1/length(tim020):-1;
        cvinddn = cvinddn(1:length(tim020));
        cvendup = cvinddn(end:-1:1);
        
% Top and bottom for CV stimuli        

    top = ones(1, length(tim002));
    bottom = -top;

% Initial stimuli - 10cm (+/- 5cm) constant velocity
% 1 V per cm.  1 is 5V, -1 is -5V.

%% 10cm / s  

    % Our desired constant velocity
    cv10up = -1:2/length(tim010):1;
        cv10up = cv10up(1:length(tim010));
    cv10dn = 1:2/-length(tim010):-1;
        cv10dn = cv10dn(1:length(tim010));

        
    ourcycle = [bottom cv10up top cv10dn];
    cv10stim(:,1) = [zipper cvinddn ourcycle ourcycle ourcycle ourcycle ourcycle bottom cvendup zipper];
    
    stimtim = 1/Fs:1/Fs:length(cv10stim)/Fs;
    
    subplot(411); plot(stimtim, cv10stim);
    
    
%% 5cm / s

    % Our desired constant velocity
    cv05up = -1:2/length(tim020):1;
        cv05up = cv05up(1:length(tim020));
    cv05dn = 1:2/-length(tim020):-1;
        cv05dn = cv05dn(1:length(tim020));

        
    ourcycle = [bottom cv05up top cv05dn];
    cv05stim(:,1) = [zipper cvinddn ourcycle ourcycle ourcycle ourcycle ourcycle bottom cvendup zipper];
    
    stimtim = 1/Fs:1/Fs:length(cv05stim)/Fs;
    subplot(412); plot(stimtim, cv05stim);
    

%% 1cm / s

    % Our desired constant velocity
    cv01up = -1:2/length(tim100):1;
        cv01up = cv01up(1:length(tim100));
    cv01dn = 1:2/-length(tim100):-1;
        cv01dn = cv01dn(1:length(tim100));

        
    ourcycle = [bottom cv01up top cv01dn];
    cv01stim(:,1) = [zipper cvinddn ourcycle ourcycle bottom cvendup zipper];
    
    stimtim = 1/Fs:1/Fs:length(cv01stim)/Fs;
    subplot(413); plot(stimtim, cv01stim);


% Tuberous, constant amplitude

tstim = ones(1, length(cv05stim));
cv05stim(:,2) = tstim;
% audiowrite('cv05tc1.wav',cv05stim,Fs);

tstim = ones(1, length(cv10stim));
cv10stim(:,2) = tstim;
% audiowrite('cv10tc1.wav',cv10stim,Fs);

tstim = ones(1, length(cv01stim));
cv01stim(:,2) = tstim;
% audiowrite('cv01tc1.wav',cv01stim,Fs);

% Tuberous, oscillation 20 Hz

AMFreq = 20;

ttim = 1/Fs:1/Fs:length(cv05stim)/Fs;
tstim = 0.175 + (sin(2*pi*ttim*AMFreq)/40); 
cv05stim(:,2) = tstim;
% audiowrite('cv05tf20.wav',cv05stim,Fs);

ttim = 1/Fs:1/Fs:length(cv10stim)/Fs;
tstim = 0.175 + (sin(2*pi*ttim*AMFreq)/40); 
cv10stim(:,2) = tstim;
audiowrite('cv10tf20.wav',cv10stim,Fs);

ttim = 1/Fs:1/Fs:length(cv01stim)/Fs;
tstim = 0.175 + (sin(2*pi*ttim*AMFreq)/40); 
cv01stim(:,2) = tstim;
% audiowrite('cv01tf20.wav',cv01stim,Fs);

