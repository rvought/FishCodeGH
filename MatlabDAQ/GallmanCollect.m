% open the session with the device NIDAQ USB-6363 works
s = daq.createSession('ni');
% Add two channels
s.addAnalogInputChannel('Dev1', 0, 'voltage');
s.addAnalogInputChannel('Dev1', 1, 'voltage');
    s.Rate = 20000;
    s.DurationInSeconds = 30;

% This is when Matlab activates the listener
s.NotifyWhenDataAvailableExceeds = s.Rate * s.DurationInSeconds;

% The listener does stuff with the data
lh = s.addlistener('DataAvailable', @listentothis);

% Schedule the events.
for j = 1:2
    s.startBackground(); % Starts data collection
    pause(60); % RACE CONDITION!
end
