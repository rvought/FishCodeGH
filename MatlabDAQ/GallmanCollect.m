s = daq.createSession('ni');
s.addAnalogInputChannel('Dev1', 0, 'voltage');
s.addAnalogInputChannel('Dev1', 1, 'voltage');
    s.Rate = 20000;
    s.DurationInSeconds = 30;
s.NotifyWhenDataAvailableExceeds = s.Rate * s.DurationInSeconds;

lh = s.addlistener('DataAvailable', @listentothis);

for j = 1:2
    s.startBackground();
    pause(60); % RACE CONDITION!
end
