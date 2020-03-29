
%% Set up the DAQ
s = daq.createSession('ni');
s.addAnalogInputChannel('Dev1', 0, 'voltage');
s.addAnalogInputChannel('Dev1', 1, 'voltage');
    s.Rate = 20000;
    s.DurationInSeconds = 60;
s.NotifyWhenDataAvailableExceeds = s.Rate * s.DurationInSeconds;

lh = s.addlistener('DataAvailable', @listentothis);

%% Set up the camera
    NET.addAssembly('C:\Program Files\Thorlabs\Scientific Imaging\DCx Camera Support\Develop\DotNet\uc480DotNet.dll');
    cam = uc480.Camera;
    cam.Init(0);
    cam.Display.Mode.Set(uc480.Defines.DisplayMode.DiB);
    cam.PixelFormat.Set(uc480.Defines.ColorMode.RGBA8Packed);
    cam.Trigger.Set(uc480.Defines.TriggerMode.Software);
    [~, MemId] = cam.Memory.Allocate(true); 
    [~, Width, Height, Bits, ~] = cam.Memory.Inquire(MemId); 
    cam.Timing.Exposure.Set(8);

%% Loop to collect data

for j = 1:144
    
    % Start the DAQ in the background
    s.startBackground();

    % 60 seconds of collection
    for k = 1:6
        %pause(5)
        % Acquire and process the photo data
        cam.Acquisition.Freeze(uc480.Defines.DeviceParameter.Wait);
        [~, tmp] = cam.Memory.CopyToArray(MemId);
        vData = reshape(uint8(tmp), [Bits/8, Width, Height]);
        vData = vData(1:3, 1:Width, 1:Height);
        vData = permute(vData, [3,2,1]);
        % himg = imshow(vData); 
        imageFileName = sprintf('GallmanImage_%s.mat', datestr(now, 'mm-dd-yyyy_HH-MM-SS'));
        save(imageFileName, 'vData');        
        pause(10);
    end

    % Here is the pause between samples
    pause(540);
    
end
