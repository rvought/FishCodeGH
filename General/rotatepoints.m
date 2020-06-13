function out = rotatepoints(dat, deg)
% dat are x,y values
% deg is rotation in degrees

%% Prepare data
% STRATEGY 1: Calculate centroid for center of rotation

    poly = polyshape(dat);
    [xC,yC] = centroid(poly);
    
% STRATEGY 2: Take the third point as the center of rotation

    tt = dat(3,:);    
    
% Convert degrees to radians 

    rad = deg2rad(deg);

%% Rotate!

    out.deg = deg; % Copy of rotation for your records
    out.centroidrotate = rotatorcuff(dat, [xC,yC], rad);
    out.thirdrotate = rotatorcuff(dat, tt, rad);
    
%% Plot
figure(1); clf;

ax(1) = subplot(121); title('Centroid'); hold on; 
    plot(dat(:,1), dat(:,2), 'b-*'); % Original data
    plot(out.centroidrotate(:,1), out.centroidrotate(:,2), 'g-*'); 
    plot(xC, yC, 'r*');

ax(2) = subplot(122); title('Third'); hold on; 
    plot(dat(:,1), dat(:,2), 'b-*'); % Original data
    plot(out.thirdrotate(:,1), out.thirdrotate(:,2), 'g-*'); 
    plot(tt(1),tt(2), 'r*');

linkaxes(ax, 'xy'); 


%% Embedded rotation function
function rot = rotatorcuff(data, cent, degR)
% data is the x,y values
% cent is the rotation center [x,y]
% degR is the rotation angle in radians

    center = repmat([cent(1); cent(2)], 1, length(data(:,1))); % Make a matrix for the center of rotation

    R = [cos(degR) -sin(degR); sin(degR) cos(degR)]; % Create rotation matrix
    
    rot = (R * (data' - center)) + center; % Rotate points
    rot = rot'; % Because I am lousy at coding

end

end