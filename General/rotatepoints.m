function out = rotatepoints(dat, deg)
% Usage: out = rotatepoints(dat, deg)
% where dat is x,y pairs (Column 1 is X, Column 2 is Y)
% deg is desired rotation in degrees
% out is a structure with the rotated XY coordinates using two center points
% out.centroidrotate is the data rotated around the centroid of the original data
% out.thirdrotate is the rotated data using the third point as the center
% There is no check to make sure that you have more than three points. So don't f*ck that up.


%% Prepare data

% STRATEGY 1: Calculate centroid for center of rotation

    poly = polyshape(dat); % Change the data into a Matlab object known as a polyshape
    [xC,yC] = centroid(poly); % centroid calculates centroids.
    
% STRATEGY 2: Take the third point as the center of rotation

    tt = dat(3,:);  % Just take the third XY value in the array. Why not(?), unless you have fewer than three points.
    
% Convert degrees to radians 

    rad = deg2rad(deg);

%% Rotate! Using the embedded function below, rotatorcuff

    out.deg = deg; % Copy of rotation for your records
    out.centroidrotate = rotatorcuff(dat, [xC,yC], rad); % Rotation around centroid 
    out.thirdrotate = rotatorcuff(dat, tt, rad); % Rotation around 3rd XY value
    
%% Plot the results

figure(1); clf;

ax(1) = subplot(121); title('Centroid'); hold on; 
    plot(dat(:,1), dat(:,2), 'b-*'); % Original data in blue
    plot(out.centroidrotate(:,1), out.centroidrotate(:,2), 'g-*'); % Rotated data in green
    % out.centroidrotate(:,1) are the X values, out.centroidrotate(:,2) are Ys
    plot(xC, yC, 'r*'); % Red dot for the center

ax(2) = subplot(122); title('Third'); hold on; 
    plot(dat(:,1), dat(:,2), 'b-*'); % Original data in blue
    plot(out.thirdrotate(:,1), out.thirdrotate(:,2), 'g-*'); % Rotated data in green
    plot(tt(1),tt(2), 'r*'); % Red dot for the center

linkaxes(ax, 'xy'); 


%% Embedded rotation function (You could make this a separate function)

function rot = rotatorcuff(data, cent, degR)
% data is the x,y values
% cent is the rotation center [x,y]
% degR is the rotation angle in radians

    centermatrix = repmat([cent(1); cent(2)], 1, length(data(:,1))); % Make a matrix for the center of rotation

    R = [cos(degR) -sin(degR); sin(degR) cos(degR)]; % Create rotation matrix
    
    rot = (R * (data' - centermatrix)) + centermatrix; % Rotate points
    rot = rot'; % Because I am lousy at coding

end

end