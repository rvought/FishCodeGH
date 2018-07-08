% ESTIMATEAMPLITUDES
%
% Code to estimate amplitudes from spatial tracks.
%
% Manu S. Madhav
% 19-May-2017

[particleFileName, pathName] = uigetfile('Choose particle file');

% particleFileName = 'TerraRonca_Calibration_01_100s_particle.mat';
% pathName = '.';

load(fullfile(pathName,particleFileName));

%% First method

K1 = zeros(particle.nFish,1);
for k = 1:particle.nFish

    A = abs(particle.fish(k).ampAct');

     G = [particle.gridCoord, zeros(size(particle.gridCoord,1),1)];
    %G = [particle.gridCoord];
    P = [particle.fish(k).x, particle.fish(k).y, particle.fish(k).z];

    % Pairwise distance between grid electrodes and fish locations 
    R = pdist2(P,G);

    % Fish angle vectors
    F = exp(1i*particle.fish(k).theta);

    % Position vector of fish location relative to grid in the X-Y plane
    V = zeros(size(P,1),size(G,1));
    for p = 1:size(P,1)
        for g = 1:size(G,1)
            V(p,g) = (P(p,1:2) - G(g,1:2))*[1;1i];
        end
    end

    % Relative angle between grid electrodes and fish angles
    T = angle(V./repmat(F,1,size(V,2)));

    % Theoretical amplitude
    B = abs(cos(T)./(R.^2));

    % Linear regression problem
    X = B(:);
    Y = A(:);

    % Remove all the Nans
    idx = ~isnan(X) & ~isnan(Y);
    X = X(idx);
    Y = Y(idx);

    K1(k) = X\Y;
end

%% Second method

K2 = zeros(particle.nFish,1);
for k = 1:particle.nFish
    % Second method
    Y = abs(particle.fish(k).ampAct(:));
    X = abs(particle.fish(k).ampTheor(:));

    % Remove all the Nans
    idx = ~isnan(X) & ~isnan(Y);
    X = X(idx);
    Y = Y(idx);

    K2(k) = X\Y;
end