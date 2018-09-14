clear
clc
% Givens
lat       = 40.01 * pi / 180; % rads
lon       = 254.83 * pi / 180; % rads
alt       = 1.615; % km
gst       = 103 * pi / 180; % rads
rE        = 6378.1363; % km

% Calculating Coordinates in ECEF
r = alt + rE; % km
x = r * cos(lat) * cos(lon); % km
y = r * cos(lat) * sin(lon); % km
z = r * sin(lat); % km

r_ECEF = [x; y; z]; % km

% Rotate to get r_ECI
r_ECI = ECEF2ECI(r_ECEF, gst); %km

% Print results
fprintf('ECI_x = %3.5f km\n', r_ECI(1))
fprintf('ECI_y = %3.5f km\n', r_ECI(2))
fprintf('ECI_z = %3.5f km\n', r_ECI(3))