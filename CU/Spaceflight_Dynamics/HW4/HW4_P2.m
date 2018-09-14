clear
clc
% Givens
r_ECI = [-5634; -2645; 2834]; % km
gst   = 82.75 * pi / 180; % rads
rE    = 6378.1363; % km

% Rotate to get r_ECEF
r_ECEF = ECI2ECEF(r_ECI, gst) % km

% Compute coordinates and r
x = r_ECEF(1); % km
y = r_ECEF(2); % km
z = r_ECEF(3); % km
r = norm(r_ECEF); % km

% Compute lat, lon, and altitude 
lat = asin(z / r); % rads
lon = pi + atan(y/x); % (added pi to correct for sign) rads
alt = r - rE; % km

% Print results
fprintf('Latitude  = %3.5f deg\n', lat * 180 / pi)
fprintf('Longitude = %3.5f deg\n', lon * 180 / pi)
fprintf('Altitude  = %3.5f km\n', alt)