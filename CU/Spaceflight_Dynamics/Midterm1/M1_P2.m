clear
clc
addpath('../../bin')
rad2deg = 180/pi;
deg2rad = pi/180;
% ------------------------------------------------------------------------
%%% Givens
% ------------------------------------------------------------------------
r_ECI=[8800;-1100;5500]; % km
dt = 6*3600; % sec
% Earth angular velocity
wE = 15*deg2rad/3600; % rad/s
% Earth radius
rE = 6378.1363; % km

%%% Boulder
B_lat = 40.015 * deg2rad; % rad
B_lon = -105.270 * deg2rad; % rad
B_r = 6379.77; % km
B_alt = B_r - rE; % km

% ------------------------------------------------------------------------
%%% Calculating lat/lon at time t
% ------------------------------------------------------------------------
% Rotation angle (GMT at t = 18:00:00)
dGST = wE * dt;%(dt + 12*3600); % rads
r_ECEF = ECI2ECEF(r_ECI, dGST);

[lat, lon] = ECEF2latlon(r_ECEF); % rad
latd = lat*rad2deg; % deg
lond = lon*rad2deg; % deg

% ------------------------------------------------------------------------
%%% Calculating azimuth and range at time t
% ------------------------------------------------------------------------
az_el_p = ECEF2az_el_p(r_ECEF, B_lat, B_lon, B_alt);

%%% Printing results
fprintf('At time t (18:00:00 GMT), the satellite has the following properties\n')
table(latd, lond)
fprintf('Azimuth = %f deg\n', az_el_p(1)*rad2deg)
fprintf('Elevation = %f deg\n', az_el_p(2)*rad2deg)
