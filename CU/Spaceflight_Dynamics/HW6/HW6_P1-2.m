clear
clc
close all
addpath('../../bin')
% ------------------------------------------------------------------------
%%% Givens
% ------------------------------------------------------------------------
u = 398600; % km^3 / s^2
wE = 7.2921158553e-5; % rad/s

% ------------------------------------------------------------------------
%%% Determining Orbital Elements
% ------------------------------------------------------------------------
[a, e, i, raan, w, M, n, year, day] = tleParser('ISS_Zarya.txt');

% ------------------------------------------------------------------------
%%% #1 - epoch
% ------------------------------------------------------------------------
% 279 = 10-06-2016
% .9209 = 22:06:5.76
%
% 10-06-2016 22:06:5.76 UTC
% 10-06-2016 16:06:5.76 MT
%
% Must propogate for 4 days and 02:53:54.24
% = 10434.24 sec
%
% ------------------------------------------------------------------------
%%% Preparing to propogate
% ------------------------------------------------------------------------
%%% Uniquely storing the initial true anomaly (from 10-10-2016 19:00:00 MT)
Mi = M + n*10434.24;

%%% Setting initial gst 
% Initial time in UTC = 10-06-2016 22:06:5.76
y = 2016;
m = 10;
d = 6.9209;
% Calculating Julian epoch
J0 = 367*y - round(7*(y+round((m+9)/12))/4) + round(275*m/9) + d +...
    1721013.5;
T_ut1 = (J0-2451545.0)/36525;
% Converting to GST (deg)
gst0_deg = 100.4606184 + 36000.77005361*T_ut1 + 0.00038793*T_ut1^2 -...
    (2.6e-8)*T_ut1^3;
% Converting to GST (rad)
gst0 = gst0_deg * pi / 180;
%%% Set a 3 hour time
ti = 0;
tf = 3*60*60; % sec

%%% Create structures to hold lat/lon data
lat_prop = zeros(tf/60,1);
lon_prop = zeros(tf/60,1);

% ------------------------------------------------------------------------
%%% Propogate state over time and determine lat/lon along the way
% ------------------------------------------------------------------------
% Iteration counter
k = 1;
for t = ti:60:tf
    % Calculate new Mean, Eccentric, and True anomalies, and gst
    M_t = Mi + n*t; % rads
    E_t = M2E(M_t,e); % rads
    ta_t = E2T(E_t,e); % rads
    gst_t = gst0 + wE * t; % rads
    
    % Calculate ECI state from orbital elements
    [rECI_t, v_t] = OE2ECI(a, e, i, raan, w, ta_t, u); % km, km/s
    
    % Convert ECI position to ECEF
    [r_t] = ECI2ECEF(rECI_t, gst_t); % km
    
    % Calculate latitude and longitude from state
    [lat_t, lon_t] = r2latlon(r_t); % rads
    
    % Convert geodetic to geocentric latitude
    lat_t = atan(tan(lat_t)*(1-0.081819221456^2)); % rads
    
    % Store propogated latitude and longitude
    lat_prop(k,1) = lat_t * 180 / pi; % deg
    lon_prop(k,1) = lon_t * 180 / pi; % deg
    k = k + 1;
end

%%% Plot longitude vs latitude
hold on
load worldmap2384.dat; % read file
x = worldmap2384(:,1);
y = worldmap2384(:,2);
plot(x,y)
plot(lon_prop,lat_prop,'x')



% ========================================================================
%%% #2
% ------------------------------------------------------------------------

% ------------------------------------------------------------------------
%%% Givens
% ------------------------------------------------------------------------
B_lat = 40.01; % deg
B_lon = 254.83-360; % deg
B_alt = 1.615; % km


% ------------------------------------------------------------------------
%%% Propogate orbit and determine azimuth and elevation
% ------------------------------------------------------------------------
az_prop = zeros(tf/60,1);
el_prop = zeros(tf/60,1);
% Iteration counter
k = 1;
for t = ti:60:tf
    % Calculate new Mean, Eccentric, and True anomalies, and gst
    M_t = Mi + n*t; % rads
    E_t = M2E(M_t,e); % rads
    ta_t = E2T(E_t,e); % rads
    gst_t = gst0 + wE * t; % rads
    
    % Calculate ECI state from orbital elements
    [rECI_t, v_t] = OE2ECI(a, e, i, raan, w, ta_t, u); % km, km/s
    
    % Convert ECI position to ECEF
    [r_t] = ECI2ECEF(rECI_t, gst_t); % km
    
    % Find elevation of ISS relative to Boulder
    [az_el_p] = ECEF2az_el_p(r_t, B_lat, B_lon, B_alt)
    
    % Store elevation and azimuth
    az_prop(k,1) = az_el_p(1);
    el_prop(k,1) = az_el_p(2);
end

figure
polarplot(az_prop,el_prop)













