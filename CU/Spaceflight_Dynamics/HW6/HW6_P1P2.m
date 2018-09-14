clear
clc
close all
addpath('../../bin')
% ------------------------------------------------------------------------
%%% #1 Givens
% ------------------------------------------------------------------------
u = 398600; % km^3 / s^2
wE = 7.2921158553e-5; % rad/s

% ------------------------------------------------------------------------
%%% #2 Givens
% ------------------------------------------------------------------------
B_lat = 40.01 * pi / 180; % rad
B_lon = (254.83) * pi / 180; % rad
B_alt = 1.615; % km

% ------------------------------------------------------------------------
%%% Determining Orbital Elements
% ------------------------------------------------------------------------
[a, e, i, raan, w, M, n, year, day] = tleParser('ISS_Zarya.txt');

% ------------------------------------------------------------------------
%%% #1
% ------------------------------------------------------------------------
% epoch1 = 279.9209336 GMT
% 279 = 10-05-2016
% .9209336 = 22:06.1144384   
%
% epoch2 = 10/11/2016 01:00:00 GMT

%%% Start Propogation
JD1 = JD(2016,10,5,22,6.1144384);

%%% End Propogation
JD2 = JD(2016,10,11,1,0); 
dt = (JD2 - JD1) * 60*60*24; % sec
% ------------------------------------------------------------------------
%%% Preparing to propagate
% ------------------------------------------------------------------------
%%% Uniquely storing the initial true anomaly (from 10-11-2016 01:00:00GMT)
Mi = M + n*dt;

%%% Setting initial gst 
T_ut1 = (JD2-2451545.0)/36525;

%%% Converting to GST (deg)
gst0_deg = (67310.54841 + (876600*60*60 + 8640184.812866)*T_ut1 +...
    .093104*T_ut1^2 - 6.2*10^(-6)*T_ut1^3)/240;

%%% Converting to GST (rad)
gst0 = gst0_deg * pi / 180;



% ------------------------------------------------------------------------
%%% #1 Propagate state over time and determine lat/lon along the way
% ------------------------------------------------------------------------
%%% Set a 3 hour time
ti = 0;
tf = 3*60*60; % sec

%%% Create structures to hold lat/lon data
lat_prop = zeros(tf/60,1);
lon_prop = zeros(tf/60,1);

%%% Iteration counter
k = 1;

for t = ti:60:tf
    % Calculate new Mean, Eccentric, and True anomalies
    M_t = Mi + n*t; % rads
    E_t = M2E(M_t,e); % rads
    ta_t = E2T(E_t,e); % rads
    
    % Calculate new gst_t
    gst_t = gst0 + wE * t; % rads
    
    % Calculate ECI state from orbital elements
    [rECI_t, v_t] = OE2ECI(a, e, i, raan, w, ta_t, u); % km, km/s
    
    % Convert ECI position to ECEF
    [r_t] = ECI2ECEF(rECI_t, gst_t); % km
  
    % Calculate latitude and longitude from state
    [lat_t, lon_t] = ECEF2latlon(r_t); % rads
    
    % Convert geodetic to geocentric latitude
    lat_t = atan(tan(lat_t)*(1-0.081819221456^2)); % rads
    
    % Store propagated latitude and longitude
    lat_prop(k,1) = lat_t * 180 / pi; % deg
    lon_prop(k,1) = lon_t * 180 / pi; % deg
    
    % Advance index
    k = k + 1;
end

% ------------------------------------------------------------------------
%%% Plotting Groundtrack
% ------------------------------------------------------------------------
%%% Plot longitude vs latitude
hold on
load worldmap2384.dat; % read file
x = worldmap2384(:,1);
y = worldmap2384(:,2);
plot(x,y)
plot(lon_prop,lat_prop,'rx')
xlim([-180,180])
ylim([-90, 90])
title('ISS Groundtrack')
PlotBoi('Longitude, \circ','Latitude, \circ')

% ------------------------------------------------------------------------
%%% #2 Propogate and analyze azimuths and zeniths of first pass
% ------------------------------------------------------------------------
%%% Iteration counter
k = 1;
%%% Visible moment counter
c = 1;
%%% Logical triggered during first pass
pass_started = false;
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
    [az_el_p] = ECEF2az_el_p(r_t, B_lat, B_lon, B_alt); % km, rad, rad, km
    
    % If the pass has previously started, but now visibility is lost
    if pass_started == true && az_el_p(2) < 0
        break
    end
    % If visible
    if az_el_p(2) > 0
        % Mark that the pass has started
        pass_started = true;
        
        % Keeping azimuth positive
        if az_el_p(1) < 0
            az_el_p(1) = az_el_p(1) + 2*pi;
        end
        
        % Store azimuth, elevation, and zenith
        az_prop(c) = az_el_p(1); % rads
        el_prop(c) = az_el_p(2); % rads
        zenith_prop(c) = pi/2 - el_prop(c); % rads
        visibility(c) = t; % sec
        c = c + 1;
    end
    k = k+1;
end

% ------------------------------------------------------------------------
%%% Plotting polar azimuth/zenith and cartesian visibility (of 1st pass)
% ------------------------------------------------------------------------

figure
polar(az_prop,zenith_prop*(180/pi),'-rx')
set(gcf,'color','white')
title('1st Pass Azimuth and Zenith Angles')


figure
plot([visibility(1):60:visibility(end)]./60,el_prop*180/pi,'-ro',...
    'linewidth',2,'markersize',8)
title('Elevation vs Time of First Pass')
PlotBoi('Time After Epoch, min','Elevation, \circ')
grid on
