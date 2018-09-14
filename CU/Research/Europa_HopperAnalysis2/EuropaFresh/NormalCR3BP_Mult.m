clear
clc
close all
addpath('../ProjectBin')

% ------------------------------------------------------------------------
%%% System Parameters
% ------------------------------------------------------------------------
a = 671100; % semimajor axis (distance-normalizing value)
wz = 2.0472003e-05; % rad/s

%%% Normalizing Factors (real / xN = normalized)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Never felt awesome about these
wNorm = 1/wz; % Angular Velocity (rad/s)
rNorm = a; % Position (km)
vNorm = rNorm/wNorm; % Velocity (km/s)
tNorm = wNorm;

%%% Gravitational Parameters
u1 = 126672520; % km^3 / s^2
u2 = 3203.413216; % km^3 / s^2
u = u2/u1;

%%% Body radii
rad1_n = 69911/rNorm; % radius of primary body (km normalized)
rad2_n = 1560.8/rNorm; % radius of secondary body (km normalized)
% rad1_n = 6371/rNorm; % Earth
% rad2_n = 1737/rNorm; % Moon

%%% Rotating frame cooridinates
rB1_BCR_n = [-u, 0, 0];
rB2_BCR_n = [1-u, 0, 0];

% ------------------------------------------------------------------------
%%% Preparing for Simulation
% ------------------------------------------------------------------------
%%% Setting initial position
lat0 = 60; % deg (-90:90)
lon0 = -45; % deg (-180:180)

%%% Setting initial velocity
v_mag_n = 1/vNorm;

%%% Setting azimuth and elevation ranges
az_range = linspace(0,360,36*2);
el_range = linspace(60,60,1);

%%% Setting simulation preferences
% az = 350; % deg (0:360)
% el = 45; % deg (0:90)

%%% Setting normalized time vector
ti = 0/tNorm;
dt = 1/tNorm;
tf = 50000/tNorm;
time = ti:dt:tf;

%%% Creating storage matrices
States_2_n_MULT = zeros(length(time),6,length(el_range),length(az_range));

% ------------------------------------------------------------------------
%%% Running Simulation
% ------------------------------------------------------------------------
lonMin = 181;
lonMax = -181;
latMin = 91;
latMax = -91;

figure; hold all
for el_i = 1:length(el_range)
    el = el_range(el_i);
    for az_i = 1:length(az_range)
        az = az_range(az_i);
        [Times_n, States_BCR_n] = normalCR3BP_getTraj(time, rad2_n, u, lat0, lon0, az, el, v_mag_n);
%         plot3(States_BCR_n(:,1),States_BCR_n(:,2), States_BCR_n(:,3),'b')
        latlon_2 = zeros(length(Times_n),2);
        States_2_n = [States_BCR_n(:,1:3)-rB2_BCR_n, States_BCR_n(:,4:6)];
        for k = 1:length(Times_n)
            [lat, lon] = ECEF2latlon(States_2_n(k,1:3));
            lat = lat*180/pi;
            lon = lon*180/pi;
            if lat < latMin
                latMin = lat;
            end
            if lat > latMax
                latMax = lat;
            end
            if lon < lonMin
                lonMin = lon;
            end
            if lon > lonMax
                lonMax = lon;
            end
            latlon_2(k,:) = [lat, lon];
        end
        plot(latlon_2(:,2),latlon_2(:,1),'b','linewidth',1.5)
%         if size(States_BCR_n,1) < Times_n
%             d = Times_n - size(States_BCR_n,1);
%             States_BCR_n = [States_BCR_n; zeros(d,6)];
%         end
%         States_2_n_MULT(:,:,el_i,az_i) = States_BCR_n;
    end
end
% xlim([-180, 180]);
% ylim([-90, 90]);
xlim([lonMin, lonMax]);
ylim([latMin, latMax]);
grid on
axis equal
PlotBoi2('Longitude, \circ','Latitude, \circ',16)
989
return
%%% Calculating dimensional times
Times = Times_n.*tNorm;

%%% Calculating Body-2-Centered Normalized States
States_2_n = [States_BCR_n(:,1:3)-rB2_BCR_n, States_BCR_n(:,4:6)];

% ------------------------------------------------------------------------
%%% Post Processing
% ------------------------------------------------------------------------
latlonalt_2 = zeros(length(Times),3);
for k = 1:length(Times)
    [lat, lon] = ECEF2latlon(States_2_n(k,1:3));
    lat = lat*180/pi;
    lon = lon*180/pi;
    alt_n = (norm(States_2_n(k,1:3)) - rad2_n);
    latlonalt_2(k,:) = [lat, lon, alt_n];
end
clear lat lon alt

% figure; hold all
plot(latlonalt_2(:,2),latlonalt_2(:,1))
xlim([-180, 180]);
ylim([-90, 90]);
grid on
axis equal








