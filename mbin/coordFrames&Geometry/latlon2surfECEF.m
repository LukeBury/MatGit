%%% Given latitude (deg), longitude (deg), and body radius (km)
%%% Returns surface coordinates in ECEF
%%% 
%%% supply 'extras.stupidMoon' to set 0-longitude to -x axis
function [rECEF] = latlon2surfECEF(lat, lon, rad, extras)

%%% Convert to radians
lat = lat * pi/180;
lon = lon * pi/180;

if nargin == 4
    if isfield(extras,'stupidMoon') == 1
        lon = lon + pi;
    end
end

%%% Assign x coordinate
x = rad * cos(lat) * cos(lon);

%%% Assign y coordinate
y = rad * cos(lat) * sin(lon);

%%% Assign z coordinate
z = sin(lat) * rad;

rECEF = [x, y, z];
end
