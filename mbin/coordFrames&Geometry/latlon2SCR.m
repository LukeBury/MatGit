%%% Given latitude (deg), longitude (deg), and body radius (km)
%%% Returns surface coordinates in SCR frame where 0-lon on secondary body
%%% points towards the primary (-x)
%%%
%%% For use in CR3BP
function [rECEF] = latlon2SCR(lat, lon, rad)

%%% Convert to radians
lat = lat * pi/180;
lon = lon * pi/180 + pi;

%%% Assign x coordinate
x = rad * cos(lat) * cos(lon);

%%% Assign y coordinate
y = rad * cos(lat) * sin(lon);

%%% Assign z coordinate
z = sin(lat) * rad;

rECEF = [x, y, z];
end