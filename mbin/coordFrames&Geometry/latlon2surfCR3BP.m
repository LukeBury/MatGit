function [rBCR] = latlon2surfCR3BP(body, lat, lon, MR, radius)
%%% Given latitude and longitude, gives coordinates on surface of body in
%%% CR3BP
%%% 
%%% Inputs:
%          1) body   - 'primary' or 'secondary'
%          2) lat    - latitude (deg)
%          3) lon    - longitude (deg)
%          4) MR     - mass ratio of 3-body system
%          5) radius - radius of body
%%% Outputs:
%          1) rBCR   - requested position vector in BCR coordinates [3x1]
% ========================================================================

%%% Convert to radians
lat = lat * pi/180;
lon = lon * pi/180;

if isequal(body,'secondary')
    lon = lon + pi;
end

%%% Assign x coordinate
x = radius * cos(lat) * cos(lon);

%%% Assign y coordinate
y = radius * cos(lat) * sin(lon);

%%% Assign z coordinate
z = sin(lat) * radius;

rBodyCR = [x; y; z];


if isequal(body,'primary')
    rBCR = rBodyCR + [-MR; 0; 0];
elseif isequal(body,'secondary')
    rBCR = rBodyCR + [1 - MR; 0; 0];
end

end
