function [lat_deg, lon_deg] = BCR2latlon(r_BCR, body, MR)
%%% Given position vector, returns latitdue and longitude in RADIANS
%%% Inputs
%         (1) r    - [3x1] body-centered-rotating position coordinates for s/c
%         (2) body - 'primary' or 'secondary'
%         (3) MR   - mass ratio of 3-body system
%%% Outputs
%         (1) lat_deg  - latitude (deg)
%         (2) lon_deg  - longitude (deg)
% ========================================================================
%%% Dimension check
if size(r_BCR) == [1,3]
    r_BCR = r_BCR';
end
if size(r_BCR) ~= [3,1]
    warning('Incorrect input dimensions')
    return
end

if isequal(body,'primary')
    r_BodyCR = r_BCR - [-MR; 0; 0];
elseif isequal(body,'secondary')
    r_BodyCR = r_BCR - [1 - MR; 0; 0];
end


%%% Define axes
x = [1 0 0]';
z = [0 0 1];

%%% Compute latitudes and longitudes
%%% Latitude
lat = pi/2 - acos(dot(z,r_BodyCR)/norm(r_BodyCR));

%%% Longitude
lon = atan(r_BodyCR(2)/r_BodyCR(1));

%%% Quadrant corrections
% -x, +y
if r_BodyCR(1) < 0 && r_BodyCR(2) >= 0
    lon = lon + pi;
% -x, -y
elseif r_BodyCR(1) < 0 && r_BodyCR(2) < 0
    lon = lon - pi;
end

%%% Moon convention
if isequal(body,'secondary')
    lon = lon + pi;
    if lon > pi
        lon = lon - 2*pi;
    end
end

%%% Format units
lat_deg = lat .* 180/pi;
lon_deg = lon .* 180/pi;

end