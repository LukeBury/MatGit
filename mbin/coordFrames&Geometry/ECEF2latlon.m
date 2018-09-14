function [lat, lon] = ECEF2latlon(r,unit,stupidMoonFlag)
%%% Given position vector, returns latitdue and longitude in RADIANS
%%% Inputs
%         (1) r - [nx3] body-centered-rotating position coordinates for s/c
%         (2) unit - [str] 'degrees' or 'radians' for outputs
%         (3) stupidMoonFlag - if this string == 'stupidMoon' then 
%         longitude will be corrected for secondary body in 3B system
%%% Outputs
%         (1) lat - [nx1] latitudes
%         (2) lon - [nx1] longitudes
% ========================================================================
%%% Dimension check
if size(r,2) ~= 3
    warning('Incorrect input dimensions')
    return
end

%%% Preallocate outputs
n = size(r,1);
lat = zeros(n,1);
lon = zeros(n,1);

%%% Define axes
x = [1 0 0]';
z = [0 0 1];

%%% Compute latitudes and longitudes
for kk = 1:n
    %%% Latitude
    lat(kk) = pi/2 - acos(dot(z,r(kk,:))/norm(r(kk,:)));

    %%% Longitude
    lon(kk) = atan(r(kk,2)/r(kk,1));

    %%% Quadrant corrections
    % -x, +y
    if r(kk,1) < 0 && r(kk,2) >= 0
        lon(kk) = lon(kk) + pi;
    % -x, -y
    elseif r(kk,1) < 0 && r(kk,2) < 0
        lon(kk) = lon(kk) - pi;
    end
    if nargin == 3 && isequal(stupidMoonFlag,'stupidMoon')
        lon(kk) = lon(kk) + pi;
        if lon(kk) > pi
            lon(kk) = lon(kk) - 2*pi;
        end
    end
end

%%% Format units
if strcmp(unit,'degrees')
    lat = lat .* 180/pi;
    lon = lon .* 180/pi;
elseif strcmp(unit,'radians')
    
else
    warning('Incorrect unit input')
    return
end

end