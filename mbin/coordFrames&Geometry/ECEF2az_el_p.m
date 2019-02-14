% Compute the range (km), elevation (rads), and azimuth (rads) of a s/c
% given:
    % 1) r_ECEF of satellite (km)
    % 2-4) Tracking station latitude (rads), longitude (rads), and 
    %      altitude of station (km)
% Output:
    % az_el_p = [azimuth; elevation; range] (rads, rads, km)
function [az, el, p] = ECEF2az_el_p(r_ECEF, lat, lon, alt)
if size(r_ECEF) ~= [3 1]
    warning('r_ECEF must be 3x1')
end
% Computing ECEF position of tracking site
rE = 6378.1363; % km
r = alt + rE; % km
site_x = r * cos(lat) * cos(lon); % km
site_y = r * cos(lat) * sin(lon); % km
site_z = r * sin(lat); % km
site_ECEF = [site_x; site_y; site_z]; % km

% Computing slant-range vector
p_ECEF = r_ECEF - site_ECEF; % km

% Creating rotation matrices
t2 = pi/2 - lat; % theta R2, rads
t3 = lon; % theta R3, rads

R2 = [cos(t2) 0 -sin(t2); 0 1 0; sin(t2) 0 cos(t2)];
R3 = [cos(t3) sin(t3) 0;-sin(t3) cos(t3) 0; 0 0 1];

% Rotating p_ECEF to p_SEZ
p_SEZ = R2 * R3 * p_ECEF;

% Calculating azimuth, elevation, and range
p = norm(p_ECEF); % range, km
el = atan2(p_SEZ(3), sqrt(p_SEZ(1)^2 + p_SEZ(2)^2)); % rads
% warning('Double check this elevation equation sometime')
az = atan2(p_SEZ(2),-p_SEZ(1)); % rads
end