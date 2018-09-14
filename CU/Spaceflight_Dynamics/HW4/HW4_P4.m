clear
clc
% Givens
boulder_lat = 40.01 * pi / 180; % rads
boulder_lon = 254.83 * pi / 180; % rads
boulder_alt = 1.615; % km
r_ECEF      = [-1681; -5173; 4405]; % km

% Compute new values
[az_el_p] = ECEF2az_el_p(r_ECEF, boulder_lat, boulder_lon, boulder_alt);

% Print results
fprintf('Azimuth   = %3.5f deg\n', az_el_p(1) * 180 / pi)
fprintf('Elevation = %3.5f deg\n', az_el_p(2) * 180 / pi)
fprintf('Range     = %3.5f km\n', az_el_p(3))