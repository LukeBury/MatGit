clear
clc
% ------------------------------------------------------------------------
%%% Givens
% ------------------------------------------------------------------------
u = 398600; % km^3 / s^2
Re = 6378.1363; % km
rp1 = 250 + Re; % km
ra1 = 600 + Re; % km
rp2 = 2000 + Re; % km
ra2 = 5000 + Re; % km

% ------------------------------------------------------------------------
%%% 1a
% ------------------------------------------------------------------------
% determining semimajor axes for initial and final orbits
a1 = (rp1 + ra1) / 2; % km
a2 = (rp2 + ra2) / 2; % km

% determining periapse and apoapse velocities for initial and final orbits
vp1 = sqrt(2*u/rp1 - u/a1); % km/s
va1 = sqrt(2*u/ra1 - u/a1); % km/s
vp2 = sqrt(2*u/rp2 - u/a2); % km/s
va2 = sqrt(2*u/ra2 - u/a2); % km/s

% determining rp, ra, and a for transfer orbit
rpt = rp1; % km
rat = ra2; % km
at  = (rpt + rat)/2; % km

% determining periapse and apoapse velocities for transfer orbit
vpt = sqrt(2*u/rpt - u/at); % km/s
vat = sqrt(2*u/rat - u/at); % km/s

dv1_a = vpt - vp1; % km/s
dv2_a = va2 - vat; % km/s
dv_total_a = dv1_a + dv2_a; % km/s

fprintf('dv1_a      = %1.3f m/s\n', dv1_a * 1000)
fprintf('dv2_a      = %1.3f m/s\n', dv2_a * 1000)
fprintf('dv_total_a = %1.3f m/s\n', dv_total_a * 1000)
fprintf('--------------------------\n')
% ------------------------------------------------------------------------
%%% 1b
% ------------------------------------------------------------------------
% determining rp, ra, and a for transfer orbit
rpt = rp1; % km
rat = rp2; % km
at  = (rpt + rat)/2; % km

% determining periapse and apoapse velocities for transfer orbit
vpt = sqrt(2*u/rpt - u/at); % km/s
vat = sqrt(2*u/rat - u/at); % km/s

dv1_b = vpt - vp1; % km/s
dv2_b = vp2 - vat; % km/s
dv_total_b = dv1_b + dv2_b; % km/s

fprintf('dv1_b      = %1.3f m/s\n', dv1_b * 1000)
fprintf('dv2_b      = %1.3f m/s\n', dv2_b * 1000)
fprintf('dv_total_b = %1.3f m/s\n', dv_total_b * 1000)
fprintf('--------------------------\n')
% ------------------------------------------------------------------------
%%% 1c
% ------------------------------------------------------------------------
% determining rp, ra, and a for transfer orbit
rpt = ra1; % km
rat = ra2; % km
at  = (rpt + rat)/2; % km

% determining periapse and apoapse velocities for transfer orbit
vpt = sqrt(2*u/rpt - u/at); % km/s
vat = sqrt(2*u/rat - u/at); % km/s

dv1_c = vpt - va1; % km/s
dv2_c = va2 - vat; % km/s
dv_total_c = dv1_c + dv2_c; % km/s

fprintf('dv1_c      = %1.3f m/s\n', dv1_c * 1000)
fprintf('dv2_c      = %1.3f m/s\n', dv2_c * 1000)
fprintf('dv_total_c = %1.3f m/s\n', dv_total_c * 1000)
fprintf('--------------------------\n')
% ------------------------------------------------------------------------
%%% 1d
% ------------------------------------------------------------------------
% determining rp, ra, and a for transfer orbit
rpt = ra1; % km
rat = rp2; % km
at  = (rpt + rat)/2; % km

% determining periapse and apoapse velocities for transfer orbit
vpt = sqrt(2*u/rpt - u/at); % km/s
vat = sqrt(2*u/rat - u/at); % km/s

dv1_d = vpt - va1; % km/s
dv2_d = vp2 - vat; % km/s
dv_total_d = dv1_d + dv2_d; % km/s

fprintf('dv1_d      = %1.3f m/s\n', dv1_d * 1000)
fprintf('dv2_d      = %1.3f m/s\n', dv2_d * 1000)
fprintf('dv_total_d = %1.3f m/s\n', dv_total_d * 1000)
fprintf('--------------------------\n')
% ------------------------------------------------------------------------
%%% 1e
% ------------------------------------------------------------------------
fprintf('periapse to apoapse (a) had the least dv')

