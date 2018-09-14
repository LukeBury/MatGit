clear
clc
addpath('../../bin')

% ------------------------------------------------------------------------
%%% P1 Givens
% ------------------------------------------------------------------------
u_E = 398600.4415;   % km^3/s^2, Earth
u_S = 1.32712428e11; % km^3/s^2, Sun
a_E = 149598023; % km, Earth semimajor axis
r_E = 6378.1363; % km, Earth radius
alt = 185; % km

% ------------------------------------------------------------------------
%%% a
% ------------------------------------------------------------------------
fprintf('---------------- 1a ----------------\n')
r = r_E + alt; % km
v1 = sqrt(u_E/r); % km/s
dv = 7; % km/s

v = v1 + dv; % km/s, Earth-relative velocity after burn

vinf = sqrt(v^2 - 2*u_E/r); % km/s, vinf relative to Earth

vE = sqrt(u_S/a_E); % km/s, Earth circular velocity

vsc = vE + vinf; % km/s, heliocentric sc velocity

at = 1 / (2/a_E - (vsc^2)/u_S); % km, new orbit semimajor axis

rp = a_E; % km, new orbit perihelion
ra = 2*at - rp % km, new orbit apohelion

% ------------------------------------------------------------------------
%%% b
% ------------------------------------------------------------------------
fprintf('---------------- 1b ----------------\n')
vsc2 = vE - vinf; % km/s, heliocentric sc velocity
at = 1 / (2/a_E - (vsc2^2)/u_S) % km, new orbit semimajor axis

ra2 = a_E; % km, new orbit perihelion
rp2 = 2*at - ra2 % km, new orbit apohelion



% ------------------------------------------------------------------------
%%% c
% ------------------------------------------------------------------------
fprintf('---------------- 1c ----------------\n')

r2d = 180/pi;

i = atan(vinf/vE)
i*r2d
















