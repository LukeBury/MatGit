clear
clc
addpath('../../bin')

% ------------------------------------------------------------------------
%%% P4 Givens
% ------------------------------------------------------------------------
alt = 300; % km

u_E = 398600.4415;   % km^3/s^2, Earth
u_S = 1.32712428e11; % km^3/s^2, Sun
u_Sat = 3.794e7; % km^3/s^2, Saturn

a_E = 149598023; % km, Earth semimajor axis
a_Sat = 1429394133; % km, Saturn semimajor axis

r_E = 6378.1363; % km, Earth radius
r_Sat = 60268; % km, Saturn radius
% ------------------------------------------------------------------------
%%% 4a
% ------------------------------------------------------------------------
fprintf('---------------- 4a ----------------\n')
r = alt + r_E; % km

at = (a_E + a_Sat)/2; % km, transfer orbit semimajor axis

vcE = sqrt(u_S/a_E); % km/s, Circular Earth heliocentric velocity
vt1 = sqrt(2*u_S/a_E - u_S/at); % km/s, transfer orbit peri
vt2 = sqrt(2*u_S/a_Sat - u_S/at); % km/s, transfer orbit apo
vcSat = sqrt(u_S/a_Sat); % km/s, Circular Saturn helio velocity

vinf = vcSat - vt2; % km/s, vinf_minus approaching saturn

vsc_p = sqrt(2*u_S/a_Sat) % km/s, heliocentri velocity to escape solar sys

%%% Numerically finding turning angle
for d = 0:.01:180
    vec1 = [-vcSat, 0, 0];
    vec2 = [vinf*cosd(d), vinf*sind(d), 0];
    
    if abs(norm(vec1 + vec2) - vsc_p) < .0001
        theta = d % degrees
    end
end
theta = theta*pi/180; % rad

rp = (u_Sat/(vinf^2))*(1/(cos((pi-theta)/2)) - 1); % km, flyby radius

clearance = rp - r_Sat % km, distance from Saturn surface during flyby

% ------------------------------------------------------------------------
%%% 4b
% ------------------------------------------------------------------------
fprintf('---------------- 4b ----------------\n')
%%% Determining dV at Earth
vinf_E = vt1 - vcE; % km/s, vinf leaving Earth
vE2 = sqrt(vinf_E^2 + 2*u_E/r); % km/s, v required to achieve vinf
vE1 = sqrt(u_E/r); % km/s, intial velocity of circular orbit about Earth
dvE = vE2 - vE1 % km/s, dV required at Earth

%%% Determining unpowered dV from Saturn
dvSat = vsc_p - vt2 % km/s, unpowered dV provided by Saturn GA

% ------------------------------------------------------------------------
%%% 4c
% ------------------------------------------------------------------------
fprintf('---------------- 4c ----------------\n')
at

et = (a_Sat - a_E) / (a_Sat + a_E) % transfer orbit eccentricity

T_hohmann = pi * sqrt((at^3)/u_S); % sec

T_hohmann = T_hohmann / (3600*24*365.25) % years













