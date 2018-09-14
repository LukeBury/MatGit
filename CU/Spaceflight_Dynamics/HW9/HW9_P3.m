clear
clc

% ------------------------------------------------------------------------
%%% Setup
% ------------------------------------------------------------------------
%%% Gravitational Parameters
u_E = 398600.4415;   % km^3/s^2, Earth
u_S = 1.32712428e11; % km^3/s^2, Sun
u_J = 1.268e8; % km^3/s^2, Jupiter

%%% Semimajor Axis of Circular Orbits
a_E = 149598023; % km, Earth
a_J = 778298361; % km, Jupiter

% ------------------------------------------------------------------------
%%% Hohmann Transfer
% ------------------------------------------------------------------------
a_t = (a_E + a_J)/2; % km, Transfer orbit

%%% General Hohmann transfer numbers
vt2 = sqrt(2*u_S/a_J - u_S/a_t); % km/s, Transfer orbit apoapse speed
vcJ = sqrt(u_S/a_J); % km/s, Circular Jupiter heliocentric velocity

%%% vInf at arrival
vInf = vcJ - vt2; % km/s, Mars-arrival vInf magnitude

% ------------------------------------------------------------------------
%%% Flyby
% ------------------------------------------------------------------------
%%% Jupiter flyby dimensions
r_J = 71492; % km
alt = 71492; % km
r_fb = r_J + alt; % km

%%% Flyby Turning angle
d_fb = pi - 2*acos(1/((r_fb*(vInf^2)/u_J) + 1)); % rad

%%% Creating velocity vectors
vInf_minus = [vInf 0]; % km/s
vJS = [-vcJ 0]; % km/s, Jupiter velocity wrt Sun
vSC_minus = vt2;

%%% Post flyby
vInf_plus = [cos(d_fb)*norm(vInf_minus) sin(d_fb)*norm(vInf_minus)]
vSC_plus = vInf_plus + vJS % km/s, SC departing velocity wrt Sun
norm(vSC_plus) % km/s, SC departing velocity wrt Sun

%%% Energy
E = (norm(vSC_plus)^2)/2 - u_S/a_J; % km^2/s^2
E_jkg = E*1000000
fprintf('Not possible b/c total energy is negative\n')



