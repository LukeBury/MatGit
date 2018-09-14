clear
clc
addpath('../../bin')

% ------------------------------------------------------------------------
%%% P5 Givens
% ------------------------------------------------------------------------
m_ast = 1e15; % kg
m_S = 1.9891e30; % kg, Sun
a_ast = 3; % AU
a_ast = a_ast * 149597870; % km

% ------------------------------------------------------------------------
%%% 5a
% ------------------------------------------------------------------------
fprintf('---------------- 5a ----------------\n')
rSOI = ((m_ast/m_S)^(2/5))*a_ast % km, asteroid SOI

% ------------------------------------------------------------------------
%%% 5b
% ------------------------------------------------------------------------
fprintf('---------------- 5b ----------------\n')
a_sat = 40; % km, sc semimajor axis about asteroid
G = 6.673e-20; % km^3/(kg*s^2)

u_ast = m_ast * G; % km^3/s^2

Tp_sat = 2*pi*sqrt((a_sat^3)/u_ast) % sec

% ------------------------------------------------------------------------
%%% 5c
% ------------------------------------------------------------------------
fprintf('---------------- 5c ----------------\n')
dV1 = 0.001; % km/s

v1c = sqrt(u_ast/a_sat) % km/s, sc initial circular-orbit velocity

v2a = v1c - dV1 % km/s, apo velocity of new elliptical orbit

%%% Finding eccentricity of new elliptical orbit
a_sat2 = 1/(2/a_sat - (v2a^2)/u_ast); % km
ra2 = a_sat; % km
rp2 = 2*a_sat2 - ra2;
e2 = (ra2 - rp2)/(ra2 + rp2);

%%% Calculating mean motion
n2 = sqrt(u_ast / (a_sat2^3)); % rad/s

%%% To find r2: M1->M2, M2->E2->ta2, ta2->r2
M1 = pi; % rad
time = 6*3600; % sec
M2 = M1 + n2*time; % rad
E2 = M2E(M2,e2); % rad
ta2 = E2T(E2,e2); % rad
r2 = a_sat2 * ((1 - e2^2)/(1 + e2*cos(ta2)))

% ------------------------------------------------------------------------
%%% 5d
% ------------------------------------------------------------------------
fprintf('---------------- 5d ----------------\n')

v2 = sqrt(2*u_ast/r2 - u_ast/a_sat2); % velocity at r2 (on elliptical)
a_sat3 = r2; % km, semimajor axis of 3rd orbit (circular)
v3 = sqrt(u_ast / a_sat3); % km/s, velocity for circular orbit

fpa = acos((1 + e2*cos(ta2))/sqrt(1 + e2^2 + 2*e2*cos(ta2)));

v3_vec = [0, v3, 0]'; % km/s, circ vel in RSW
v2_vec = [-sin(fpa)*v2, cos(fpa)*v2, 0]'; % km/s, v2 in RSW

dV2_vec = v3_vec - v2_vec; % km/s, dV in RSW
norm(dV2_vec)

% ------------------------------------------------------------------------
%%% 5e
% ------------------------------------------------------------------------
fprintf('---------------- 5e ----------------\n')






















