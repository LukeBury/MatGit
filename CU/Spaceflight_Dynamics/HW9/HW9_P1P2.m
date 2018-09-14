clear
clc
close all
% ------------------------------------------------------------------------
%%% P1 Givens
% ------------------------------------------------------------------------
u_E = 398600.4415;   % km^3/s^2, Earth
u_S = 1.32712428e11; % km^3/s^2, Sun
u_M = 4.305e4;       % km^3/s^2, Mars
r_E = 6378.1363; % km, Earth radius
r_M = 3397.2; % km, Mars radius
alt_E = 185; % km
alt_M = 300; % km

%%% Semimajor Axis of Circular Orbits
a_E = 149598023; % km, Earth
a_M = 227939186; % km, Mars

%%% Masses
m_E = 5.9742e24; % kg, Earth
m_M = 6.4191e23; % kg, Mars
m_S = 1.9891e30; % kg, Sun
m_Moon = 7.3483e22; % kg, Moon

% ------------------------------------------------------------------------
%%% SOI Calculations
% ------------------------------------------------------------------------
fprintf('------- a -------\n')
rSOI_E = ((m_E/m_S)^(2/5))*a_E % km, Earth
rSOI_M = ((m_M/m_S)^(2/5))*a_M % km, Mars
rSOI_EM = (((m_Moon+m_E)/m_S)^(2/5))*a_E % km, Earth-Moon

% ------------------------------------------------------------------------
%%% Earth Mars Hohmann Transfer
% ------------------------------------------------------------------------
fprintf('------- b -------\n')
a_t = (a_E + a_M)/2; % km, Transfer orbit

%%% General Hohmann transfer numbers
vcE = sqrt(u_S/a_E); % km/s, Circular Earth heliocentric velocity
vt1 = sqrt(2*u_S/a_E - u_S/a_t) % km/s, Transfer orbit periapse speed
vt2 = sqrt(2*u_S/a_M - u_S/a_t) % km/s, Transfer orbit apoapse speed
vcM = sqrt(u_S/a_M) % km/s, Circular Mars heliocentric velocity

fprintf('------- c -------\n')
%%% vInf at departure and arrival
vInf1 = vt1 - vcE % km/s, Earth-departure vInf magnitude
vInf2 = vcM - vt2 % km/s, Mars-arrival vInf magnitude

fprintf('------- d -------\n')
%%% dV1
vParkE = sqrt(u_E/(alt_E+r_E)); % km/s, Velocity in parking orbit
vNeedE = sqrt(vInf1^2 + 2*u_E/(alt_E+r_E)); % km/s, Vel to achieve vInf1
dV1 = vNeedE - vParkE % km/s

fprintf('------- e -------\n')
%%% dV2
vNeedM = sqrt(u_M/(alt_M+r_M)); % km/s, Velocity in circular Mars orbit
vArrM = sqrt(vInf2^2 + 2*u_M/(alt_M+r_M)); % km/s, Arrival velocity at Mars
dV2 = vArrM - vNeedM % km/s

fprintf('------- f -------\n')
%%% Transfer Time
T_trans = pi * sqrt(a_t^3 / u_S); % seconds
T_trans/86400 % days

% ------------------------------------------------------------------------
%%% P2 - flyby
% ------------------------------------------------------------------------
fprintf('------- 2 -------\n')
%%% Flyby radius
r_fb = 200 + r_M % km

%%% Max Turning angle
d_fb = pi - 2*acos(1/((r_fb*(vInf2^2)/u_M) + 1)) % rad

vInf_minus = [vInf2 0]; % km/s
vMS = [-vcM 0] % km/s, Mars velocity wrt Sun

%%% Fast flyby
vInf_plus = [cos(d_fb)*norm(vInf_minus) sin(d_fb)*norm(vInf_minus)]
vt2
vSC_plus = vInf_plus + vMS % km/s, SC departing velocity wrt Sun
norm(vSC_plus) % km/s, SC departing velocity wrt Sun







