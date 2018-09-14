clear
clc

%%% Study values
vs = 2130.552; % m/s
ms = 0.023; % kg

%%% Assume half velocity
v1 = vs/2; % m/s

%%% Necessary momentum
momNec = ms * v1; % kg*m/s

%%% Copper Density
pCopper = 8960; % kg/m^3

%%% Necessary radius for sphere of mass ... V = (4/3)*pi*r^3 = mass/dens
% r = ((mass*3) / (4*dens*pi))^(1/3)
r = ((ms*3) / (4*pCopper*pi))^(1/3);

%%% Area of Chamber
A = pi*r*r;

%%% Assume time in Chamber
t = 0.01; % sec

%%% Necessary pressure
% P = mom/(A*t)
P_Pa = momNec / (A*t) % Pa

%%% Pa to psi
P_psi = P_Pa * 0.000145038
