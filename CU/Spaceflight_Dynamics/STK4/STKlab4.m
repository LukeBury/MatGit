clear
clc

a_E = 147000000; % km
a_V = 108208601; % km

a_t = (a_E + a_V)/2 % km

u_S = 1.32712428e11;

vc_E = sqrt(u_S/a_E); % km/s
vc_V = sqrt(u_S/a_V); % km/s

vt1 = sqrt(2*u_S/a_E - u_S/a_t); % km/s
vt2 = sqrt(2*u_S/a_V - u_S/a_t); % km/s

dv1 = vc_E - vt1;
dv2 = vt2 - vc_V;