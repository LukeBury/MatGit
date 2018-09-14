clear
clc
addpath('../../bin')

% ------------------------------------------------------------------------
%%% P6 Givens
% ------------------------------------------------------------------------
alt = 300; % km
r_E = 6378.1363; % km, Earth radius
u_E = 398600.4415;   % km^3/s^2, Earth
a = alt + r_E

x0 = 0; y0 = 0; z0 = 0; % m
x0dot = 0; y0dot = 1; z0dot = 0; % m/s

% ------------------------------------------------------------------------
%%% 6a
% ------------------------------------------------------------------------
fprintf('---------------- 6a ----------------\n')
t5 = 5 * 60; % sec
w = sqrt(u_E / (a^3)); % rad/s
wt5 = w*t5; % rad

%%% Compute experiment position at t1
x5 = 2*y0dot*cos(wt5)/w + 2*y0dot/w % m
y5 = 4*y0dot*sin(wt5)/w - 3*y0dot*t5 % m
z5 = 0 % m

% ------------------------------------------------------------------------
%%% 6b
% ------------------------------------------------------------------------
fprintf('---------------- 6b ----------------\n')

%%% Compute experiment velocity at t1
x5dot = 2*y0dot*sin(wt5); % m/s
y5dot = 4*y0dot*cos(wt5) - 3*y0dot; % m/s
z5dot = 0; % m/s

%%% Compute necessary experiment velocity to rendezvous at t2 (from t1
%%% position)
t10 = 10 * 60; % sec
wt10 = w*t10; % rad, relative to t1

ydot_r = ((6*x5*(wt10 - sin(wt10)) - y5)*w*sin(wt10) - 2*w*x5*(4-3*cos(wt10))*(1-cos(wt10)))/...
    ((4*sin(wt10)-3*wt10)*sin(wt10) + 4*((1-cos(wt10))^2));
xdot_r = -(w*x5*(4-3*cos(wt10)) + 2*(1-cos(wt10))*ydot_r)/sin(wt10);
zdot_r = 0; % m/s

%%% dV for experiment
dV_r_e = [xdot_r-x5dot; ydot_r-y5dot; zdot_r-z5dot]
%%% dV for shuttle
dV_r_s = -dV_r_e
dV_r = norm(dV_r_s)

% ------------------------------------------------------------------------
%%% 6c
% ------------------------------------------------------------------------
fprintf('---------------- 6c ----------------\n')

%%% Final velocity of experiment at rendezvous
xdot_f = xdot_r*cos(wt10) + (3*w*x5 + 2*ydot_r)*sin(wt10)
ydot_f = (6*w*x5 + 4*ydot_r)*cos(wt10) - 2*xdot_r*sin(wt10) - (6*w*x5 + 3*ydot_r)

%%% above vector is also the appropriate thrust vector for the shuttle
norm([xdot_f ydot_f]) % m/s


























