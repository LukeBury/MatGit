clear
clc
% ------------------------------------------------------------------------
%%% Givens
% ------------------------------------------------------------------------
u_Earth = 398600; % km^3 / s^2
Re = 6378.1363; % km
a_ISS = 400 + Re; % km
t = 30 * 60; % seconds

% ------------------------------------------------------------------------
%%% Define state of target (ISS) relative to interceptor (Progress)
% ------------------------------------------------------------------------
x0 = 200; % m
y0 = 300; % m
z0 = 50; % m

xdot0 = 0.1; % m/s
ydot0 = 0.0; % m/s
zdot0 = -0.1; % m/s

w = sqrt(u_Earth / (a_ISS^3)); % rad/s

% ------------------------------------------------------------------------
%%% Calculating dV1
% ------------------------------------------------------------------------
%%% Calculating initial velocity of interceptor on rendezvous trajectory
% Taking x(t) eqns and setting x(t) = 0, solving for xdot0)
wt = w*t; % rads
ydot_i = ((6*x0*(wt - sin(wt)) - y0)*w*sin(wt) - 2*w*x0*(4-3*cos(wt))*(1-cos(wt)))/...
    ((4*sin(wt) - 3*wt)*sin(wt) + 4*(1-cos(wt))^2);
xdot_i = -(w*x0*(4-3*cos(wt)) + 2*(1 - cos(wt))*ydot_i)/sin(wt);
zdot_i = -z0*w*cot(wt);
%%% Calculating dV1
dV1_x = abs(xdot_i - xdot0);
dV1_y = abs(ydot_i - ydot0);
dV1_z = abs(zdot_i - zdot0);

dV1 = norm([dV1_x, dV1_y, dV1_z])


% ------------------------------------------------------------------------
%%% Calculating dV2
% ------------------------------------------------------------------------
%%% Calculating target speed at t2 on rendezvous trajectory
% Using the ()_dot equations with initial velocities equal to the ()dot_i
% velocities calculated above
xdot_tgt = xdot_i*cos(w*t) + (3*w*x0 + 2*ydot_i)*sin(w*t);
ydot_tgt = (6*w*x0 + 4*ydot_i)*cos(w*t) - 2*xdot_i*sin(w*t) - (6*w*x0 + 3*ydot_i);
zdot_tgt = -z0*w*sin(w*t) + zdot_i*cos(w*t);

% x_t = (xdot0/w)*sin(w*t1) - (3*x0 + 2*ydot0/w)*cos(w*t1) + (4*x0 + 2*ydot0/w)
% y_t = (6*x0 + 4*ydot0/w)*sin(w*t1) + 2*xdot0/w*cos(w*t1) - (6*w*x0 + 3*ydot0)*t1 +...
%     (y0 - 2*xdot0/w)
% z_t = z0*cos(w*t1) + zdot0/w*sin(w*t1)
% 
% xdot_t = xdot0*cos(w*t1) + (3*w*x0 + 2*ydot0)*sin(w*t1)
% ydot_t = (6*w*x0 + 4*ydot0)*cos(w*t1) - 2*xdot0*sin(w*t1) - (6*w*x0 + 3*ydot0)
% zdot_t = -z0*w*sin(w*t1) + zdot0*cos(w*t1)

dV2 = norm([xdot_tgt, ydot_tgt, zdot_tgt])

dVtotal = dV1 + dV2

%ans=table(a_tgt, omega_trgt, a_phase, tao_trans);







