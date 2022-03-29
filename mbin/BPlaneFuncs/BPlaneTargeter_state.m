function [BT, BR, B, theta_rad, S_hat, T_hat, R_hat] = BPlaneTargeter_state(X, mu)
%%% Description
%       Calculates B-Plane targeting parameters from initial state
%     
% --------------------------------------------------------------
%%% Inputs
%       X  - Initial state vector (WRT FLYBY BODY) (km, km/s) [6x1]
%       mu - Gravitational parameter of central body (km^3/s^2)
% --------------------------------------------------------------
%%% Outputs
%       BT        - [3x1] B-dot-T vector (km)
%       BR        - [3x1] B-dot-R vector (km)
%       B         - [3x1] B vector (km)
%       theta_rad - [1x1] Theta vector (angle of B vector from B_t), (rad)
%       S_hat     - [3x1] S_hat vector (km)
%       R_hat     - [3x1] R_hat vector (km)
%       T_hat     - [3x1] T_hat vector (km)
% --------------------------------------------------------------
%%% Author
%       Luke Bury, luke.bury@colorado.edu
% ===============================================================

%%% Breaking down state
r = X(1:3);
v = X(4:6);

%%% State magnitudes
R = norm(r)
V = norm(v)

%%% Necessary parameters
h_hat = cross(r,v)/norm(cross(r,v))
e_vec = (1/mu)*((V^2 - mu/R)*r - dot(r,v)*v)
e = norm(e_vec)
p = acos(1/e)

%%% Defining axes
S_hat = cos(p)*(e_vec/e) + sin(p)*cross(h_hat,e_vec)/norm(cross(h_hat,e_vec));
k_hat = [0; 0; 1];
T_hat = cross(S_hat, k_hat)/norm(cross(S_hat,k_hat));
R_hat = cross(S_hat, T_hat);
B_hat = cross(S_hat, h_hat);

%%% |Semimajor Axis|
a = abs(1 / (2/R - (V^2)/mu));

%%% 'b' vector
b = a*sqrt(e^2 - 1);

%%% 'B' vector
B = B_hat.*b;

%%% B-dot-T and B-dot-R
BT = dot(B,T_hat);
BR = dot(B,R_hat);

%%% Angle from B-dot-T to B
theta_rad = acos(dot(T_hat,B_hat));

end



