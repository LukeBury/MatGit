function [dX] = Int_2BI_STM(t, X, mu)
%%% Description
%       For use with ODE integrators
%       
% --------------------------------------------------------------
%%% Inputs
%       t -  [1x1] time (s) 
%       X -  [6x1] State (km; km/s)
%       mu - [1x1] gravitational parameter of primary (km^3/s^2)
% --------------------------------------------------------------
%%% Outputs
%       dX - [6x1] time derivative of state (km/s; km/s^2)
% --------------------------------------------------------------
%%% Author
%       Luke Bury, luke.bury@colorado.edu
% ===============================================================
% ------------------------------------
%%% To get A and EQM
% ------------------------------------
% syms x y z dx dy dz mu
% r = sqrt(x^2 + y^2 + z^2);
% Utot_uJ2J3 = -mu/r; % Total potential
% EQM = [dx; dy; dz; diff(Utot_uJ2J3, x); diff(Utot_uJ2J3, y); diff(Utot_uJ2J3,z)];
% state = [x; y; z; dx; dy; dz];
% Asym = jacobian(EQM, state);
% 
% ------------------------------------
%%% Equations of motion
% ------------------------------------

%%% Size dY to fit state and all of reshaped (n^2,1) STM
dX = zeros(6+6^2,1);

%%% Unpack state
x = X(1);
y = X(2);
z = X(3);
dx = X(4);
dy = X(5);
dz = X(6);

%%% Reshape (n^2,1) stm to (n,n)
stm = reshape(X(7:end),6,6);

%%% Build A matrix and evaluate at current state
A = zeros(6,6);
A(1:3,4:6) = eye(3,3);
A(4,1) = -(mu/(x^2 + y^2 + z^2)^(3/2) - (3*mu*x^2)/(x^2 + y^2 + z^2)^(5/2));
A(4,2) = -(-(3*mu*x*y)/(x^2 + y^2 + z^2)^(5/2));
A(4,3) = -(-(3*mu*x*z)/(x^2 + y^2 + z^2)^(5/2));
A(5,1) = -(-(3*mu*x*y)/(x^2 + y^2 + z^2)^(5/2));
A(5,2) = -(mu/(x^2 + y^2 + z^2)^(3/2) - (3*mu*y^2)/(x^2 + y^2 + z^2)^(5/2));
A(5,3) = -(-(3*mu*y*z)/(x^2 + y^2 + z^2)^(5/2));
A(6,1) = -(-(3*mu*x*z)/(x^2 + y^2 + z^2)^(5/2));
A(6,2) = -(-(3*mu*y*z)/(x^2 + y^2 + z^2)^(5/2));
A(6,3) = -(mu/(x^2 + y^2 + z^2)^(3/2) - (3*mu*z^2)/(x^2 + y^2 + z^2)^(5/2));

%%% Calculate new STM
stm_dot = A*stm;

%%% Creating X-dot
% Spacecraft velocities
dX(1:3) = [dx; dy; dz];

%%% Using u, J2, and J3 terms as dynamics
% -EQM(4)
dX(4) = -(mu*x)/(x^2 + y^2 + z^2)^(3/2);
% -EQM(5)
dX(5) = -(mu*y)/(x^2 + y^2 + z^2)^(3/2);
% -EQM(6)
dX(6) = -(mu*z)/(x^2 + y^2 + z^2)^(3/2);

% Filling in reshaped (6^2,1) STM to state
dX(7:end) = reshape(stm_dot,36,1);


end