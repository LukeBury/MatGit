function [dX] = Int_2BI_J2(t,X,u,R,J2)
%%% For numerical integration of ECI state under influence of standard
%%% 2-body gravity and J2
%%% Inputs:
%          t - time vector
%          X - initial state [6x1]
%          u - gravitational parameter of primary body
%          R1 - radius of primary body
%          J2 - J2 of primary body

%%% Preallocate state output
dX = zeros(6,1);

%%% Unpack the barycentric state vector
x = X(1); y = X(2); z = X(3); % Position
dx = X(4); dy = X(5); dz = X(6); % Velocity

%%% Distances to primary (1) and secondary (2) bodies
r = sqrt(x^2 + y^2 + z^2);

%%% Define Gamma values for shortening J2 terms
gamma1 = 3*u*R^2*J2*(5*z^2-r^2)/(2*r^7);
gamma1_z = 3*u*R^2*J2*(5*z^2-3*r^2)/(2*r^7);

%%% Equations of Motion
ddx = x*(-u/r^3 + gamma1);
ddy = y*(-u/r^3 + gamma1);
ddz = z*(-u/r^3 + gamma1_z);

%%% Output the derivative of the state
dX(1:3) = [dx; dy; dz]; % km/s
dX(4:6) = [ddx; ddy; ddz]; % km/s^2

end
