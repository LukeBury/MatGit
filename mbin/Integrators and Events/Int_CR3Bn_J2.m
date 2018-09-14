function [dX] = Int_CR3Bn_J2(t,X,u,R1,R2,J21,J22)
%%% For numerical integration in the normalized CR3BP with J2 of each body
%%% Inputs:
%          t - normalized time vector
%          X - initial state [6x1]
%          u - mass ratio of CR3BP system
%          R1 - radius of primary body
%          R2 - radius of secondary body
%          J21 - J2 of primary body
%          J22 - J2 of secondary body
% =======================================================================
%%% Preallocate state output
dX = zeros(6,1);

%%% Unpack the barycentric state vector
x = X(1); y = X(2); z = X(3); % Position
dx = X(4); dy = X(5); dz = X(6); % Velocity

%%% Define position of bodies
x1 = -u;
x2 = 1-u;

%%% Distances to primary (1) and secondary (2) bodies
r1 = sqrt((x+u)^2 + y^2 + z^2);
r2 = sqrt((x+u-1)^2 + y^2 + z^2);

%%% Define Gamma values for shortening J2 terms
gamma1 = 3*(1-u)*R1^2*J21*(5*z^2-r1^2)/(2*r1^7);
gamma2 = 3*u*R2^2*J22*(5*z^2-r2^2)/(2*r2^7);
gamma1_z = 3*(1-u)*R1^2*J21*(5*z^2-3*r1^2)/(2*r1^7);
gamma2_z = 3*u*R2^2*J22*(5*z^2-3*r2^2)/(2*r2^7);

%%% Equations of Motion
ddx = 2*dy + x + ((1-u)/r1^3-gamma1)*(x1-x) + (u/r2^3 - gamma2)*(x2-x);
ddy = -2*dx + y*(-(1-u)/r1^3 - u/r2^3 + gamma1 + gamma2 + 1);
ddz = z*(-(1-u)/r1^3 - u/r2^3 + gamma1_z + gamma2_z);

%%% Output the derivative of the state
dX(1:3) = [dx; dy; dz]; % km/s
dX(4:6) = [ddx; ddy; ddz]; % km/s^2

end
