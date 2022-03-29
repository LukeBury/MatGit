function [dX] = Int_CR3Bn_J2(t,X,prms)
%%% For numerical integration in the normalized CR3BP with J2 of each body
%%% Inputs:
%          t - normalized time vector
%          X - initial state [6x1]
%          prms - (u, R1_n, R2_n, J21, J22)
% =======================================================================
%%% Preallocate state output
dX = zeros(6,1);

%%% Unpack the barycentric state vector
x = X(1); y = X(2); z = X(3); % Position
dx = X(4); dy = X(5); dz = X(6); % Velocity

%%% Define position of bodies
x1 = -prms.u;
x2 = 1-prms.u;

%%% Distances to primary (1) and secondary (2) bodies
r1 = sqrt((x+prms.u)^2 + y^2 + z^2);
r2 = sqrt((x+prms.u-1)^2 + y^2 + z^2);

%%% Define Gamma values for shortening J2 terms
gamma1 = 3*(1-prms.u)*prms.R1_n^2*prms.J21*(5*z^2-r1^2)/(2*r1^7);
gamma2 = 3*prms.u*prms.R2_n^2*prms.J22*(5*z^2-r2^2)/(2*r2^7);
gamma1_z = 3*(1-prms.u)*prms.R1_n^2*prms.J21*(5*z^2-3*r1^2)/(2*r1^7);
gamma2_z = 3*prms.u*prms.R2_n^2*prms.J22*(5*z^2-3*r2^2)/(2*r2^7);

%%% Equations of Motion
ddx = 2*prms.n*dy + (prms.n^2)*x + ((1-prms.u)/r1^3-gamma1)*(x1-x) + (prms.u/r2^3 - gamma2)*(x2-x);
ddy = -2*prms.n*dx + y*(-(1-prms.u)/r1^3 - prms.u/r2^3 + gamma1 + gamma2 + prms.n^2);
ddz = z*(-(1-prms.u)/r1^3 - prms.u/r2^3 + gamma1_z + gamma2_z);

%%% Output the derivative of the state
dX(1:3) = [dx; dy; dz]; % km/s
dX(4:6) = [ddx; ddy; ddz]; % km/s^2

end
