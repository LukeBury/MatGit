function [dX] = Int_CR3BnSTM_J2(t,X,prms)
%%% For numerical integration in the normalized CR3BP
%%% Inputs:
%          t - normalized time vector
%          X - initial state [42x1]
%          prms.u - mass ratio of CR3BP system
%          prms.R1 - radius of primary body
%          prms.R2 - radius of secondary body
%          prms.J2p - J2 of primary body
%          prms.J2s - J2 of secondary body
%          prms.n   - normalized mean motion
% =======================================================================
%%% Preallocate state output
dX = zeros(42,1);

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
gamma1 = 3*(1-prms.u)*prms.R1^2*prms.J2p*(5*z^2-r1^2)/(2*r1^7);
gamma2 = 3*prms.u*prms.R2^2*prms.J2s*(5*z^2-r2^2)/(2*r2^7);
gamma1_z = 3*(1-prms.u)*prms.R1^2*prms.J2p*(5*z^2-3*r1^2)/(2*r1^7);
gamma2_z = 3*prms.u*prms.R2^2*prms.J2s*(5*z^2-3*r2^2)/(2*r2^7);

%%% Equations of Motion
ddx = 2*prms.n*dy + (prms.n^2)*x + ((1-prms.u)/r1^3-gamma1)*(x1-x) + (prms.u/r2^3 - gamma2)*(x2-x);
ddy = -2*prms.n*dx + y*(-(1-prms.u)/r1^3 - prms.u/r2^3 + gamma1 + gamma2 + prms.n^2);
ddz = z*(-(1-prms.u)/r1^3 - prms.u/r2^3 + gamma1_z + gamma2_z);

%%% Output the derivative of the state
dX(1:3) = [dx; dy; dz]; % km/s
dX(4:6) = [ddx; ddy; ddz]; % km/s^2

%%% Reshaping STM
stm = reshape(X(7:end),6,6);

%%% Evaluate A matrix
% A = zeros(6,6);
% A(1:3,4:6) = eye(3);
% A(4,5) = 2;
% A(5,4) = -2;
% A(4,1) = (prms.u + x - 1)*((3*prms.u*(2*prms.u + 2*x - 2))/(2*((prms.u + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*prms.J2s*prms.R2^2*prms.u*(2*prms.u + 2*x - 2))/(2*((prms.u + x - 1)^2 + y^2 + z^2)^(7/2)) + (21*prms.J2s*prms.R2^2*prms.u*(2*prms.u + 2*x - 2)*((prms.u + x - 1)^2 + y^2 - 4*z^2))/(4*((prms.u + x - 1)^2 + y^2 + z^2)^(9/2))) + (prms.u - 1)/((prms.u + x)^2 + y^2 + z^2)^(3/2) - (prms.u + x)*((3*(2*prms.u + 2*x)*(prms.u - 1))/(2*((prms.u + x)^2 + y^2 + z^2)^(5/2)) - (prms.J2p*prms.R1^2*(2*prms.u + 2*x)*(3*prms.u - 3))/(2*((prms.u + x)^2 + y^2 + z^2)^(7/2)) + (7*prms.J2p*prms.R1^2*(2*prms.u + 2*x)*(3*prms.u - 3)*((prms.u + x)^2 + y^2 - 4*z^2))/(4*((prms.u + x)^2 + y^2 + z^2)^(9/2))) - prms.u/((prms.u + x - 1)^2 + y^2 + z^2)^(3/2) - (3*prms.J2s*prms.R2^2*prms.u*((prms.u + x - 1)^2 + y^2 - 4*z^2))/(2*((prms.u + x - 1)^2 + y^2 + z^2)^(7/2)) + (prms.J2p*prms.R1^2*(3*prms.u - 3)*((prms.u + x)^2 + y^2 - 4*z^2))/(2*((prms.u + x)^2 + y^2 + z^2)^(7/2)) + 1;
% A(4,2) = (prms.u + x - 1)*((3*prms.u*y)/((prms.u + x - 1)^2 + y^2 + z^2)^(5/2) - (3*prms.J2s*prms.R2^2*prms.u*y)/((prms.u + x - 1)^2 + y^2 + z^2)^(7/2) + (21*prms.J2s*prms.R2^2*prms.u*y*((prms.u + x - 1)^2 + y^2 - 4*z^2))/(2*((prms.u + x - 1)^2 + y^2 + z^2)^(9/2))) - (prms.u + x)*((3*y*(prms.u - 1))/((prms.u + x)^2 + y^2 + z^2)^(5/2) - (prms.J2p*prms.R1^2*y*(3*prms.u - 3))/((prms.u + x)^2 + y^2 + z^2)^(7/2) + (7*prms.J2p*prms.R1^2*y*(3*prms.u - 3)*((prms.u + x)^2 + y^2 - 4*z^2))/(2*((prms.u + x)^2 + y^2 + z^2)^(9/2)));
% A(4,3) = (prms.u + x - 1)*((3*prms.u*z)/((prms.u + x - 1)^2 + y^2 + z^2)^(5/2) + (12*prms.J2s*prms.R2^2*prms.u*z)/((prms.u + x - 1)^2 + y^2 + z^2)^(7/2) + (21*prms.J2s*prms.R2^2*prms.u*z*((prms.u + x - 1)^2 + y^2 - 4*z^2))/(2*((prms.u + x - 1)^2 + y^2 + z^2)^(9/2))) - (prms.u + x)*((3*z*(prms.u - 1))/((prms.u + x)^2 + y^2 + z^2)^(5/2) + (4*prms.J2p*prms.R1^2*z*(3*prms.u - 3))/((prms.u + x)^2 + y^2 + z^2)^(7/2) + (7*prms.J2p*prms.R1^2*z*(3*prms.u - 3)*((prms.u + x)^2 + y^2 - 4*z^2))/(2*((prms.u + x)^2 + y^2 + z^2)^(9/2)));
% A(5,1) = -y*((3*(2*prms.u + 2*x)*(prms.u - 1))/(2*((prms.u + x)^2 + y^2 + z^2)^(5/2)) - (3*prms.u*(2*prms.u + 2*x - 2))/(2*((prms.u + x - 1)^2 + y^2 + z^2)^(5/2)) - (prms.J2p*prms.R1^2*(2*prms.u + 2*x)*(3*prms.u - 3))/(2*((prms.u + x)^2 + y^2 + z^2)^(7/2)) + (3*prms.J2s*prms.R2^2*prms.u*(2*prms.u + 2*x - 2))/(2*((prms.u + x - 1)^2 + y^2 + z^2)^(7/2)) - (21*prms.J2s*prms.R2^2*prms.u*(2*prms.u + 2*x - 2)*((prms.u + x - 1)^2 + y^2 - 4*z^2))/(4*((prms.u + x - 1)^2 + y^2 + z^2)^(9/2)) + (7*prms.J2p*prms.R1^2*(2*prms.u + 2*x)*(3*prms.u - 3)*((prms.u + x)^2 + y^2 - 4*z^2))/(4*((prms.u + x)^2 + y^2 + z^2)^(9/2)));
% A(5,2) = (prms.u - 1)/((prms.u + x)^2 + y^2 + z^2)^(3/2) + y*((3*prms.u*y)/((prms.u + x - 1)^2 + y^2 + z^2)^(5/2) - (3*y*(prms.u - 1))/((prms.u + x)^2 + y^2 + z^2)^(5/2) + (prms.J2p*prms.R1^2*y*(3*prms.u - 3))/((prms.u + x)^2 + y^2 + z^2)^(7/2) - (3*prms.J2s*prms.R2^2*prms.u*y)/((prms.u + x - 1)^2 + y^2 + z^2)^(7/2) - (7*prms.J2p*prms.R1^2*y*(3*prms.u - 3)*((prms.u + x)^2 + y^2 - 4*z^2))/(2*((prms.u + x)^2 + y^2 + z^2)^(9/2)) + (21*prms.J2s*prms.R2^2*prms.u*y*((prms.u + x - 1)^2 + y^2 - 4*z^2))/(2*((prms.u + x - 1)^2 + y^2 + z^2)^(9/2))) - prms.u/((prms.u + x - 1)^2 + y^2 + z^2)^(3/2) - (3*prms.J2s*prms.R2^2*prms.u*((prms.u + x - 1)^2 + y^2 - 4*z^2))/(2*((prms.u + x - 1)^2 + y^2 + z^2)^(7/2)) + (prms.J2p*prms.R1^2*(3*prms.u - 3)*((prms.u + x)^2 + y^2 - 4*z^2))/(2*((prms.u + x)^2 + y^2 + z^2)^(7/2)) + 1;
% A(5,3) = y*((3*prms.u*z)/((prms.u + x - 1)^2 + y^2 + z^2)^(5/2) - (3*z*(prms.u - 1))/((prms.u + x)^2 + y^2 + z^2)^(5/2) - (4*prms.J2p*prms.R1^2*z*(3*prms.u - 3))/((prms.u + x)^2 + y^2 + z^2)^(7/2) + (12*prms.J2s*prms.R2^2*prms.u*z)/((prms.u + x - 1)^2 + y^2 + z^2)^(7/2) - (7*prms.J2p*prms.R1^2*z*(3*prms.u - 3)*((prms.u + x)^2 + y^2 - 4*z^2))/(2*((prms.u + x)^2 + y^2 + z^2)^(9/2)) + (21*prms.J2s*prms.R2^2*prms.u*z*((prms.u + x - 1)^2 + y^2 - 4*z^2))/(2*((prms.u + x - 1)^2 + y^2 + z^2)^(9/2)));
% A(6,1) = -z*((3*(2*prms.u + 2*x)*(prms.u - 1))/(2*((prms.u + x)^2 + y^2 + z^2)^(5/2)) - (3*prms.u*(2*prms.u + 2*x - 2))/(2*((prms.u + x - 1)^2 + y^2 + z^2)^(5/2)) - (prms.J2p*prms.R1^2*(6*prms.u + 6*x)*(3*prms.u - 3))/(2*((prms.u + x)^2 + y^2 + z^2)^(7/2)) + (3*prms.J2s*prms.R2^2*prms.u*(6*prms.u + 6*x - 6))/(2*((prms.u + x - 1)^2 + y^2 + z^2)^(7/2)) + (7*prms.J2p*prms.R1^2*(2*prms.u + 2*x)*(3*prms.u - 3)*(3*(prms.u + x)^2 + 3*y^2 - 2*z^2))/(4*((prms.u + x)^2 + y^2 + z^2)^(9/2)) - (21*prms.J2s*prms.R2^2*prms.u*(2*prms.u + 2*x - 2)*(3*(prms.u + x - 1)^2 + 3*y^2 - 2*z^2))/(4*((prms.u + x - 1)^2 + y^2 + z^2)^(9/2)));
% A(6,2) = z*((3*prms.u*y)/((prms.u + x - 1)^2 + y^2 + z^2)^(5/2) - (3*y*(prms.u - 1))/((prms.u + x)^2 + y^2 + z^2)^(5/2) + (3*prms.J2p*prms.R1^2*y*(3*prms.u - 3))/((prms.u + x)^2 + y^2 + z^2)^(7/2) - (9*prms.J2s*prms.R2^2*prms.u*y)/((prms.u + x - 1)^2 + y^2 + z^2)^(7/2) + (21*prms.J2s*prms.R2^2*prms.u*y*(3*(prms.u + x - 1)^2 + 3*y^2 - 2*z^2))/(2*((prms.u + x - 1)^2 + y^2 + z^2)^(9/2)) - (7*prms.J2p*prms.R1^2*y*(3*prms.u - 3)*(3*(prms.u + x)^2 + 3*y^2 - 2*z^2))/(2*((prms.u + x)^2 + y^2 + z^2)^(9/2)));
% A(6,3) = z*((3*prms.u*z)/((prms.u + x - 1)^2 + y^2 + z^2)^(5/2) - (3*z*(prms.u - 1))/((prms.u + x)^2 + y^2 + z^2)^(5/2) - (2*prms.J2p*prms.R1^2*z*(3*prms.u - 3))/((prms.u + x)^2 + y^2 + z^2)^(7/2) + (6*prms.J2s*prms.R2^2*prms.u*z)/((prms.u + x - 1)^2 + y^2 + z^2)^(7/2) - (7*prms.J2p*prms.R1^2*z*(3*prms.u - 3)*(3*(prms.u + x)^2 + 3*y^2 - 2*z^2))/(2*((prms.u + x)^2 + y^2 + z^2)^(9/2)) + (21*prms.J2s*prms.R2^2*prms.u*z*(3*(prms.u + x - 1)^2 + 3*y^2 - 2*z^2))/(2*((prms.u + x - 1)^2 + y^2 + z^2)^(9/2))) + (prms.u - 1)/((prms.u + x)^2 + y^2 + z^2)^(3/2) - prms.u/((prms.u + x - 1)^2 + y^2 + z^2)^(3/2) + (prms.J2p*prms.R1^2*(3*prms.u - 3)*(3*(prms.u + x)^2 + 3*y^2 - 2*z^2))/(2*((prms.u + x)^2 + y^2 + z^2)^(7/2)) - (3*prms.J2s*prms.R2^2*prms.u*(3*(prms.u + x - 1)^2 + 3*y^2 - 2*z^2))/(2*((prms.u + x - 1)^2 + y^2 + z^2)^(7/2));

A = zeros(6,6);
A(1:3,4:6) = eye(3);
A(4,5) = 2*prms.n;
A(5,4) = -2*prms.n;
A(4,1) = (prms.u - 1)/((prms.u + x)^2 + y^2 + z^2)^(3/2) + prms.n^2 - prms.u/((prms.u + x - 1)^2 + y^2 + z^2)^(3/2) + (3*prms.u*(2*prms.u + 2*x - 2)^2)/(4*((prms.u + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*(2*prms.u + 2*x)^2*(prms.u - 1))/(4*((prms.u + x)^2 + y^2 + z^2)^(5/2)) - (prms.J2p*prms.R1^2*(prms.u - 1))/((prms.u + x)^2 + y^2 + z^2)^(5/2) + (prms.J2s*prms.R2^2*prms.u)/((prms.u + x - 1)^2 + y^2 + z^2)^(5/2) - (5*prms.J2s*prms.R2^2*prms.u*((prms.u + x - 1)^2 + y^2 - 2*z^2))/(2*((prms.u + x - 1)^2 + y^2 + z^2)^(7/2)) - (5*prms.J2s*prms.R2^2*prms.u*(2*prms.u + 2*x - 2)^2)/(2*((prms.u + x - 1)^2 + y^2 + z^2)^(7/2)) + (5*prms.J2p*prms.R1^2*(prms.u - 1)*((prms.u + x)^2 + y^2 - 2*z^2))/(2*((prms.u + x)^2 + y^2 + z^2)^(7/2)) + (5*prms.J2p*prms.R1^2*(2*prms.u + 2*x)^2*(prms.u - 1))/(2*((prms.u + x)^2 + y^2 + z^2)^(7/2)) + (35*prms.J2s*prms.R2^2*prms.u*(2*prms.u + 2*x - 2)^2*((prms.u + x - 1)^2 + y^2 - 2*z^2))/(8*((prms.u + x - 1)^2 + y^2 + z^2)^(9/2)) - (35*prms.J2p*prms.R1^2*(2*prms.u + 2*x)^2*(prms.u - 1)*((prms.u + x)^2 + y^2 - 2*z^2))/(8*((prms.u + x)^2 + y^2 + z^2)^(9/2));
A(4,2) = (3*prms.u*y*(2*prms.u + 2*x - 2))/(2*((prms.u + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*y*(2*prms.u + 2*x)*(prms.u - 1))/(2*((prms.u + x)^2 + y^2 + z^2)^(5/2)) - (5*prms.J2s*prms.R2^2*prms.u*y*(2*prms.u + 2*x - 2))/((prms.u + x - 1)^2 + y^2 + z^2)^(7/2) + (5*prms.J2p*prms.R1^2*y*(2*prms.u + 2*x)*(prms.u - 1))/((prms.u + x)^2 + y^2 + z^2)^(7/2) + (35*prms.J2s*prms.R2^2*prms.u*y*(2*prms.u + 2*x - 2)*((prms.u + x - 1)^2 + y^2 - 2*z^2))/(4*((prms.u + x - 1)^2 + y^2 + z^2)^(9/2)) - (35*prms.J2p*prms.R1^2*y*(2*prms.u + 2*x)*(prms.u - 1)*((prms.u + x)^2 + y^2 - 2*z^2))/(4*((prms.u + x)^2 + y^2 + z^2)^(9/2));
A(4,3) = (3*prms.u*z*(2*prms.u + 2*x - 2))/(2*((prms.u + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*z*(2*prms.u + 2*x)*(prms.u - 1))/(2*((prms.u + x)^2 + y^2 + z^2)^(5/2)) + (5*prms.J2s*prms.R2^2*prms.u*z*(2*prms.u + 2*x - 2))/(2*((prms.u + x - 1)^2 + y^2 + z^2)^(7/2)) - (5*prms.J2p*prms.R1^2*z*(2*prms.u + 2*x)*(prms.u - 1))/(2*((prms.u + x)^2 + y^2 + z^2)^(7/2)) + (35*prms.J2s*prms.R2^2*prms.u*z*(2*prms.u + 2*x - 2)*((prms.u + x - 1)^2 + y^2 - 2*z^2))/(4*((prms.u + x - 1)^2 + y^2 + z^2)^(9/2)) - (35*prms.J2p*prms.R1^2*z*(2*prms.u + 2*x)*(prms.u - 1)*((prms.u + x)^2 + y^2 - 2*z^2))/(4*((prms.u + x)^2 + y^2 + z^2)^(9/2));
A(5,1) = (3*prms.u*y*(2*prms.u + 2*x - 2))/(2*((prms.u + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*y*(2*prms.u + 2*x)*(prms.u - 1))/(2*((prms.u + x)^2 + y^2 + z^2)^(5/2)) - (5*prms.J2s*prms.R2^2*prms.u*y*(2*prms.u + 2*x - 2))/((prms.u + x - 1)^2 + y^2 + z^2)^(7/2) + (5*prms.J2p*prms.R1^2*y*(2*prms.u + 2*x)*(prms.u - 1))/((prms.u + x)^2 + y^2 + z^2)^(7/2) + (35*prms.J2s*prms.R2^2*prms.u*y*(2*prms.u + 2*x - 2)*((prms.u + x - 1)^2 + y^2 - 2*z^2))/(4*((prms.u + x - 1)^2 + y^2 + z^2)^(9/2)) - (35*prms.J2p*prms.R1^2*y*(2*prms.u + 2*x)*(prms.u - 1)*((prms.u + x)^2 + y^2 - 2*z^2))/(4*((prms.u + x)^2 + y^2 + z^2)^(9/2));
A(5,2) = (prms.u - 1)/((prms.u + x)^2 + y^2 + z^2)^(3/2) + prms.n^2 - prms.u/((prms.u + x - 1)^2 + y^2 + z^2)^(3/2) - (3*y^2*(prms.u - 1))/((prms.u + x)^2 + y^2 + z^2)^(5/2) + (3*prms.u*y^2)/((prms.u + x - 1)^2 + y^2 + z^2)^(5/2) - (prms.J2p*prms.R1^2*(prms.u - 1))/((prms.u + x)^2 + y^2 + z^2)^(5/2) + (prms.J2s*prms.R2^2*prms.u)/((prms.u + x - 1)^2 + y^2 + z^2)^(5/2) - (5*prms.J2s*prms.R2^2*prms.u*((prms.u + x - 1)^2 + y^2 - 2*z^2))/(2*((prms.u + x - 1)^2 + y^2 + z^2)^(7/2)) + (10*prms.J2p*prms.R1^2*y^2*(prms.u - 1))/((prms.u + x)^2 + y^2 + z^2)^(7/2) + (5*prms.J2p*prms.R1^2*(prms.u - 1)*((prms.u + x)^2 + y^2 - 2*z^2))/(2*((prms.u + x)^2 + y^2 + z^2)^(7/2)) - (10*prms.J2s*prms.R2^2*prms.u*y^2)/((prms.u + x - 1)^2 + y^2 + z^2)^(7/2) - (35*prms.J2p*prms.R1^2*y^2*(prms.u - 1)*((prms.u + x)^2 + y^2 - 2*z^2))/(2*((prms.u + x)^2 + y^2 + z^2)^(9/2)) + (35*prms.J2s*prms.R2^2*prms.u*y^2*((prms.u + x - 1)^2 + y^2 - 2*z^2))/(2*((prms.u + x - 1)^2 + y^2 + z^2)^(9/2));
A(5,3) = (3*prms.u*y*z)/((prms.u + x - 1)^2 + y^2 + z^2)^(5/2) - (3*y*z*(prms.u - 1))/((prms.u + x)^2 + y^2 + z^2)^(5/2) - (5*prms.J2p*prms.R1^2*y*z*(prms.u - 1))/((prms.u + x)^2 + y^2 + z^2)^(7/2) + (5*prms.J2s*prms.R2^2*prms.u*y*z)/((prms.u + x - 1)^2 + y^2 + z^2)^(7/2) - (35*prms.J2p*prms.R1^2*y*z*(prms.u - 1)*((prms.u + x)^2 + y^2 - 2*z^2))/(2*((prms.u + x)^2 + y^2 + z^2)^(9/2)) + (35*prms.J2s*prms.R2^2*prms.u*y*z*((prms.u + x - 1)^2 + y^2 - 2*z^2))/(2*((prms.u + x - 1)^2 + y^2 + z^2)^(9/2));
A(6,1) = (3*prms.u*z*(2*prms.u + 2*x - 2))/(2*((prms.u + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*z*(2*prms.u + 2*x)*(prms.u - 1))/(2*((prms.u + x)^2 + y^2 + z^2)^(5/2)) + (5*prms.J2s*prms.R2^2*prms.u*z*(2*prms.u + 2*x - 2))/(2*((prms.u + x - 1)^2 + y^2 + z^2)^(7/2)) - (5*prms.J2p*prms.R1^2*z*(2*prms.u + 2*x)*(prms.u - 1))/(2*((prms.u + x)^2 + y^2 + z^2)^(7/2)) + (35*prms.J2s*prms.R2^2*prms.u*z*(2*prms.u + 2*x - 2)*((prms.u + x - 1)^2 + y^2 - 2*z^2))/(4*((prms.u + x - 1)^2 + y^2 + z^2)^(9/2)) - (35*prms.J2p*prms.R1^2*z*(2*prms.u + 2*x)*(prms.u - 1)*((prms.u + x)^2 + y^2 - 2*z^2))/(4*((prms.u + x)^2 + y^2 + z^2)^(9/2));
A(6,2) = (3*prms.u*y*z)/((prms.u + x - 1)^2 + y^2 + z^2)^(5/2) - (3*y*z*(prms.u - 1))/((prms.u + x)^2 + y^2 + z^2)^(5/2) - (5*prms.J2p*prms.R1^2*y*z*(prms.u - 1))/((prms.u + x)^2 + y^2 + z^2)^(7/2) + (5*prms.J2s*prms.R2^2*prms.u*y*z)/((prms.u + x - 1)^2 + y^2 + z^2)^(7/2) - (35*prms.J2p*prms.R1^2*y*z*(prms.u - 1)*((prms.u + x)^2 + y^2 - 2*z^2))/(2*((prms.u + x)^2 + y^2 + z^2)^(9/2)) + (35*prms.J2s*prms.R2^2*prms.u*y*z*((prms.u + x - 1)^2 + y^2 - 2*z^2))/(2*((prms.u + x - 1)^2 + y^2 + z^2)^(9/2));
A(6,3) = (prms.u - 1)/((prms.u + x)^2 + y^2 + z^2)^(3/2) - prms.u/((prms.u + x - 1)^2 + y^2 + z^2)^(3/2) - (3*z^2*(prms.u - 1))/((prms.u + x)^2 + y^2 + z^2)^(5/2) + (3*prms.u*z^2)/((prms.u + x - 1)^2 + y^2 + z^2)^(5/2) + (2*prms.J2p*prms.R1^2*(prms.u - 1))/((prms.u + x)^2 + y^2 + z^2)^(5/2) - (2*prms.J2s*prms.R2^2*prms.u)/((prms.u + x - 1)^2 + y^2 + z^2)^(5/2) - (5*prms.J2s*prms.R2^2*prms.u*((prms.u + x - 1)^2 + y^2 - 2*z^2))/(2*((prms.u + x - 1)^2 + y^2 + z^2)^(7/2)) - (20*prms.J2p*prms.R1^2*z^2*(prms.u - 1))/((prms.u + x)^2 + y^2 + z^2)^(7/2) + (5*prms.J2p*prms.R1^2*(prms.u - 1)*((prms.u + x)^2 + y^2 - 2*z^2))/(2*((prms.u + x)^2 + y^2 + z^2)^(7/2)) + (20*prms.J2s*prms.R2^2*prms.u*z^2)/((prms.u + x - 1)^2 + y^2 + z^2)^(7/2) - (35*prms.J2p*prms.R1^2*z^2*(prms.u - 1)*((prms.u + x)^2 + y^2 - 2*z^2))/(2*((prms.u + x)^2 + y^2 + z^2)^(9/2)) + (35*prms.J2s*prms.R2^2*prms.u*z^2*((prms.u + x - 1)^2 + y^2 - 2*z^2))/(2*((prms.u + x - 1)^2 + y^2 + z^2)^(9/2));


%%% Calculate stmDot and output result
stmDot = A*stm;
dX(7:end) = reshape(stmDot,36,1);

end

