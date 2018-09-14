function [dX] = int_CR3BnSTM_J2(t,X,u,R1,R2,J21,J22)
%%% For numerical integration in the normalized CR3BP
%%% Inputs:
%          t - normalized time vector
%          X - initial state [42x1]
%          u - mass ratio of CR3BP system
%          R1 - radius of primary body
%          R2 - radius of secondary body
%          J21 - J2 of primary body
%          J22 - J2 of secondary body
% =======================================================================
%%% Preallocate state output
dX = zeros(6+36,1);

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

%%% Reshaping STM
stm = reshape(X(7:end),6,6);

%%% Evaluate A matrix
A = zeros(6,6);
A(1:3,4:6) = eye(3);
A(4,5) = 2;
A(5,4) = -2;
A(4,1) = (u + x - 1)*((3*u*(2*u + 2*x - 2))/(2*((u + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*J22*R2^2*u*(2*u + 2*x - 2))/(2*((u + x - 1)^2 + y^2 + z^2)^(7/2)) + (21*J22*R2^2*u*(2*u + 2*x - 2)*((u + x - 1)^2 + y^2 - 4*z^2))/(4*((u + x - 1)^2 + y^2 + z^2)^(9/2))) + (u - 1)/((u + x)^2 + y^2 + z^2)^(3/2) - (u + x)*((3*(2*u + 2*x)*(u - 1))/(2*((u + x)^2 + y^2 + z^2)^(5/2)) - (J21*R1^2*(2*u + 2*x)*(3*u - 3))/(2*((u + x)^2 + y^2 + z^2)^(7/2)) + (7*J21*R1^2*(2*u + 2*x)*(3*u - 3)*((u + x)^2 + y^2 - 4*z^2))/(4*((u + x)^2 + y^2 + z^2)^(9/2))) - u/((u + x - 1)^2 + y^2 + z^2)^(3/2) - (3*J22*R2^2*u*((u + x - 1)^2 + y^2 - 4*z^2))/(2*((u + x - 1)^2 + y^2 + z^2)^(7/2)) + (J21*R1^2*(3*u - 3)*((u + x)^2 + y^2 - 4*z^2))/(2*((u + x)^2 + y^2 + z^2)^(7/2)) + 1;
A(4,2) = (u + x - 1)*((3*u*y)/((u + x - 1)^2 + y^2 + z^2)^(5/2) - (3*J22*R2^2*u*y)/((u + x - 1)^2 + y^2 + z^2)^(7/2) + (21*J22*R2^2*u*y*((u + x - 1)^2 + y^2 - 4*z^2))/(2*((u + x - 1)^2 + y^2 + z^2)^(9/2))) - (u + x)*((3*y*(u - 1))/((u + x)^2 + y^2 + z^2)^(5/2) - (J21*R1^2*y*(3*u - 3))/((u + x)^2 + y^2 + z^2)^(7/2) + (7*J21*R1^2*y*(3*u - 3)*((u + x)^2 + y^2 - 4*z^2))/(2*((u + x)^2 + y^2 + z^2)^(9/2)));
A(4,3) = (u + x - 1)*((3*u*z)/((u + x - 1)^2 + y^2 + z^2)^(5/2) + (12*J22*R2^2*u*z)/((u + x - 1)^2 + y^2 + z^2)^(7/2) + (21*J22*R2^2*u*z*((u + x - 1)^2 + y^2 - 4*z^2))/(2*((u + x - 1)^2 + y^2 + z^2)^(9/2))) - (u + x)*((3*z*(u - 1))/((u + x)^2 + y^2 + z^2)^(5/2) + (4*J21*R1^2*z*(3*u - 3))/((u + x)^2 + y^2 + z^2)^(7/2) + (7*J21*R1^2*z*(3*u - 3)*((u + x)^2 + y^2 - 4*z^2))/(2*((u + x)^2 + y^2 + z^2)^(9/2)));
A(5,1) = -y*((3*(2*u + 2*x)*(u - 1))/(2*((u + x)^2 + y^2 + z^2)^(5/2)) - (3*u*(2*u + 2*x - 2))/(2*((u + x - 1)^2 + y^2 + z^2)^(5/2)) - (J21*R1^2*(2*u + 2*x)*(3*u - 3))/(2*((u + x)^2 + y^2 + z^2)^(7/2)) + (3*J22*R2^2*u*(2*u + 2*x - 2))/(2*((u + x - 1)^2 + y^2 + z^2)^(7/2)) - (21*J22*R2^2*u*(2*u + 2*x - 2)*((u + x - 1)^2 + y^2 - 4*z^2))/(4*((u + x - 1)^2 + y^2 + z^2)^(9/2)) + (7*J21*R1^2*(2*u + 2*x)*(3*u - 3)*((u + x)^2 + y^2 - 4*z^2))/(4*((u + x)^2 + y^2 + z^2)^(9/2)));
A(5,2) = (u - 1)/((u + x)^2 + y^2 + z^2)^(3/2) + y*((3*u*y)/((u + x - 1)^2 + y^2 + z^2)^(5/2) - (3*y*(u - 1))/((u + x)^2 + y^2 + z^2)^(5/2) + (J21*R1^2*y*(3*u - 3))/((u + x)^2 + y^2 + z^2)^(7/2) - (3*J22*R2^2*u*y)/((u + x - 1)^2 + y^2 + z^2)^(7/2) - (7*J21*R1^2*y*(3*u - 3)*((u + x)^2 + y^2 - 4*z^2))/(2*((u + x)^2 + y^2 + z^2)^(9/2)) + (21*J22*R2^2*u*y*((u + x - 1)^2 + y^2 - 4*z^2))/(2*((u + x - 1)^2 + y^2 + z^2)^(9/2))) - u/((u + x - 1)^2 + y^2 + z^2)^(3/2) - (3*J22*R2^2*u*((u + x - 1)^2 + y^2 - 4*z^2))/(2*((u + x - 1)^2 + y^2 + z^2)^(7/2)) + (J21*R1^2*(3*u - 3)*((u + x)^2 + y^2 - 4*z^2))/(2*((u + x)^2 + y^2 + z^2)^(7/2)) + 1;
A(5,3) = y*((3*u*z)/((u + x - 1)^2 + y^2 + z^2)^(5/2) - (3*z*(u - 1))/((u + x)^2 + y^2 + z^2)^(5/2) - (4*J21*R1^2*z*(3*u - 3))/((u + x)^2 + y^2 + z^2)^(7/2) + (12*J22*R2^2*u*z)/((u + x - 1)^2 + y^2 + z^2)^(7/2) - (7*J21*R1^2*z*(3*u - 3)*((u + x)^2 + y^2 - 4*z^2))/(2*((u + x)^2 + y^2 + z^2)^(9/2)) + (21*J22*R2^2*u*z*((u + x - 1)^2 + y^2 - 4*z^2))/(2*((u + x - 1)^2 + y^2 + z^2)^(9/2)));
A(6,1) = -z*((3*(2*u + 2*x)*(u - 1))/(2*((u + x)^2 + y^2 + z^2)^(5/2)) - (3*u*(2*u + 2*x - 2))/(2*((u + x - 1)^2 + y^2 + z^2)^(5/2)) - (J21*R1^2*(6*u + 6*x)*(3*u - 3))/(2*((u + x)^2 + y^2 + z^2)^(7/2)) + (3*J22*R2^2*u*(6*u + 6*x - 6))/(2*((u + x - 1)^2 + y^2 + z^2)^(7/2)) + (7*J21*R1^2*(2*u + 2*x)*(3*u - 3)*(3*(u + x)^2 + 3*y^2 - 2*z^2))/(4*((u + x)^2 + y^2 + z^2)^(9/2)) - (21*J22*R2^2*u*(2*u + 2*x - 2)*(3*(u + x - 1)^2 + 3*y^2 - 2*z^2))/(4*((u + x - 1)^2 + y^2 + z^2)^(9/2)));
A(6,2) = z*((3*u*y)/((u + x - 1)^2 + y^2 + z^2)^(5/2) - (3*y*(u - 1))/((u + x)^2 + y^2 + z^2)^(5/2) + (3*J21*R1^2*y*(3*u - 3))/((u + x)^2 + y^2 + z^2)^(7/2) - (9*J22*R2^2*u*y)/((u + x - 1)^2 + y^2 + z^2)^(7/2) + (21*J22*R2^2*u*y*(3*(u + x - 1)^2 + 3*y^2 - 2*z^2))/(2*((u + x - 1)^2 + y^2 + z^2)^(9/2)) - (7*J21*R1^2*y*(3*u - 3)*(3*(u + x)^2 + 3*y^2 - 2*z^2))/(2*((u + x)^2 + y^2 + z^2)^(9/2)));
A(6,3) = z*((3*u*z)/((u + x - 1)^2 + y^2 + z^2)^(5/2) - (3*z*(u - 1))/((u + x)^2 + y^2 + z^2)^(5/2) - (2*J21*R1^2*z*(3*u - 3))/((u + x)^2 + y^2 + z^2)^(7/2) + (6*J22*R2^2*u*z)/((u + x - 1)^2 + y^2 + z^2)^(7/2) - (7*J21*R1^2*z*(3*u - 3)*(3*(u + x)^2 + 3*y^2 - 2*z^2))/(2*((u + x)^2 + y^2 + z^2)^(9/2)) + (21*J22*R2^2*u*z*(3*(u + x - 1)^2 + 3*y^2 - 2*z^2))/(2*((u + x - 1)^2 + y^2 + z^2)^(9/2))) + (u - 1)/((u + x)^2 + y^2 + z^2)^(3/2) - u/((u + x - 1)^2 + y^2 + z^2)^(3/2) + (J21*R1^2*(3*u - 3)*(3*(u + x)^2 + 3*y^2 - 2*z^2))/(2*((u + x)^2 + y^2 + z^2)^(7/2)) - (3*J22*R2^2*u*(3*(u + x - 1)^2 + 3*y^2 - 2*z^2))/(2*((u + x - 1)^2 + y^2 + z^2)^(7/2));


%%% Calculate stmDot and output result
stmDot = A*stm;
dX(7:end) = reshape(stmDot,36,1);

end

