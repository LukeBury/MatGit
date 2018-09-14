function [dX] = int_CR3BnSTM(t,X,u)
%%% For numerical integration in the normalized CR3BP
%%% Inputs:
%          t - normalized time vector
%          X - initial state [42x1]
%          u - mass ratio of CR3BP system

%%% Preallocate state output
dX = zeros(6+36,1);

%%% Unpack the barycentric state vector
x = X(1); y = X(2); z = X(3); % Position
dx = X(4); dy = X(5); dz = X(6); % Velocity

%%% Distances to primary (1) and secondary (2) bodies
r1 = sqrt((x+u)^2 + y^2 + z^2);
r2 = sqrt((x+u-1)^2 + y^2 + z^2);

%%% Equations of Motion
ddx = 2*dy + x - (1-u)*(x+u)/(r1^3) - u*(x+u-1)/(r2^3);
ddy = -2*dx + y -((1-u)/(r1^3) + u/(r2^3))*y;
ddz = -((1-u)/(r1^3) + u/(r2^3))*z;

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
A(4,1) = (u - 1)/((u + x)^2 + y^2 + z^2)^(3/2) - u/((u + x - 1)^2 + y^2 + z^2)^(3/2) + (3*u*(2*u + 2*x - 2)*(u + x - 1))/(2*((u + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*(2*u + 2*x)*(u + x)*(u - 1))/(2*((u + x)^2 + y^2 + z^2)^(5/2)) + 1;
A(4,2) = (3*u*y*(u + x - 1))/((u + x - 1)^2 + y^2 + z^2)^(5/2) - (3*y*(u + x)*(u - 1))/((u + x)^2 + y^2 + z^2)^(5/2);
A(4,3) = (3*u*z*(u + x - 1))/((u + x - 1)^2 + y^2 + z^2)^(5/2) - (3*z*(u + x)*(u - 1))/((u + x)^2 + y^2 + z^2)^(5/2);
A(5,1) = -y*((3*(2*u + 2*x)*(u - 1))/(2*((u + x)^2 + y^2 + z^2)^(5/2)) - (3*u*(2*u + 2*x - 2))/(2*((u + x - 1)^2 + y^2 + z^2)^(5/2)));
A(5,2) = (u - 1)/((u + x)^2 + y^2 + z^2)^(3/2) + y*((3*u*y)/((u + x - 1)^2 + y^2 + z^2)^(5/2) - (3*y*(u - 1))/((u + x)^2 + y^2 + z^2)^(5/2)) - u/((u + x - 1)^2 + y^2 + z^2)^(3/2) + 1;
A(5,3) = y*((3*u*z)/((u + x - 1)^2 + y^2 + z^2)^(5/2) - (3*z*(u - 1))/((u + x)^2 + y^2 + z^2)^(5/2));
A(6,1) = -z*((3*(2*u + 2*x)*(u - 1))/(2*((u + x)^2 + y^2 + z^2)^(5/2)) - (3*u*(2*u + 2*x - 2))/(2*((u + x - 1)^2 + y^2 + z^2)^(5/2)));
A(6,2) = z*((3*u*y)/((u + x - 1)^2 + y^2 + z^2)^(5/2) - (3*y*(u - 1))/((u + x)^2 + y^2 + z^2)^(5/2));
A(6,3) = (u - 1)/((u + x)^2 + y^2 + z^2)^(3/2) + z*((3*u*z)/((u + x - 1)^2 + y^2 + z^2)^(5/2) - (3*z*(u - 1))/((u + x)^2 + y^2 + z^2)^(5/2)) - u/((u + x - 1)^2 + y^2 + z^2)^(3/2);


%%% Calculate stmDot and output result
stmDot = A*stm;
dX(7:end) = reshape(stmDot,36,1);

end

