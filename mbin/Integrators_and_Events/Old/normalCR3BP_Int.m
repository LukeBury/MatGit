function [dX] = normalCR3BP_Int(t,X,u,R2)
%%% For numerical integration in the normalized CR3BP
%%% Inputs:
%          t - normalized time vector
%          X - initial state [6x1]
%          u - mass ratio of CR3BP system
%          R2 - radius of secondary body

%%% Preallocate state output
dX = zeros(6,1);

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
end

