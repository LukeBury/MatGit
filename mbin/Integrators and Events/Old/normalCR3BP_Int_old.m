function [dX] = normalCR3BP_Int(t,X,u,rB1,rB2,rad2)
%%% Inputs:
%          t - normalized time vector
%          X - initial state [6x1]
%          u - mass ratio of CR3BP system
%          rB1 - radius of primary body
%          rB2 - radius of secondary body
%          t - normalized time vector

dX = zeros(6,1);

%%% Unpack the barycentric state vector
x = X(1); y = X(2); z = X(3); % Position
dx = X(4); dy = X(5); dz = X(6); % Velocity

%%% Distances to primary (1) and secondary (2) bodies
r1 = sqrt((x-rB1(1))^2 + (y-rB1(2))^2 + (z-rB1(3))^2);
r2 = sqrt((x-rB2(1))^2 + (y-rB2(2))^2 + (z-rB2(3))^2);

%%% Equations of Motion
ddx = 2*dy + x - (1-u)*(x+u)/(r1^3) - u*(x+u-1)/(r2^3);
ddy = -2*dx + y -((1-u)/(r1^3) + u/(r2^3))*y;
ddz = -((1-u)/(r1^3) + u/(r2^3))*z;

%%% Storing EOMs
ddxP = [ddx; ddy; ddz];

%%% Output the derivative of the state
dX(1:3) = [dx; dy; dz]; % km/s
dX(4:6) = ddxP; % km/s^2
end

