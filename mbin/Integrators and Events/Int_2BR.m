function [dX] = Int_2BR(t,X,u,dTheta)
%%% For numerical integration in the normalized CR3BP
%%% Inputs:
%          t - time vector
%          X - initial state [6x1]
%          u - gravitational parameter of primary
%          dTheta - ang. vel. of rotating frame wrt inertial (rad/s)

%%% Preallocate state output
dX = zeros(6,1);

%%% Unpack the barycentric state vector
x = X(1); y = X(2); z = X(3); % Position
dx = X(4); dy = X(5); dz = X(6); % Velocity

%%% Distances to primary (1) and secondary (2) bodies
r = sqrt(x^2 + y^2 + z^2);

%%% Equations of Motion
ddx = 2*dy*dTheta + x*dTheta^2 - u*x/(r^3);
ddy = -2*dx*dTheta + y*dTheta^2 - u*y/(r^3);
ddz = -u*z/(r^3);

%%% Output the derivative of the state
dX(1:3) = [dx; dy; dz]; % km/s
dX(4:6) = [ddx; ddy; ddz]; % km/s^2
end

