function [dX] = Int_2BI_finite(t,X,u,tMag)
%%% For numerical integration of ECI state under influence of standard
%%% 2-body gravity and constant in-track acceleration
%%% Inputs:
%          t - time vector
%          X - initial state [6x1]
%          u - gravitational parameter of primary body
%          tMag - in-track thrust value

%%% Preallocate state output
dX = zeros(6,1);

%%% Unpack the state vector
x = X(1); y = X(2); z = X(3); % Position
dx = X(4); dy = X(5); dz = X(6); % Velocity

%%% Distances from body to spacecraft
r = sqrt(x^2 + y^2 + z^2);

%%% finding direction of motion
vHat = [dx, dy, dz] ./ norm([dx, dy, dz]);

%%% Finding acceleration from thrust
aThrust = vHat .* tMag;

%%% Equations of Motion
ddx = x*(-u/r^3) + aThrust(1);
ddy = y*(-u/r^3) + aThrust(2);
ddz = z*(-u/r^3) + aThrust(3);

%%% Output the derivative of the state
dX(1:3) = [dx; dy; dz]; % km/s
dX(4:6) = [ddx; ddy; ddz]; % km/s^2

end
