function [dX] = Int_CR3BP_val(t,X,u)

% ************************************************************************
% ************************************************************************
% IMPORTANT NOTE: this integrator is simply a normalized version of the
% 2-body rotating problem. It was designed for the express purpose of
% validating my implementation of J2 in the rotating frame and shouldn't be
% used for further application.
% ************************************************************************
% ************************************************************************

%%% For numerical integration in the normalized CR3BP
%%% Inputs:
%          t - normalized time vector
%          X - initial state [6x1]
%          u - mass ratio of CR3BP system

%%% Defining gravitational strengths
u1 = 1-u;

%%% Preallocate state output
dX = zeros(6,1);

%%% Unpack the barycentric state vector
x = X(1); y = X(2); z = X(3); % Position
dx = X(4); dy = X(5); dz = X(6); % Velocity

%%% Distances to primary (1) and secondary (2) bodies
r1 = sqrt(x^2 + y^2 + z^2);

%%% Equations of Motion
ddx = 2*dy + x - u1*x/(r1^3);
ddy = -2*dx + y - u1*y/(r1^3);
ddz = -u1*z/(r1^3);

%%% Output the derivative of the state
dX(1:3) = [dx; dy; dz]; % km/s
dX(4:6) = [ddx; ddy; ddz]; % km/s^2
end

