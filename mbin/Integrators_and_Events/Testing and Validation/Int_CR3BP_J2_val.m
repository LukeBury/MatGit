function [dX] = Int_CR3BP_J2_val(t,X,u,R1,J21)

% ************************************************************************
% ************************************************************************
% IMPORTANT NOTE: this integrator is simply a normalized version of the
% 2-body +J2 rotating problem. It was designed for the express purpose of
% validating my implementation of J2 in the rotating frame and shouldn't be
% used for further application.
% ************************************************************************
% ************************************************************************



%%% For numerical integration in the normalized CR3BP
%%% Inputs:
%          t - normalized time vector
%          X - initial state [6x1]
%          u - mass ratio of CR3BP system
%          R1 - radius of primary body
%          J21 - J2 of primary body

%%% Defining gravitational strengths
u1 = 1-u;
u2 = 0;

%%% Preallocate state output
dX = zeros(6,1);

%%% Unpack the barycentric state vector
x = X(1); y = X(2); z = X(3); % Position
dx = X(4); dy = X(5); dz = X(6); % Velocity

%%% Define position of bodies
x1 = 0;
x2 = 0;

%%% Distances to primary (1) and secondary (2) bodies
r1 = sqrt((x)^2 + y^2 + z^2);
r2 = sqrt((x)^2 + y^2 + z^2);

%%% Define Gamma values for shortening J2 terms
gamma1 = 3*(u1)*R1^2*J21*(5*z^2-r1^2)/(2*r1^7);
gamma2 = 0;
gamma1_z = 3*(u1)*R1^2*J21*(5*z^2-3*r1^2)/(2*r1^7);
gamma2_z = 0;

%%% Equations of Motion
ddx = 2*dy + x + ((u1)/r1^3 - gamma1)*(-x);
ddy = -2*dx + y*(-(u1)/r1^3 + gamma1 + 1);
ddz = z*(-(u1)/r1^3 + gamma1_z);

%%% Output the derivative of the state
dX(1:3) = [dx; dy; dz]; % km/s
dX(4:6) = [ddx; ddy; ddz]; % km/s^2

end
