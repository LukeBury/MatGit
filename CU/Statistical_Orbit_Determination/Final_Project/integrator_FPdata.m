function [dX] = integrator_FPdata(t,X,u,rB1,rB2,rad2,tNorm,rNorm)
dX = zeros(6,1);

%%% Unpack the Hopper state vector
rP = X(1:3); % Hopper Position
vP = X(4:6); % Hopper Velocity

%%% Assigning Variables
x = rP(1); y = rP(2); z = rP(3);
dx = vP(1); dy = vP(2); dz = vP(3);

%%% Creating hopper distances to bodies
r1 = sqrt((rP(1)-rB1(1))^2 + (rP(2)-rB1(2))^2 + (rP(3)-rB1(3))^2);
r2 = sqrt((rP(1)-rB2(1))^2 + (rP(2)-rB2(2))^2 + (rP(3)-rB2(3))^2);

%%% Equations of Motion (CR3BP)
ddx_CR3BP = 2*dy + x - (1-u)*(x+u)/(r1^3) - u*(x+u-1)/(r2^3);
ddy_CR3BP = -2*dx + y -((1-u)/(r1^3) + u/(r2^3))*y;
ddz_CR3BP = -((1-u)/(r1^3) + u/(r2^3))*z;

%%% Equations of Motion (J2)
% ddx_J2 = (3*J2*RE^2*u*x)/(2*(x^2 + y^2 + z^2)^(5/2)) - (15*J2*RE^2*u*x*z^2)/(2*(x^2 + y^2 + z^2)^(7/2))
% ddy_J2 =(3*J2*RE^2*u*y)/(2*(x^2 + y^2 + z^2)^(5/2)) - (15*J2*RE^2*u*y*z^2)/(2*(x^2 + y^2 + z^2)^(7/2))
% ddz_J2 =(9*J2*RE^2*u*z)/(2*(x^2 + y^2 + z^2)^(5/2)) - (15*J2*RE^2*u*z^3)/(2*(x^2 + y^2 + z^2)^(7/2))

%%% Storing EOMs
ddxP = [ddx_CR3BP; ddy_CR3BP; ddz_CR3BP];

%%% Output the derivative of the state
dX(1:3) = vP; % km/s
dX(4:6) = ddxP; % km/s^2
end

