function [dX] = CR3BP_norm(~, X)
% Don Kuettel

% This function is used to complete the numerical integration of the 
% normalized CR3BP.

v = X(7);

dX = zeros(7 + 6^2,1);

x = X(1);           % x-position
y = X(2);           % y-position
z = X(3);           % z-position
vx = X(4);          % x-velocity
vy = X(5);          % y-velocity
vz = X(6);          % z-velocity

r1 = sqrt((x + v)^2 + y*y + z*z);
r2 = sqrt((x - 1 + v)^2 + y*y + z*z);

% Derivatives
dX(1) = vx;
dX(2) = vy;
dX(3) = vz;
dX(4) = 2*vy + x - (1 - v)/r1^3*(x + v) - v/r2^3*(x - 1 + v);
dX(5) = -2*vx + y - (1 - v)/r1^3*y - v/r2^3*y;
dX(6) = -(1 - v)/r1^3*z - v/r2^3*z;

% Keeping v constant
dX(7) = 0;

% Here is where we reshape the STM from a 1x6^2 to an 6x6 matrix
phi = reshape(X(8:end),6,6);

%  Evaluate the partial matrix A(t)
A = zeros(6,6);

A(1,4) = 1; A(2,5) = 1; A(3,6) = 1; A(4,5) = 2; A(5,4) = -2;

A(4,1) = (v - 1)/((v + x)^2 + y^2 + z^2)^(3/2) - v/((v + x - 1)^2 + y^2 + z^2)^(3/2) + (3*v*(2*v + 2*x - 2)^2)/(4*((v + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*(2*v + 2*x)^2*(v - 1))/(4*((v + x)^2 + y^2 + z^2)^(5/2)) + 1;

A(4,2) = (3*v*y*(2*v + 2*x - 2))/(2*((v + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*y*(2*v + 2*x)*(v - 1))/(2*((v + x)^2 + y^2 + z^2)^(5/2));

A(4,3) = (3*v*z*(2*v + 2*x - 2))/(2*((v + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*z*(2*v + 2*x)*(v - 1))/(2*((v + x)^2 + y^2 + z^2)^(5/2));

A(5,1) = (3*v*y*(2*v + 2*x - 2))/(2*((v + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*y*(2*v + 2*x)*(v - 1))/(2*((v + x)^2 + y^2 + z^2)^(5/2));

A(5,2) = (v - 1)/((v + x)^2 + y^2 + z^2)^(3/2) - v/((v + x - 1)^2 + y^2 + z^2)^(3/2) - (3*y^2*(v - 1))/((v + x)^2 + y^2 + z^2)^(5/2) + (3*v*y^2)/((v + x - 1)^2 + y^2 + z^2)^(5/2) + 1;

A(5,3) = (3*v*y*z)/((v + x - 1)^2 + y^2 + z^2)^(5/2) - (3*y*z*(v - 1))/((v + x)^2 + y^2 + z^2)^(5/2);

A(6,1) = (3*v*z*(2*v + 2*x - 2))/(2*((v + x - 1)^2 + y^2 + z^2)^(5/2)) - (3*z*(2*v + 2*x)*(v - 1))/(2*((v + x)^2 + y^2 + z^2)^(5/2));

A(6,2) = (3*v*y*z)/((v + x - 1)^2 + y^2 + z^2)^(5/2) - (3*y*z*(v - 1))/((v + x)^2 + y^2 + z^2)^(5/2);

A(6,3) = (v - 1)/((v + x)^2 + y^2 + z^2)^(3/2) - v/((v + x - 1)^2 + y^2 + z^2)^(3/2) - (3*z^2*(v - 1))/((v + x)^2 + y^2 + z^2)^(5/2) + (3*v*z^2)/((v + x - 1)^2 + y^2 + z^2)^(5/2);

% Solve for the time derivative of the STM using A(t)
phiDot = A*phi;

% Reshape the output and write it to the output array
dX(8:end) = reshape(phiDot,6*6,1);

end
