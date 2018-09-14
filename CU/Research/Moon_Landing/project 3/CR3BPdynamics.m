function [ f ] = CR3BPdynamics(X, v)

f = zeros(6,1);

x = X(1);           % x-position
y = X(2);           % y-position
z = X(3);           % z-position
vx = X(4);          % x-velocity
vy = X(5);          % y-velocity
vz = X(6);          % z-velocity

r1 = sqrt((x + v)^2 + y*y + z*z);
r2 = sqrt((x - 1 + v)^2 + y*y + z*z);

% Derivatives
f(1) = vx;
f(2) = vy;
f(3) = vz;
f(4) = 2*vy + x - (1 - v)/r1^3*(x + v) - v/r2^3*(x - 1 + v);
f(5) = -2*vx + y - (1 - v)/r1^3*y - v/r2^3*y;
f(6) = -(1 - v)/r1^3*z - v/r2^3*z;

end

