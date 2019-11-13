% Integrates concentric 4BP, where m3 is in orbit m1 (not barycenter)
% T2 and T3 are periods of m2 and m3 in same units
% r13 is distance from m1 to m3 normatlized by distance from m1 to m2
% mu3 is m3/m1+m2, theta0 is initial phase angle (m3 leading m2 for
% positive angle)
% Assumes m2,m3 << m1. m1 and m2 orbit in circles about their barycenter,
% and M3 ORBITS M1, NOT M1,M2 BARYCENTER. According to previous assumption,
% this affect is small. TO INSTEAD ORBIT BARYCENTER, CHANGE DEFINITION OF
% "a" TO BE a = r13*cos(theta).

function xdot_CR4BP = xdot_CR4BP(t,x, mu, mu3, r13, theta0, T2, T3)

theta = theta0 + (T2/T3 - 1)*t;
a = r13*cos(theta) - mu;
b = r13*sin(theta);
r1 = sqrt((x(1)+mu)^2+x(2)^2+x(3)^2);
r2 = sqrt((x(1)-1+mu)^2+x(2)^2+x(3)^2);
r3 = sqrt((x(1)-a)^2 + (x(2)-b)^2 + x(3)^2);
Omega_x = x(1)-(1-mu)*(x(1)+mu)/r1^3-mu*(x(1)-1+mu)/r2^3 - mu3*(x(1)-a)/r3^3;
Omega_y = x(2)-(1-mu)*x(2)/r1^3-mu*x(2)/r2^3 - mu3*(x(2)-b)/r3^3;
Omega_z = -x(3)*(1-mu)/r1^3 - x(3)*mu/r2^3 - mu3*x(3)/r3^3;

xdot_CR4BP = zeros(6,1);
xdot_CR4BP(1) = x(4);
xdot_CR4BP(2) = x(5);
xdot_CR4BP(3) = x(6);
xdot_CR4BP(4) = 2*x(5) + Omega_x;
xdot_CR4BP(5) = -2*x(4) + Omega_y;
xdot_CR4BP(6) = Omega_z;
end