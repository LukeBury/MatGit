function [dX] = Int_CR3Bn_perturb(t,X,prms)
%%% For numerical integration in the normalized CR3BP with J2 of each body
%%% Inputs:
%          t - normalized time vector
%          X - initial state [6x1]
%          prms - (u, R1_n, R2_n, J21, J22, J31, J32, J41, J42)
%
%
%          u - mass ratio of CR3BP system
%          R1 - radius of primary body
%          R2 - radius of secondary body
%          J21 - J2 of primary body
%          J22 - J2 of secondary body
%          J31 - J3 of primary body
%          J32 - J3 of secondary body
%          J41 - J4 of primary body
%          J42 - J4 of secondary body
% =======================================================================
%%% Preallocate state output
dX = zeros(6,1);

%%% Unpacking parameters (for brevity/clarity)
u   = prms.u;
R1  = prms.R1_n;
R2  = prms.R2_n;
J21 = prms.J21;
J22 = prms.J22;
J31 = prms.J31;
J32 = prms.J32;
J41 = prms.J41;
J42 = prms.J42;

%%% Distances to primary (1) and secondary (2) bodies
r1 = sqrt((X(1)+u)^2 + X(2)^2 + X(3)^2);
r2 = sqrt((X(1)+u-1)^2 + X(2)^2 + X(3)^2);

%%% Calculating effects of J2, J3, J4
gamma21 = 3*(1-u)*R1^2*J21*(5*X(3)^2-r1^2)/(2*r1^7);
gamma21_z = 3*(1-u)*R1^2*J21*(5*X(3)^2-3*r1^2)/(2*r1^7);

gamma22 = 3*u*R2^2*J22*(5*X(3)^2-r2^2)/(2*r2^7);
gamma22_z = 3*u*R2^2*J22*(5*X(3)^2-3*r2^2)/(2*r2^7);

gamma31 = R1^3*J31*(1-u)*5*X(3)*(7*X(3)^2-3*r1^2)/(2*r1^9);
gamma31_z = R1^3*J31*(1-u)*(3*r1^4-30*r1^2*X(3)^2+35*X(3)^4)/(2*r1^9);

gamma32 = R2^3*J32*u*5*X(3)*(7*X(3)^2-3*r2^2)/(2*r2^9);
gamma32_z = R2^3*J32*u*(3*r2^4-30*r2^2*X(3)^2+35*X(3)^4)/(2*r2^9);

gamma41 = R1^4*J41*(1-u)*15*(r1^4-14*r1^2*X(3)^2+21*X(3)^4)/(8*r1^11);
gamma41_z = R1^4*J41*(1-u)*5*(15*r1^4-70*r1^2*X(3)^2+63*X(3)^4)/(8*r1^11);

gamma42 = R2^4*J42*u*15*(r2^4-14*r2^2*X(3)^2+21*X(3)^4)/(8*r2^11);
gamma42_z = R2^4*J42*u*5*(15*r2^4-70*r2^2*X(3)^2+63*X(3)^4)/(8*r2^11);

%%% Equations of Motion
ddx = 2*X(5) + X(1) + ((1-u)/r1^3-gamma21 - gamma31 - gamma41)*((-u)-X(1)) + (u/r2^3 - gamma22 - gamma32 - gamma42)*((1-u)-X(1));
ddy = -2*X(4) + X(2) - X(2)*((1-u)/r1^3 + u/r2^3 - gamma21 - gamma22 - gamma31 - gamma32 - gamma41 - gamma42);
ddz = X(3)*(-(1-u)/r1^3 - u/r2^3 + gamma21_z + gamma22_z + gamma41_z + gamma42_z) + (gamma31_z+gamma32_z);

%%% Output the derivative of the state
dX(1:3) = [X(4); X(5); X(6)]; % km/s
dX(4:6) = [ddx; ddy; ddz]; % km/s^2

end
