function [dX] = Int_CR3Bn_ZH(t,X,prms)
%%% For numerical integration in the normalized CR3BP with capabilities for
%%% zonal harmonics terms from either body (p-primary, s-secondary)
%%% Inputs:
%          t - normalized time vector
%          X - initial state [6x1]
%          prms - (u, R1, R2, J2p, J2s, J3p, J3s, J4p, J4s)
%
%
%          u - mass ratio of CR3BP system
%          R1 - radius of primary body
%          R2 - radius of secondary body
%          J2p - J2 of primary body
%          J2s - J2 of secondary body
%          J3p - J3 of primary body
%          J3s - J3 of secondary body
%          J4p - J4 of primary body
%          J4s - J4 of secondary body
%          J5p - J5 of primary body
%          J5s - J5 of secondary body
%          J6p - J6 of primary body
%          J6s - J6 of secondary body
% =======================================================================
%%% Preallocate output as column vector
dX = zeros(6,1);

%%% Unpacking some parameters (for brevity/clarity)
mu = prms.u;

%%% Distances to primary (1) and secondary (2) bodies
r1 = sqrt((X(1)+mu)^2 + X(2)^2 + X(3)^2);
r2 = sqrt((X(1)+mu-1)^2 + X(2)^2 + X(3)^2);

if isfield(prms, 'J2p')
    gamma_J2p = 3*(1-mu)*(prms.R1^2)*prms.J2p*(5*X(3)^2-r1^2) / (2*r1^7);
    az_J2p    = 3*(1-mu)*(prms.R1^2)*prms.J2p*(5*X(3)^2-3*r1^2)*X(3) / (2*r1^7);
else
    gamma_J2p = 0;
    az_J2p    = 0;
end

if isfield(prms, 'J2s')
    gamma_J2s = 3*mu*prms.R2^2*prms.J2s*(5*X(3)^2-r2^2) / (2*r2^7);
    az_J2s    = 3*mu*prms.R2^2*prms.J2s*(5*X(3)^2-3*r2^2)*X(3) / (2*r2^7);
else
    gamma_J2s = 0;
    az_J2s    = 0;
end

if isfield(prms, 'J3p')
    gamma_J3p = 5*(1-mu)*(prms.R1^3)*prms.J3p*X(3)*(7*X(3)^2 - 3*r1^2) / (2*r1^9);
    az_J3p    = (1-mu)*(prms.R1^3)*prms.J3p*(3*r1^4 - 30*X(3)^2*r1^2 + 35*X(3)^4) / (2*r1^9);
else
    gamma_J3p = 0;
    az_J3p    = 0;
end

if isfield(prms, 'J3s')
    gamma_J3s = 5*mu*(prms.R2^3)*prms.J3s*X(3)*(7*X(3)^2 - 3*r2^2) / (2*r2^9);
    az_J3s    = mu*(prms.R2^3)*prms.J3s*(3*r2^4 - 30*X(3)^2*r2^2 + 35*X(3)^4) / (2*r2^9);
else
    gamma_J3s = 0;
    az_J3s    = 0;
end

if isfield(prms, 'J4p')
    gamma_J4p = 15*(1-mu)*(prms.R1^4)*prms.J4p*(21*X(3)^4 - 14*r1^2*X(3)^2 + r1^4) / (8*r1^11);
    az_J4p    = 5*(1-mu)*(prms.R1^4)*prms.J4p*X(3)*(15*r1^4 - 70*X(3)^2*r1^2 + 63*X(3)^4) / (8*r1^11);
else
    gamma_J4p = 0;
    az_J4p    = 0;
end

if isfield(prms, 'J4s')
    gamma_J4s = 15*mu*(prms.R2^4)*prms.J4s*(21*X(3)^4 - 14*r2^2*X(3)^2 + r2^4) / (8*r2^11);
    az_J4s    = 5*mu*(prms.R2^4)*prms.J4s*X(3)*(15*r2^4 - 70*X(3)^2*r2^2 + 63*X(3)^4) / (8*r2^11);
else
    gamma_J4s = 0;
    az_J4s    = 0;
end

if isfield(prms, 'J5p')
    gamma_J5p = 7*(1-mu)*(prms.R1^5)*prms.J5p*X(3)*(15*r1^4 - 90*X(3)^2*r1^2 + 99*X(3)^4) / (8*r1^13);
    az_J5p    = 3*(1-mu)*(prms.R1^5)*prms.J5p*(-5*r1^6 + 105*X(3)^2*r1^4 - 315*X(3)^4*r1^2 + 231*X(3)^6) / (8*r1^13);
else
    gamma_J5p = 0;
    az_J5p    = 0;
end

if isfield(prms, 'J5s')
    gamma_J5s = 7*mu*(prms.R2^5)*prms.J5s*X(3)*(15*r2^4 - 90*X(3)^2*r2^2 + 99*X(3)^4) / (8*r2^13);
    az_J5s    = 3*mu*(prms.R2^5)*prms.J5s*(-5*r2^6 + 105*X(3)^2*r2^4 - 315*X(3)^4*r2^2 + 231*X(3)^6) / (8*r2^13);
else
    gamma_J5s = 0;
    az_J5s    = 0;
end

if isfield(prms, 'J6p')
    gamma_J6p = 7*(1-mu)*(prms.R1^6)*prms.J6p*(-5*r1^6 + 135*X(3)^2*r1^4 - 495*X(3)^4*r1^2 + 429*X(3)^6) / (16*r1^15);
    az_J6p    = 7*(1-mu)*(prms.R1^6)*prms.J6p*X(3)*(-35*r1^6 + 315*X(3)^2*r1^4 - 693*X(3)^4*r1^2 + 429*X(3)^6) / (16*r1^15);
else
    gamma_J6p = 0;
    az_J6p    = 0;
end

if isfield(prms, 'J6s')
    gamma_J6s = 7*mu*(prms.R2^6)*prms.J6s*(-5*r2^6 + 135*X(3)^2*r2^4 - 495*X(3)^4*r2^2 + 429*X(3)^6) / (16*r2^15);
    az_J6s    = 7*mu*(prms.R2^6)*prms.J6s*X(3)*(-35*r2^6 + 315*X(3)^2*r2^4 - 693*X(3)^4*r2^2 + 429*X(3)^6) / (16*r2^15);
else
    gamma_J6s = 0;
    az_J6s    = 0;
end

%%% Equations of Motion
ddx = 2*prms.n*X(5) + (prms.n^2)*X(1) + (-(1-mu)/r1^3 + gamma_J2p + gamma_J3p + gamma_J4p + gamma_J5p + gamma_J6p)*(X(1) + mu) + (-mu/r2^3 + gamma_J2s + gamma_J3s + gamma_J4s + gamma_J5s + gamma_J6s)*(X(1) - 1 + mu);
ddy = -2*prms.n*X(4) + (prms.n^2)*X(2) + X(2)*(-(1-mu)/r1^3 - mu/r2^3 + gamma_J2p + gamma_J2s + gamma_J3p + gamma_J3s + gamma_J4p + gamma_J4s + gamma_J5p + gamma_J5s + gamma_J6p + gamma_J6s);
ddz = X(3)*(-(1-mu)/r1^3 - mu/r2^3) + az_J2p + az_J2s + az_J3p + az_J3s + az_J4p + az_J4s + az_J5p + az_J5s + az_J6p + az_J6s;

%%% Output the derivative of the state
dX(1:3) = [X(4); X(5); X(6)]; 
dX(4:6) = [ddx; ddy; ddz]; 


end
