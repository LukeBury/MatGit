function [dX] = Int_CR3Bn_J2pJ4pJ6pJ2s(t,X,prms)
%%% For numerical integration in the normalized CR3BP with capabilities for
%%% zonal harmonics terms from either body (p-primary, s-secondary)
%%% Inputs:
%          t - normalized time vector
%          X - initial state [6x1]
%          prms - (u, R1_n, R2_n, J2p, J2s, J3p, J3s, J4p, J4s)
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
% =======================================================================
%%% Preallocate state output
dX = zeros(6,1);

%%% Unpacking some parameters (for brevity/clarity)
mu = prms.u;

%%% Distances to primary (1) and secondary (2) bodies
r1 = sqrt((X(1)+mu)^2 + X(2)^2 + X(3)^2);
r2 = sqrt((X(1)+mu-1)^2 + X(2)^2 + X(3)^2);

% if isfield(prms, 'J2p')
    gamma_J2p = 3*(1-mu)*(prms.R1^2)*prms.J2p*(5*X(3)^2-r1^2) / (2*r1^7);
    az_J2p    = 3*(1-mu)*(prms.R1^2)*prms.J2p*(5*X(3)^2-3*r1^2)*X(3) / (2*r1^7);
% else
%     gamma_J2p = 0;
%     az_J2p    = 0;
% end

% if isfield(prms, 'J2s')
    gamma_J2s = 3*mu*prms.R2^2*prms.J2s*(5*X(3)^2-r2^2) / (2*r2^7);
    az_J2s    = 3*mu*prms.R2^2*prms.J2s*(5*X(3)^2-3*r2^2)*X(3) / (2*r2^7);
% else
%     gamma_J2s = 0;
%     az_J2s    = 0;
% end

% if isfield(prms, 'J4p')
    gamma_J4p = 15*(1-mu)*(prms.R1^4)*prms.J4p*(21*X(3)^4 - 14*r1^2*X(3)^2 + r1^4) / (8*r1^11);
    az_J4p    = 5*(1-mu)*(prms.R1^4)*prms.J4p*X(3)*(15*r1^4 - 70*X(3)^2*r1^2 + 63*X(3)^4) / (8*r1^11);
% else
%     gamma_J4p = 0;
%     az_J4p    = 0;
% end

% if isfield(prms, 'J6p')
    gamma_J6p = 7*(1-mu)*(prms.R1^6)*prms.J6p*(-5*r1^6 + 135*X(3)^2*r1^4 - 495*X(3)^4*r1^2 + 429*X(3)^6) / (16*r1^15);
    az_J6p    = 7*(1-mu)*(prms.R1^6)*prms.J6p*X(3)*(-35*r1^6 + 315*X(3)^2*r1^4 - 693*X(3)^4*r1^2 + 429*X(3)^6) / (16*r1^15);
% else
%     gamma_J6p = 0;
%     az_J6p    = 0;
% end


%%% Equations of Motion
ddx = 2*X(5) + X(1) + (-(1-mu)/r1^3 + gamma_J2p + gamma_J4p + gamma_J6p)*(X(1) + mu) + (-mu/r2^3 + gamma_J2s)*(X(1) - 1 + mu);
ddy = -2*X(4) + X(2) + X(2)*(-(1-mu)/r1^3 - mu/r2^3 + gamma_J2p + gamma_J2s + gamma_J4p + gamma_J6p);
ddz = X(3)*(-(1-mu)/r1^3 - mu/r2^3) + az_J2p + az_J2s + az_J4p + az_J6p;

%%% Output the derivative of the state
dX(1:3) = [X(4); X(5); X(6)]; 
dX(4:6) = [ddx; ddy; ddz]; 

end
