function [dX] = Int_CR3BnSTM_ZH_new(t,X,prms)
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
dX = zeros(42,1);

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
ddx = 2*X(5) + X(1) + (-(1-mu)/r1^3 + gamma_J2p + gamma_J3p + gamma_J4p + gamma_J5p + gamma_J6p)*(X(1) + mu) + (-mu/r2^3 + gamma_J2s + gamma_J3s + gamma_J4s + gamma_J5s + gamma_J6s)*(X(1) - 1 + mu);
ddy = -2*X(4) + X(2) + X(2)*(-(1-mu)/r1^3 - mu/r2^3 + gamma_J2p + gamma_J2s + gamma_J3p + gamma_J3s + gamma_J4p + gamma_J4s + gamma_J5p + gamma_J5s + gamma_J6p + gamma_J6s);
ddz = X(3)*(-(1-mu)/r1^3 - mu/r2^3) + az_J2p + az_J2s + az_J3p + az_J3s + az_J4p + az_J4s + az_J5p + az_J5s + az_J6p + az_J6s;

%%% Output the derivative of the state
dX(1:3) = [X(4); X(5); X(6)]; 
dX(4:6) = [ddx; ddy; ddz]; 






%%% Reshaping STM
stm = reshape(X(7:end),6,6);

%%% Evaluate A matrix
A = zeros(6,6);
A(1:3,4:6) = eye(3);
A(4,5) = 2;
A(5,4) = -2;


warning('Not finished - all these terms still need to be defined')
A(4,1) = daxdx_3b + daxdx_J2p + daxdx_J2s + daxdx_J3p + daxdx_J3s + daxdx_J4p + daxdx_J4s + daxdx_J5p + daxdx_J5s + daxdx_J6p + daxdx_J6s; % dxdd/dx
A(4,2) = daxdy_3b + daxdy_J2p + daxdy_J2s + daxdy_J3p + daxdy_J3s + daxdy_J4p + daxdy_J4s + daxdy_J5p + daxdy_J5s + daxdy_J6p + daxdy_J6s; % dxdd/dy
A(4,3) = daxdz_3b + daxdz_J2p + daxdz_J2s + daxdz_J3p + daxdz_J3s + daxdz_J4p + daxdz_J4s + daxdz_J5p + daxdz_J5s + daxdz_J6p + daxdz_J6s; % dxdd/dz

A(5,1) = daydx_3b + daydx_J2p + daydx_J2s + daydx_J3p + daydx_J3s + daydx_J4p + daydx_J4s + daydx_J5p + daydx_J5s + daydx_J6p + daydx_J6s; % dydd/dx
A(5,2) = daydy_3b + daydy_J2p + daydy_J2s + daydy_J3p + daydy_J3s + daydy_J4p + daydy_J4s + daydy_J5p + daydy_J5s + daydy_J6p + daydy_J6s; % dydd/dy
A(5,3) = daydz_3b + daydz_J2p + daydz_J2s + daydz_J3p + daydz_J3s + daydz_J4p + daydz_J4s + daydz_J5p + daydz_J5s + daydz_J6p + daydz_J6s; % dydd/dz

A(6,1) = dazdx_3b + dazdx_J2p + dazdx_J2s + dazdx_J3p + dazdx_J3s + dazdx_J4p + dazdx_J4s + dazdx_J5p + dazdx_J5s + dazdx_J6p + dazdx_J6s; % dzdd/dx
A(6,2) = dazdy_3b + dazdy_J2p + dazdy_J2s + dazdy_J3p + dazdy_J3s + dazdy_J4p + dazdy_J4s + dazdy_J5p + dazdy_J5s + dazdy_J6p + dazdy_J6s; % dzdd/dy
A(6,3) = dazdz_3b + dazdz_J2p + dazdz_J2s + dazdz_J3p + dazdz_J3s + dazdz_J4p + dazdz_J4s + dazdz_J5p + dazdz_J5s + dazdz_J6p + dazdz_J6s; % dzdd/dz


%%% Calculate stmDot and output result
stmDot = A*stm;
dX(7:end) = reshape(stmDot,36,1);



end
