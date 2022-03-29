function [dX] = Int_3BI_J2_testingMeanMotionDerivation(t,X,prms)
%%% 
%%% Inputs:
%          t - normalized time vector
%          X - initial state [6x1]
%          prms - (up, us, Rp, Rs, J2p, J2s, a, meanMot)
% =======================================================================
%%% Preallocate state output
dX = zeros(6,1);

%%% Unpack the primary-centered state vector
x = X(1); y = X(2); z = X(3); % Position
xd = X(4); yd = X(5); zd = X(6); % Velocity


%%% Unpack prms
up      = prms.up;
us      = prms.us;
Rp      = prms.Rp;
Rs      = prms.Rs;
J2p     = prms.J2p;
J2s     = prms.J2s;
a       = prms.a;
meanMot = prms.meanMot;
m1      = prms.m1;
m2      = prms.m2;
mTot = m1+m2;


% %%% BaCI_n
% r_BaCI_n = [x;y;z];
% 
% rp_BaCI_n = R3_rot([-prms.u;0;0], t);
% rs_BaCI_n = R3_rot([1-prms.u;0;0], t);
% 
% rp_n = r_BaCI_n - rp_BaCI_n;
% rs_n = r_BaCI_n - rs_BaCI_n;
% rps_n = rs_BaCI_n - rp_BaCI_n;
% 
% rp_norm_n = norm(rp_n);
% rs_norm_n = norm(rs_n);
% rps_norm_n = norm(rps_n);
% 
% ap_n = -((1-prms.u)/(rp_norm_n^3)).*rp_n;
% as_n = -((prms.u)/(rs_norm_n^3)).*rs_n;
% aps_n = (1-prms.u).*rps_n;
% 
% acc = ap_n + as_n;



% %%% BaCI
% r_BaCI = [x;y;z];
% rp_BaCI = R3_rot([-(a*m2/mTot);0;0], t*meanMot);
% rs_BaCI = R3_rot([a - a*m2/mTot;0;0], t*meanMot);
% rp = r_BaCI - rp_BaCI;
% rs = r_BaCI - rs_BaCI;
% rps = rs_BaCI - rp_BaCI;
% rp_norm  = norm(rp);
% rs_norm  = norm(rs);
% ap = -(up/(rp_norm^3)).*rp;
% as = -(us/(rs_norm^3)).*rs;
% acc = ap + as;


%%% BaCI w/ J2s and J2p
r_BaCI = [x;y;z];
rp_BaCI = R3_rot([-(a*m2/mTot);0;0], t*meanMot);
rs_BaCI = R3_rot([a - a*m2/mTot;0;0], t*meanMot);
rp = r_BaCI - rp_BaCI;
rs = r_BaCI - rs_BaCI;
rps = rs_BaCI - rp_BaCI;
rp_norm  = norm(rp);
rs_norm  = norm(rs);

ap = -(up/(rp_norm^3)).*rp;
as = -(us/(rs_norm^3)).*rs;

a_J2x = (3*Rp*Rp*J2p*up*(5*z*z - rp_norm^2)/(2*(rp_norm^7)))*rp(1) + (3*Rs*Rs*J2s*us*(5*z*z - rs_norm^2)/(2*(rs_norm^7)))*rs(1);
a_J2y = (3*Rp*Rp*J2p*up*(5*z*z - rp_norm^2)/(2*(rp_norm^7)))*rp(2) + (3*Rs*Rs*J2s*us*(5*z*z - rs_norm^2)/(2*(rs_norm^7)))*rs(2);
a_J2z = (3*Rp*Rp*J2p*z*up*(5*z*z - 3*(rp_norm^2))/(2*(rp_norm^7))) + (3*Rs*Rs*J2s*z*us*(5*z*z - 3*(rs_norm^2))/(2*(rs_norm^7)));
acc = ap + as + [a_J2x; a_J2y; a_J2z];




% if isfield(prms, 'J2p')
%     gamma_J2p = 3*(1-mu)*(prms.R1^2)*prms.J2p*(5*X(3)^2-r1^2) / (2*r1^7);
%     az_J2p    = 3*(1-mu)*(prms.R1^2)*prms.J2p*(5*X(3)^2-3*r1^2)*X(3) / (2*r1^7);
% else
%     gamma_J2p = 0;
%     az_J2p    = 0;
% end


%%% Output the derivative of the state
dX(1:3) = [xd; yd; zd]; % km/s
dX(4:6) = [acc(1); acc(2); acc(3)]; % km/s^2

% if t > 800000
%     989;
% end
end
