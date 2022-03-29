function [dX] = Int_3BI_J2_testingMeanMotionDerivation_PCIn(t,X,prms)
%%% For numerical integration in the normalized CR3BP with J2 of each body
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

% rp_BC = 


% % % % % % rp = [x; y; z];
% % % % % % 
% % % % % % rps = R3_rot([prms.a;0;0], prms.meanMot*t);
% % % % % % 
% % % % % % rs = rp - rps;
% % % % % % 
% % % % % % ups = up + us;
% % % % % % 
% % % % % % rp_norm = norm(rp);
% % % % % % rps_norm = norm(rps);
% % % % % % rs_norm = norm(rs);
% % % % % % 
% % % % % % % ap = -(up/(rp_norm^3))*rp - ((3/2)*J2p*up*Rp*Rp/(rp_norm^3)).*[ rp(1)/rp_norm - 5*rp(3)*rp(3)*rp(1)/(rp_norm^3); rp(2)/rp_norm - 5*rp(3)*rp(3)*rp(2)/(rp_norm^3); 3*rp(3)/rp_norm - 5*rp(3)*rp(3)*rp(3)/(rp_norm^3)];
% % % % % % % as = -(us/(rs_norm^3))*rs - ((3/2)*J2s*us*Rs*Rs/(rs_norm^3)).*[ rs(1)/rs_norm - 5*rs(3)*rs(3)*rs(1)/(rs_norm^3); rs(2)/rs_norm - 5*rs(3)*rs(3)*rs(2)/(rs_norm^3); 3*rs(3)/rs_norm - 5*rs(3)*rs(3)*rs(3)/(rs_norm^3)];
% % % % % % % aps = (ups/(rps_norm^3))*(1 + 1.5*(J2p*Rp*Rp + J2s*Rs*Rs)/(rps_norm^2)).*rps;
% % % % % % ap = -(up/(rp_norm^3))*rp;
% % % % % % as = -(us/(rs_norm^3))*rs;
% % % % % % aps = (up/(rps_norm^3))*rps;
% % % % % % 
% % % % % % acc = ap + as;
% % % % % % % acc = -us*rs./(rs_norm^3) - up*(rp./(rp_norm^3) - rps./(rps_norm^3));



% 
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
% 
% ap_n = -((1-prms.u)/(rp_norm_n^3)).*rp_n;
% as_n = -((prms.u)/(rs_norm_n^3)).*rs_n;
% aps_n = -prms.u * rps_n;
% 
% acc = ap_n + as_n + aps_n;


%%% PCI_n
r_PCI_n = [x;y;z];

rp_n = r_PCI_n;
rps_n = R3_rot([1;0;0], t);
rs_n = rp_n - rps_n;

rp_norm_n = norm(rp_n);
rs_norm_n = norm(rs_n);

ap_n = -((1-prms.u)/(rp_norm_n^3)).*rp_n;
as_n = -((prms.u)/(rs_norm_n^3)).*rs_n;
aps_n = -prms.u * rps_n;

acc = ap_n + as_n + aps_n;










%%% Output the derivative of the state
dX(1:3) = [xd; yd; zd]; % km/s
dX(4:6) = [acc(1); acc(2); acc(3)]; % km/s^2

% if t > 800000
%     989;
% end
end
