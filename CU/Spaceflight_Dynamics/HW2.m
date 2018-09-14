clear
clc
addpath('../bin/')

d2r = pi/180;
r2d = 180/pi;


%% ========================================================
% #1
%----------------------------------------------------------
% a = 12756;
% e = 0.15;
% u = 398600.4418;
% n = sqrt(u/(a^3));

%%%%%%%%% 1C %%%%%%%%%%%%%%%%%
% EC = 200*d2r; % rad
% TC = E2T(EC,e) * r2d % deg
% MC = (EC - e*sin(EC)) * r2d % deg
% TimeC = (EC - e*sin(EC))/(n*60) % minute

%%%%%%%%% 1D %%%%%%%%%%%%%%%%%
% TD = 90 * d2r;
% ED = T2E(TD,e) * r2d
% MD = (T2E(TD,e) - e * sin(T2E(TD,e)));
% TimeD = MD/(n*60)
% MD = MD * r2d

%%%%%%%%% 1E %%%%%%%%%%%%%%%%%
% ME = 270 * d2r;
% EE = M2E(ME,e);
% TE = E2T(EE,e) * r2d + 360
% EE = EE * r2d
% TimeE = ME/(n*60)

%%%%%%%%% 1F %%%%%%%%%%%%%%%%%
% TimeF = 25*60;
% MF = n*TimeF;
% EF = M2E(MF,e);
% TF = E2T(EF,e) * r2d
% EF = EF * r2d
% MF = MF * r2d

%% ========================================================
% #2
%----------------------------------------------------------
% dt = 65*60 % sec
% rp = 6378 + 321 % km
% ra = 6378 + 551 % km
% v1 = 330 * d2r % True Anomaly, rad
% e = (ra-rp)/(ra+rp)
% a = (rp+ra)/2
% u = 398600.4418
% n = sqrt(u/(a^3))
% Tp = 2*pi/n
% 
% E1 = T2E(v1,e)
% E1 * r2d
% M1 = E1 - e * sin(E1) 
% M1 * r2d
% dtp1 = M1 / n
% dtp2 = dtp1 + dt
% M2 = n * dtp2 - 2*pi
% M2 * r2d
% E2 = M2E(M2,e)
% E2 * r2d
% T2d = E2T(E2,e) * r2d

%% ========================================================
% #3
%----------------------------------------------------------
% % l = launched
% % p = planned
% u = 398600.4418;
% 
% ra_l = 817 + 6378;
% rp_l = 798 + 6378;
% 
% ra_p = 814 + 6378;
% rp_p = 794 + 6378;
% 
% a_l = (ra_l + rp_l)/2;
% a_p = (ra_p + rp_p)/2;
% a_error = abs(a_l - a_p)/a_p * 100;
% 
% e_l = (ra_l - rp_l)/(ra_l + rp_l)
% e_p = (ra_p - rp_p)/(ra_p + rp_p)
% e_error = abs(e_l - e_p)/e_p * 100
% 
% Tp_l = 2 * pi * sqrt((a_l^3)/u)
% Tp_p = 2 * pi * sqrt((a_p^3)/u)
% Tp_error = abs(Tp_l - Tp_p)/Tp_p * 100


%% ========================================================
% #4
% %----------------------------------------------------------
% r = [-5650; -2650; 2850];
% v = [2.415; -7.032; -1.796];
% u = 398600.4418;
% %n = sqrt(u/(a^3));
% R = norm(r)
% V = norm(v)
% 
% %%% a)
% a = 1/(2/R - (V^2)/u) % Rearranging Vis Viva
% 
% %%% b)
% e = ((V^2 - u/R)*r - dot(r,v)*v) / u % from pg 98
% eee = e;
% e = norm(e)
% 
% %%% c)
% ra = a*(1+e)
% rp = a*(1-e)
% ra_alt = ra - 6378
% rp_alt = rp - 6378
% 
% %%% d)
% EnergyPerUnitMass = (V^2)/2 - u/R
% 
% %%% e)
% h = cross(r,v)
% 
% %%% f)
% k = [0;0;1];
% i = acos(dot(k,h)/(norm(k)*norm(h)))
% i * r2d
% 
% %%% g)
% fpa = acos(norm(h)/(R*V))
% fpa * r2d
% 
% fpa2 = asin(dot(r,v)/(R*V)) *r2d
% acos(dot(eee,r)/(e*R)) * r2d

%% ========================================================
% #5
%----------------------------------------------------------
% R = 8200;
% u = 398600.4418;
% want e = 0;
% Vc = sqrt(u/R) % and fpa = 0 deg b/c e = 0


%% ========================================================
% #6
%----------------------------------------------------------
u = 132712440018;
a = 0.387 * 149597870 % AU to km
e = 0.205;
% Tp in years? rp and ra? vel at rp? 
% 1 year = 365.25 days
% 1 day = 86400 seconds
Tp = 2 * pi * sqrt((a^3)/u) / (86400*365.25) % Period in years
rp = a * (1 - e)
ra = a * (1 + e)
v_p = sqrt(2*u/rp - u/a)
% 
% 
% 
