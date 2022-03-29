%%% Testground for figuring out primary-centric inclination changes as a
%%% result of a flyby of the secondary
clear
clc
close all
mbinPath = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/mbin';
addpath(genpath(mbinPath))

% ========================================================================
%%% Importing Data
% ========================================================================
%%% General data on solar system bodies
bodies = getBodyData(mbinPath);

%%% Color options/schemes
colors = get_colors();

% ========================================================================
%%% Run Switches
% ========================================================================


% ========================================================================
%%% Setup
% ========================================================================
% n_thetas = 100;
% thetas = linspace(0,360,n_thetas);

% orbitColors = colorScale([colors.std.cyan; colors.std.mag],n_thetas);

% for jj = 1:n_thetas
% theta_FB_deg = thetas(jj);
% ------------------------------------
%%% Preferences
% ------------------------------------
flybyAltitde = 10000; % km
theta_FB_deg = 0; % degrees (angle from -z axis, rotated around vInf)

% ------------------------------------
%%% System
% ------------------------------------
primary   = bodies.earth; 
secondary = bodies.moon;

% primary   = bodies.mars; 
% secondary = bodies.phobos;

% primary   = bodies.jupiter; 
% secondary = bodies.europa;
% ------------------------------------
%%% Initial orbits
% ------------------------------------
%%% Time
t_i = 0; % sec
t_f = secondary.Tp; % Long bc events are watching for impact or escape
n_dt = 100;
time0 = linspace(t_i,t_f,n_dt);

%%% Choosing ode tolerance
tol = 1e-13;

%%% Setting integrator options
options = odeset('RelTol',tol,'AbsTol',tol);

%%% Initial conditions
X0_secondary = [secondary.a; 0; 0; 0; sqrt(primary.u/secondary.a); 0];
[a_secondary,e_secondary,i_secondary,raan_secondary,w_secondary,ta_secondary] = RV2COE(X0_secondary(1:3),X0_secondary(4:6),primary.u);

%%% Integrating
[time, X_secondary_orbit] = ode113(@Int_2BI, time0, X0_secondary, options, primary.u);

% ========================================================================
%%% Study
% ========================================================================
% ------------------------------------
%%% Set up a spacecraft orbit
% ------------------------------------
a0 = secondary.a;
e0 = 0.5;
i0 = 0;
raan0 = 270;
w0 = 0;
ta0 = 0;
% a0 = secondary.a;
% e0 = 0.1;
% i0 = 0;
% raan0 = 270;
% w0 = 0;
% ta0 = 0;
[r0_sc, v0_sc] = COE2RV(a0, e0, i0, raan0, w0, ta0, primary.u);
X0_sc = [r0_sc; v0_sc];
[time, X_sc0_orbit] = ode113(@Int_2BI, time0, X0_sc, options, primary.u);




% ------------------------------------
%%% Find conditions at secondary crossing
% ------------------------------------
p = a0*(1-e0^2);
ta_crossing = acosd(((p/a0)-1)/e0);
[rCross_sc, vCross_sc_m] = COE2RV(a0, e0, i0, raan0, w0, ta_crossing, primary.u);
[rCross_secondary, vCross_secondary] = COE2RV(a0, 0, 0, 270, 0, ta_crossing, primary.u);

[ fpa_rad ] = rv2fpa( rCross_sc', vCross_sc_m');
fpa_deg = fpa_rad*180/pi;

vInf_m = vCross_sc_m - vCross_secondary;
vInf_m_hat = vInf_m./norm(vInf_m);

rInf0_m = [0; 0; -1];
theta_FB_rad = theta_FB_deg*pi/180;
rInf_m = rotVecAroundArbAxis(rInf0_m',vInf_m_hat',theta_FB_rad);
rInf_m = rInf_m';

hInf = cross(rInf_m, vInf_m);
hInf_hat = hInf ./ norm(hInf);

%%% Get rotation matrix from inertial to infinity frame
N_N = eye(3);
I_N = [rInf_m, vInf_m, hInf_hat];
[IN] = getDCM(I_N, N_N);

%%% Rotate spacecraft  into infinity frame
% vCross_sc_InfinityFrame_m = IN * vCross_sc_m; % km/s
vInf_InfinityFrame_m = IN * vInf_m;

%%% Calculate turning angle
flybyRadius = secondary.R + flybyAltitde; % km
[ turningAngle_rad ] = calcTurningAngle( flybyRadius, norm(vInf_m), secondary.u );

%%% Rotate vInf_m about hInfinity
% vCross_sc_InfinityFrame_p = R3(vCross_sc_InfinityFrame_m, turningAngle_rad);
vInf_InfinityFrame_p = R3(vInf_InfinityFrame_m, turningAngle_rad);

%%% Rotate new spacecraft velocity back into inertial frame
NI = IN';
% vCross_sc_p = NI * vCross_sc_InfinityFrame_p;
vInf_p = NI * vInf_InfinityFrame_p;

vCross_sc_p = vCross_secondary + vInf_p;


%%% Calculate new inclination
[a_new,e_new,i_new,raan_new,w_new,ta_new] = RV2COE(rCross_sc, vCross_sc_p, primary.u);

%%% Changes in orbital elements
delta_a = a_new - a0;
delta_e = e_new - e0;
delta_i = i_new - i0;

preFB_energy = 0.5*(norm(vCross_sc_m))^2 - primary.u/norm(rCross_sc);
postFB_energy = 0.5*(norm(vCross_sc_p))^2 - primary.u/norm(rCross_sc);
delta_energy = postFB_energy - preFB_energy;


%%% Time
%%% Plot the orbits
figure(1); hold all
h = gcf;
h.Position = [686 385 560 420];
% figure('position',[686 385 560 420]); hold all
p1 = plot3(X_secondary_orbit(:,1),X_secondary_orbit(:,2),X_secondary_orbit(:,3),'r','linewidth',2);
p2 = plot3(X_sc0_orbit(:,1),X_sc0_orbit(:,2),X_sc0_orbit(:,3),'m','linewidth',2);

if exist('t_f_new') == 1
    t_last = t_f_new;
end
t_f_new = 2*pi*sqrt((a_new^3)/primary.u); 
if imag(t_f_new) ~= 0
    t_f_new = t_last;
end
time0_new = linspace(t_i,t_f_new,n_dt);
X0_postFB = [rCross_sc; vCross_sc_p];
[time_new, X_sc_postFB] = ode113(@Int_2BI, time0_new, X0_postFB, options, primary.u);
p3 = plot3(X_sc_postFB(:,1),X_sc_postFB(:,2),X_sc_postFB(:,3),'--g','linewidth',2);
% plot3(X_sc_postFB(:,1),X_sc_postFB(:,2),X_sc_postFB(:,3),'--','linewidth',2,'color',orbitColors(jj,:))
p4 = plot3([0 rCross_sc(1)], [0 rCross_sc(2)], [0 rCross_sc(3)], 'k');
plotBody2(primary.R,[0,0,0],primary.color,colors.std.black,1,1)
PlotBoi2('$X_n$','$Y_n$',18,'LaTex')
axis equal
legend([p1 p2 p3 p4],'Secondary Orbit','Initial sc Orbit','Post FB sc Orbit','FB Location')

figure('position',[688 7 560 420]); hold all
p1 = quiver3(0, 0, 0, vCross_secondary(1), vCross_secondary(2), vCross_secondary(3),0,'color',colors.std.blue,'linewidth',2);
p2 = quiver3(vCross_secondary(1), vCross_secondary(2), vCross_secondary(3), vInf_m(1), vInf_m(2), vInf_m(3),0,'color',colors.std.purp,'linewidth',2);
p3 = quiver3(0, 0, 0, vCross_sc_m(1), vCross_sc_m(2), vCross_sc_m(3),0,'color',colors.std.orange,'linewidth',2);
p4 = quiver3(vCross_secondary(1), vCross_secondary(2), vCross_secondary(3), vInf_p(1), vInf_p(2), vInf_p(3),0,'--','color',colors.std.purp,'linewidth',2);
p5 = quiver3(0, 0, 0, vCross_sc_p(1), vCross_sc_p(2), vCross_sc_p(3),0,'--','color',colors.std.orange,'linewidth',2);
axis equal
legend([p1 p2 p3 p4 p5],'V_{body}','V_{\infty^-}','V_{sc^-}','V_{\infty^+}','V_{sc^+}','location','best')
PlotBoi3('$X_n$','$Y_n$','$Z_n$',18,'LaTex')


fprintf('Turning angle:\t%1.4f deg\n',turningAngle_rad*180/pi)
fprintf('delta a:\t%1.4f km\n',delta_a)
fprintf('delta e:\t%1.4f\n',delta_e)
fprintf('delta i:\t%1.4f deg\n',delta_i)
fprintf('delta Energy:\t%1.4f km^2/s^2\n',delta_energy)
% theta_FB_deg
% drawnow
% end















