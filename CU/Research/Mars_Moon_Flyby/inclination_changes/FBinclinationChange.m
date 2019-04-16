%%% Studying how inclination wrt the primary body can change as a result of
%%% flying by a secondary body
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
% ------------------------------------
%%% Preferences
% ------------------------------------
flybyAltitude = 1; % km

% ------------------------------------
%%% System
% ------------------------------------
% primary   = bodies.earth; 
% secondary = bodies.moon;

primary   = bodies.mars; 
secondary = bodies.phobos;

% primary   = bodies.jupiter; 
% secondary = bodies.europa;
% ------------------------------------
%%% Initial orbits
% ------------------------------------
%%% Time
t_i = 0; % sec
t_f = secondary.Tp; % Long bc events are watching for impact or escape
n_dt = 10000;
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
[time, X_sc_orbit] = ode113(@Int_2BI, time0, X0_sc, options, primary.u);

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

n_thetas = 100000;
thetas_rad = linspace(0,2*pi,n_thetas);

%%% Preallocating
results_delta_a = zeros(n_thetas,1);
results_delta_e = zeros(n_thetas,1);
results_delta_i = zeros(n_thetas,1);
results_delta_energy = zeros(n_thetas,1);

%%% Loop through theta values (where theta is angle from negative z-axis
%%% rotated about vInf_m
for kk = 1:n_thetas
rInf_m = rotVecAroundArbAxis(rInf0_m',vInf_m_hat',thetas_rad(kk));
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
flybyRadius = secondary.R + flybyAltitude; % km
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

%%% Storing results
results_delta_a(kk) = delta_a;
results_delta_e(kk) = delta_e;
results_delta_i(kk) = delta_i;
results_delta_energy(kk) = delta_energy;

end
thetas_deg = thetas_rad.*(180/pi);

% %%% Plotting 3D pre- and post-flyby, along with secondary's orbit
% figure; hold all
% plot3(X_secondary_orbit(:,1),X_secondary_orbit(:,2),X_secondary_orbit(:,3),'r','linewidth',2)
% plot3(X_sc_orbit(:,1),X_sc_orbit(:,2),X_sc_orbit(:,3),'m','linewidth',2)
% X0_postFB = [rCross_sc; vCross_sc_p];
% [time, X_sc_postFB] = ode113(@Int_2BI, time0, X0_postFB, options, primary.u);
% plot3(X_sc_postFB(:,1),X_sc_postFB(:,2),X_sc_postFB(:,3),'--g','linewidth',2)
% plot3([0 rCross_sc(1)], [0 rCross_sc(2)], [0 rCross_sc(3)], 'k')
% plotBody2(primary.R,[0,0,0],primary.color,colors.std.black,1,1)
% PlotBoi2('$X_n$','$Y_n$',18,'LaTex')
% axis equal


% %%% Velocity triangles
% figure; hold all
% p1 = quiver3(0, 0, 0, vCross_secondary(1), vCross_secondary(2), vCross_secondary(3),0,'color',colors.std.blue,'linewidth',2);
% p2 = quiver3(vCross_secondary(1), vCross_secondary(2), vCross_secondary(3), vInf_m(1), vInf_m(2), vInf_m(3),0,'color',colors.std.purp,'linewidth',2);
% p3 = quiver3(0, 0, 0, vCross_sc_m(1), vCross_sc_m(2), vCross_sc_m(3),0,'color',colors.std.orange,'linewidth',2);
% p4 = quiver3(vCross_secondary(1), vCross_secondary(2), vCross_secondary(3), vInf_p(1), vInf_p(2), vInf_p(3),0,'--','color',colors.std.purp,'linewidth',2);
% p5 = quiver3(0, 0, 0, vCross_sc_p(1), vCross_sc_p(2), vCross_sc_p(3),0,'--','color',colors.std.orange,'linewidth',2);
% axis equal
% legend([p1 p2 p3 p4 p5],'V_{body}','V_{\infty^-}','V_{sc^-}','V_{\infty^+}','V_{sc^+}','location','best')
% PlotBoi3('$X_n$','$Y_n$','$Z_n$',18,'LaTex')


% fprintf('Turning angle:\t%1.4f deg\n',turningAngle_rad*180/pi)
% fprintf('delta a:\t%1.4f km\n',delta_a)
% fprintf('delta e:\t%1.4f\n',delta_e)
% fprintf('delta i:\t%1.4f deg\n',delta_i)
% fprintf('delta Energy:\t%1.4f km^2/s^2\n',delta_energy)


figure('position',[440 25 609 773]); 
subplot(4,1,1); hold all
title(sprintf('Body: %s  ..... Flyby Altitude: %1.1f km  .....  vInf_m: %1.2f km/s', secondary.title, flybyAltitude, norm(vInf_m)))
plot(thetas_deg, results_delta_i,'b','linewidth',2)
plot([180 180], [min(results_delta_i) max(results_delta_i)], '--k','linewidth',2)
PlotBoi2('','$\Delta i$ wrt Primary, $^\circ$',18,'LaTex')
xlim([0 360])
ylim([min(results_delta_i) max(results_delta_i)])

subplot(4,1,2); hold all
plot(thetas_deg, results_delta_a,'b','linewidth',2)
plot([180 180], [min(results_delta_a) max(results_delta_a)], '--k','linewidth',2)
PlotBoi2('','$\Delta a$ wrt Primary, $^\circ$',18,'LaTex')
xlim([0 360])
ylim([min(results_delta_a) max(results_delta_a)])

subplot(4,1,3); hold all
plot(thetas_deg, results_delta_e,'b','linewidth',2)
plot([180 180], [min(results_delta_e) max(results_delta_e)], '--k','linewidth',2)
PlotBoi2('','$\Delta e$ wrt Primary, $^\circ$',18,'LaTex')
xlim([0 360])
ylim([min(results_delta_e) max(results_delta_e)])

subplot(4,1,4); hold all
plot(thetas_deg, results_delta_energy,'b','linewidth',2)
plot([180 180], [min(results_delta_energy) max(results_delta_energy)], '--k','linewidth',2)
PlotBoi2('$\theta$, $^\circ$','$\Delta E$ wrt Primary, $^\circ$',18,'LaTex')
xlim([0 360])
ylim([min(results_delta_energy) max(results_delta_energy)])









