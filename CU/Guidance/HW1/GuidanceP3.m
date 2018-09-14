clear
clc
close all
addpath(genpath('/Users/lukebury/Documents/MATLAB/mbin'))
tic
% ========================================================================
%%% Run/Plot Switches
% ========================================================================

% ========================================================================
%%% Importing Data
% ========================================================================
%%% General data on solar system bodies
bodies = getBodyData();

%%% Color options/schemes
colors = get_colors();

% ========================================================================
%%% Running
% ========================================================================

% -------------------------------------------------
% Hohmann transfer
% -------------------------------------------------
%%% Initial orbit
alt1  = 200; % km
r1 = alt1 + bodies.earth.R; % km
Tp_1 = 2*pi*sqrt(r1^3 / bodies.earth.u);


%%% Final orbit
alt2 = 35786; % km
r2 = alt2 + bodies.earth.R; % km
Tp_2 = 2*pi*sqrt(r2^3 / bodies.earth.u);

%%% Transfer orbit
at = (r1 + r2) / 2;
Tp_t = 2*pi*sqrt(at^3 / bodies.earth.u);

%%% Velocities
vc1   = visviva_v( r1, r1, bodies.earth.u);
v1_t  = visviva_v( r1, at, bodies.earth.u);
v2_t  = visviva_v( r2, at, bodies.earth.u);
vc2   = visviva_v( r2, r2, bodies.earth.u);

dV1 = v1_t - vc1
dV2 = vc2 - v2_t

% -------------------------------------------------
% Propagating standard hohmann transfer
% -------------------------------------------------
%%% Setting Initial states
r0_1 = [0, -r1, 0];
v0_1 = [vc1,  0, 0];
r0_t = [0, -r1, 0];
v0_t = [v1_t, 0, 0];
r0_2 = [0, -r2, 0];
v0_2 = [vc2, 0, 0];

X0_1 = [r0_1, v0_1];
X0_t = [r0_t, v0_t];
X0_2 = [r0_2, v0_2];

%%% Setting time vector
ti = 0;
dt = 1;
time_1 = ti:dt:Tp_1;
time_t = ti:dt:Tp_t / 2; 
time_2 = ti:dt:Tp_2;

%%% Choosing ode45 tolerance
tol = 1e-13;

%%% Setting integrator options
options = odeset('RelTol',tol,'AbsTol',tol);

%%% Propagating trajectories
[time_1, X_1] = ode45(@Int_2BI, time_1, X0_1, options, bodies.earth.u);
[time_t, X_t] = ode45(@Int_2BI, time_t, X0_t, options, bodies.earth.u);
[time_2, X_2] = ode45(@Int_2BI, time_2, X0_2, options, bodies.earth.u);

figure(1); hold all
plot(X_1(:,1),X_1(:,2),'color',colors.std.blue,'linewidth',1.5)
t_nom = plot(X_t(:,1),X_t(:,2),'color',colors.std.black,'linewidth',1.5);
plot(X_2(:,1),X_2(:,2),'color',colors.std.blue,'linewidth',1.5)

plotBodyTexture3(bodies.earth.R, [0,0,0], bodies.earth.img)
PlotBoi2('X', 'Y', 14, 'LaTex')
axis equal
view(0,90)

% -------------------------------------------------
% Creating orbit with error
% -------------------------------------------------
%%% Setting Initial state with 1 degree error and getting new parameters
v0_t_e1 = v0_1 + R3([dV1, 0, 0], -1 * pi/180);
X0_t_e1 = [r0_t, v0_t_e1];

%%% Setting options so the trajectory stops at apoapsis
options_Apsis = odeset('Events',@event_Apsis_2BI,'RelTol',tol,'AbsTol',tol);

%%% Propagating trajectories
[time_t_e1, X_t_e1] = ode45(@Int_2BI, time_t*2, X0_t_e1, options_Apsis, bodies.earth.u);

figure(1); hold all
t_error = plot(X_t_e1(:,1),X_t_e1(:,2),'--','color',colors.std.red,'linewidth',2.5);

legend([t_nom, t_error], 'Perfect Transfer', 'Transfer w/ 1° error')
% -------------------------------------------------
% Comparing orbital elements
% -------------------------------------------------
[a,e,i,raan,w,ta] = RV2COE(X_t(1,1:3),X_t(1,4:6),bodies.earth.u);
OEs_t = [a,e,i,raan,w,ta];

[a,e,i,raan,w,ta] = RV2COE(X_t_e1(1,1:3),X_t_e1(1,4:6),bodies.earth.u);
OEs_t_e1 = [a,e,i,raan,w,ta];

d_a_te1    = OEs_t_e1(1) - OEs_t(1)
d_e_te1    = OEs_t_e1(2) - OEs_t(2)
d_i_te1    = OEs_t_e1(3) - OEs_t(3)
d_raan_te1 = OEs_t_e1(4) - OEs_t(4)
d_w_te1    = OEs_t_e1(5) - OEs_t(5)
d_ta_te1   = OEs_t_e1(6) - OEs_t(6)


% -------------------------------------------------
% Adding dV to orbit with error at apoapsis and finding OEs
% -------------------------------------------------
%%% Finding v-hat of ending of trajectory with error
vHat_e1 = X_t_e1(end,4:6) / norm(X_t_e1(end,4:6));

%%% Finding dV
dv_apo_e1 = vHat_e1 .* dV2;

%%% Creating new state
X_postBurn_e1 = [X_t_e1(end,1:3), X_t_e1(end,4:6) + dv_apo_e1];

%%% Finding orbital elements of final orbits
X_finalNominal = [0, r2, 0, -vc2, 0, 0];

[a,e,i,raan,w,ta] = RV2COE(X_postBurn_e1(1,1:3),X_postBurn_e1(1,4:6),bodies.earth.u);
OEs_fe1 = [a,e,i,raan,w,ta];

[a,e,i,raan,w,ta] = RV2COE(X_finalNominal(1,1:3),X_finalNominal(1,4:6),bodies.earth.u);
OEs_fNom = [a,e,i,raan,w,ta];

fprintf('--------\n')
d_aF1    = OEs_fe1(1) - OEs_fNom(1)
d_eF1    = OEs_fe1(2) - OEs_fNom(2)
d_iF1    = OEs_fe1(3) - OEs_fNom(3)
d_raanF1 = OEs_fe1(4) - OEs_fNom(4)
d_wF1    = OEs_fe1(5) - OEs_fNom(5)
d_taF1   = OEs_fe1(6) - OEs_fNom(6)




% -------------------------------------------------
% Creating orbits with 0.1 (e2) and -0.3 (e3) degree pointing errors
% -------------------------------------------------
%%% Setting Initial state with 0.1 and -0.3 degree error and getting new parameters
v0_t_e2 = v0_1 + R3([dV1, 0, 0], -0.1 * pi/180);
X0_t_e2 = [r0_t, v0_t_e2];

v0_t_e3 = v0_1 + R3([dV1, 0, 0], 0.3 * pi/180);
X0_t_e3 = [r0_t, v0_t_e3];

%%% Propagating trajectories
[time_t_e1, X_t_e2] = ode45(@Int_2BI, time_t*2, X0_t_e2, options_Apsis, bodies.earth.u);

%%% Propagating trajectories
options_Apsis_noStop = odeset('Events',@event_Apsis_noStop_2BI,'RelTol',tol,'AbsTol',tol);
[time_t_e1, X_t_e3, time_event, X_event_e3, index_event] = ode45(@Int_2BI, time_t*1.1, X0_t_e3, options_Apsis_noStop, bodies.earth.u);

t_error2 = plot(X_t_e2(:,1),X_t_e2(:,2),'--','color',colors.std.mag,'linewidth',2.5);
t_error3 = plot(X_t_e3(:,1),X_t_e3(:,2),'--','color',colors.std.maglt,'linewidth',2.5);





%%% Finding v-hat of ending of trajectory with error
vHat_e2 = X_t_e2(end,4:6) / norm(X_t_e2(end,4:6));
vHat_e3 = X_event_e3(2,4:6) / norm(X_event_e3(2,4:6));

%%% Finding dV
dv_apo_e2 = vHat_e2 .* dV2;
dv_apo_e3 = vHat_e3 .* dV2;

%%% Creating new state
X_postBurn_e2 = [X_t_e2(end,1:3), X_t_e2(end,4:6) + dv_apo_e2];
X_postBurn_e3 = [X_event_e3(2,1:3), X_event_e3(2,4:6) + dv_apo_e3];

[a,e,i,raan,w,ta] = RV2COE(X_postBurn_e2(1,1:3),X_postBurn_e2(1,4:6),bodies.earth.u);
OEs_fe2 = [a,e,i,raan,w,ta];

[a,e,i,raan,w,ta] = RV2COE(X_postBurn_e3(1,1:3),X_postBurn_e3(1,4:6),bodies.earth.u);
OEs_fe3 = [a,e,i,raan,w,ta];

fprintf('--------\n')
d_aF2    = OEs_fe2(1) - OEs_fNom(1)
d_eF2    = OEs_fe2(2) - OEs_fNom(2)


fprintf('--------\n')
d_aF3    = OEs_fe3(1) - OEs_fNom(1)
d_eF3    = OEs_fe3(2) - OEs_fNom(2)


% -------------------------------------------------
% Creating nominal orbit with finite burns
% -------------------------------------------------
tMag = 30e-3; % km/s
tBurn1 = (dV1/tMag);
time_burn1 = [0 tBurn1];

[time_fb1, X_fb1] = ode45(@Int_2BI_finite, time_burn1, X0_1, options, bodies.earth.u, tMag);

X0_coast1 = X_fb1(end,:);
time_coast1 = [0 (Tp_t*2)];

[time_coast1, X_coast1] = ode45(@Int_2BI, time_coast1, X0_coast1, options_Apsis, bodies.earth.u);

X0_burn2 = X_coast1(end,:);
tBurn2 = (dV2/tMag);
time_burn2 = [0 tBurn2];

[time_fb2, X_fb2] = ode45(@Int_2BI_finite, time_burn2, X0_burn2, options, bodies.earth.u, tMag);

X0_coast2 = X_fb2(end,:);
time_coast2 = [0 Tp_2*.95];

[time_coast2, X_coast2] = ode45(@Int_2BI, time_coast2, X0_coast2, options, bodies.earth.u);


cp1 = plot(X_coast1(:,1),X_coast1(:,2),'k--','linewidth',2);
cp2 = plot(X_coast2(:,1),X_coast2(:,2),'b--','linewidth',2);
bp1 = plot(X_fb1(:,1),X_fb1(:,2),'c','linewidth',3);
bp2 = plot(X_fb2(:,1),X_fb2(:,2),'c','linewidth',3);
legend([cp1 cp2],'Finite Burn, Transfer','Finte Burn, Final')

uistack([bp1,bp2])


% -------------------------------------------------
% Orbital elements of final orbit from finite burn
% -------------------------------------------------
[a,e,i,raan,w,ta] = RV2COE(X_coast2(1,1:3),X_coast2(1,4:6),bodies.earth.u);
OEs_fb_final = [a,e,i,raan,w,ta];

fprintf('--------\n')
d_aF3    = OEs_fb_final(1) - OEs_fNom(1)
d_eF3    = OEs_fb_final(2) - OEs_fNom(2)
d_wF3    = OEs_fb_final(5) - OEs_fNom(5)
d_tF3    = OEs_fb_final(6) - OEs_fNom(6)






