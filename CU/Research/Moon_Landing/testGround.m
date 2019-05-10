clear
% clc
close all
mbinPath = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/mbin';
moonFuncsPath = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/CU/Research/Moon_Landing/Moon_Landing_funcs';
addpath(genpath(mbinPath))
addpath(genpath(moonFuncsPath))

% ========================================================================
%%% Importing Data
% ========================================================================
%%% General data on solar system bodies
bodies = getBodyData(mbinPath);

%%% Color options/schemes
colors = get_colors();

% ========================================================================
%%% Testing
% ========================================================================


% %%% Choosing ode tolerance
% tol = 1e-11;
% %%% Setting integrator options
% options = odeset('RelTol',tol,'AbsTol',tol);
% t_i = 0; % sec
% t_f = 2*365.242189*86400; % Long bc events are watching for impact or escape
% n_dt = 10000;
% time0 = linspace(t_i,t_f,n_dt);
% 
% X0 = [-1.116247042938662e+08;...
%     9.651787718164873e+07;...
%     -5.115716280996529e+03;...
%     3.056174806192775e+05;...
%     -2.642564712324188e+05;...
%     14.006328899915344];
% [~, X_out] = ode113(@Int_2BI, time0, X0, options, Sun.u);




% JD_departure = 2460545;
% JD_arrival   = 2460919;
% transferTime_sec = (JD_arrival - JD_departure)*86400;
% 
% [L_rad, a, e, i_rad, Omega_rad, Pi_rad, w_rad, M_rad, ta_rad] = ...
%     getPlanetElements_Meeus(JD_departure, 'Earth', 'radians');
% [r_departure, v_departure] = OE2ECI(a, e, i_rad, Omega_rad, w_rad, ta_rad, bodies.sun.u);
% 
% %%% Find Jupiter position at arrival time
% [L_rad, a, e, i_rad, Omega_rad, Pi_rad, w_rad, M_rad, ta_rad] = ...
%     getPlanetElements_Meeus(JD_arrival, 'Venus', 'radians');
% [r_arrival, v_arrival] = OE2ECI(a, e, i_rad, Omega_rad, w_rad, ta_rad, bodies.sun.u);
% 
% % r_departure
% % r_arrival
% 
% nOrbits = 1;
% [V1, V2,exitflag] = lambertSolver(r_departure, r_arrival, transferTime_sec, nOrbits, 0, bodies.sun.u);
% 
% transferType = 3;
% [v1,v2] = lambertTargeting(r_departure, r_arrival, transferTime_sec, transferType, 1, bodies.sun.u, 0);
% 
% %%% Creating time vector
% t0 = 0;             % sec
% tf = transferTime_sec; % sec
% time = linspace(t0,tf,10000);
% 
% 
% %%% Choosing ode tolerance
% tol = 2.22045e-14;
% 
% %%% Setting integrator options
% options = odeset('RelTol',tol,'AbsTol',tol);
% 
% 
% [time, X_Earth] = ode113(@Int_2BI,time, [r_departure; v_departure], options, bodies.sun.u);
% [time, X_Venus] = ode113(@Int_2BI,time, [r_arrival; v_arrival], options, bodies.sun.u);
% 
% [time_sc, X_sc] = ode113(@Int_2BI,time, [r_departure; V1], options, bodies.sun.u);
% 
% [time_sc, X_sc_m] = ode113(@Int_2BI,time, [r_departure; v1], options, bodies.sun.u);
% 
% figure; hold all
% p1 = plot3(X_Earth(:,1),X_Earth(:,2),X_Earth(:,3),'color',colors.std.blue,'linewidth',2);
% plot3(X_Earth(1,1),X_Earth(1,2),X_Earth(1,3),'ro','linewidth',2)
% p2 = plot3(X_Venus(:,1),X_Venus(:,2),X_Venus(:,3),'color',colors.std.orange,'linewidth',2);
% plot3(X_Venus(1,1),X_Venus(1,2),X_Venus(1,3),'go','linewidth',2)
% p3 = plot3(X_sc(:,1),X_sc(:,2),X_sc(:,3),'--k','linewidth',2);
% p4 = plot3(X_sc_m(:,1),X_sc_m(:,2),X_sc_m(:,3),'--c','linewidth',2);
% PlotBoi3('x','y','z',15)
% axis equal
% legend([p1 p2 p3 p4],{'Earth','Venus','Old Lambert','My Lambert'})



primary = bodies.jupiter;
secondary = bodies.europa;

L123 = EquilibriumPoints(secondary.MR,1:3); % [3x3] of L1, L2, L3 normalized BCR coordinates

%%% Normalizing constants
rNorm = secondary.a;         % n <-> km
tNorm = 1/secondary.meanMot; % n <-> sec
vNorm = rNorm / tNorm;       % n <-> km/sec



%%% Choosing ode tolerance
tol = 1e-13;

%%% Setting integrator options
options = odeset('Event',@event_ImpactEscape_CR3Bn, 'RelTol',tol,'AbsTol',tol);

%%% Setting necessary parameters for integration
prms.u = secondary.MR;
prms.R2_n = secondary.R_n;
prms.L1x = L123(1,1);
prms.L2x = L123(2,1);

% X0 = [1.0204617015266166;-0.0012084061620286;0.0000000000000000;-0.0028596241729016;0.0009291482175995;0.0000000000000000]';
% X0 = [1.0204617015266166,-0.0019850725823437,0.0000000000000000,-0.0010976413130697,-0.0030157447223195,0.0098771743854318];
% X0 = [1.020461701526617; 0.000037454199667; 0.000074923811876; -0.010853703861100; 0.001140770244121; 0]';
X0 = [1.020461701526617; 0.000410765286112; 0.003440050933435; -0.003449156947614; -0.001991371692182; 0.012257623746852];
t_f = 12.49307;

t_i = 0; % sec
% t_f = 4*pi; % Long bc events are watching for impact or escape
n_dt = 10000;
time0_n = linspace(t_i,t_f,n_dt);

[time_n, X_BCR_n] = ode113(@Int_CR3Bn, time0_n, X0, options, prms);

figure; hold all
plot3(X_BCR_n(:,1),X_BCR_n(:,2),X_BCR_n(:,3),'r','linewidth',1.5)
[JC_scInitial] = JacobiConstantCalculator(secondary.MR,X0(1:3)',X0(4:6)');
plotCR3BP_YZNeck( JC_scInitial, secondary.MR , 2, 0, prms, colors.std.black, 1.5)
plotCR3BP_YZNeck( JC_scInitial, secondary.MR , 1, 0, prms, colors.std.black, 1.5)
plotCR3BP_Neck(secondary,L123,JC_scInitial,600,200,colors.std.black,1.5)
% plotBody2(secondary.R_n,[1-secondary.MR,0,0],colors.std.blue,colors.std.black,1,0)
plotBodyTexture3(secondary.R_n, [1-secondary.MR, 0, 0], secondary.img)
view(0,90)
axis equal
PlotBoi3('$X_n$','$Y_n$','$Z_n$',18,'LaTex')

