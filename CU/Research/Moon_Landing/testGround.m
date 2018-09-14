clear
clc
close all
addpath(genpath('/Users/CU_Google_Drive/lukebury/Documents/MatGit/mbin'))

% ========================================================================
%%% Importing Data
% ========================================================================
%%% General data on solar system bodies
bodies = getBodyData();

%%% Color options/schemes
colors = get_colors();

% ========================================================================
%%% Converting Schaub's zonal harmonics accelerations to nicer equations
% ========================================================================
% syms J2 J3 J4 u R r x y z rs1 rs2 r1 r2 real
% % J2 = .3;
% % u = .2;
% % R = .6;
% % x = 2;
% % y = 0;
% % z = 0;
%  
% % r = sqrt(x^2 + y^2 + z^2);
% % rs1 = x+u;     % for sch_x1
% % rs2 = x-(1-u); % for sch_x2
% % r1 = sqrt((-u-x)^2 + y^2 + z^2);
% % r2 = sqrt(((1-u)-x)^2 + y^2 + z^2);
% 
% % % schx  = -(3/2)*J*u*(R*R/(r^4))*(1-5*(z/r)^2)*x/r
% sch_x1  = -(3/2)*J2*(1-u)*(R*R/(rs1^4))*(1-5*(z/rs1)^2)*(x+u)/rs1;
% % sch_x2  = -(3/2)*J2*u*(R*R/(rs2^4))*(1-5*(z/rs2)^2)*(x-(1-u))/rs2
% % 
% % 
% my_x1  = -3*R*R*J2*(1-u)*(-u-x)*(5*z^2-r1^2)/(2*r1^7);
% % my_x2 = 3*R*R*J2*u*((1-u)-x)*(5*z^2-r2^2)/(2*r2^7)
% % % eq2z = 
% 
% aj2_x = -(3/2)*J2*u*(R^2/r^4)*(1-5*(z/r)^2)*(x/r);
% aj2_y = -(3/2)*J2*u*(R^2/r^4)*(1-5*(z/r)^2)*(y/r);
% aj2_z = -(3/2)*J2*u*(R^2/r^4)*(3-5*(z/r)^2)*(z/r);
% 
% aj3_x = -.5*J3*u*(R*R*R/r^5)*(5*(3*(z/r)-7*(z/r)^3)*(x/r));
% aj3_y = -.5*J3*u*(R*R*R/r^5)*(5*(3*(z/r)-7*(z/r)^3)*(y/r));
% aj3_z = -.5*J3*u*(R*R*R/r^5)*(3*(-1+10*(z/r)^2-(35/3)*(z/r)^4));
% 
% aj4_x = -(5/8)*J4*u*(R^4/r^6)*(-3+42*(z/r)^2-63*(z/r)^4)*(x/r);
% aj4_y = -(5/8)*J4*u*(R^4/r^6)*(-3+42*(z/r)^2-63*(z/r)^4)*(y/r);
% aj4_z = -(5/8)*J4*u*(R^4/r^6)*(-15+70*(z/r)^2-63*(z/r)^4)*(z/r);


% % ========================================================================
% %%% J2 vs J4?
% % ========================================================================
% %%% Selecting bodies
% primary   = bodies.saturn;
% secondary = bodies.enceladus;
% 
% %%% Normalizing constants
% rNorm = secondary.a;         % n <-> km
% tNorm = 1/secondary.meanMot; % n <-> sec
% vNorm = rNorm / tNorm;       % n <-> km/sec

% % %%% Acquire Lagrange points
% % L123_noJ2n = EquilibriumPoints(secondary.MR,1:3); % [3x3] of L1, L2, L3 normalized BCR coordinates
% % L123_J2n = EquilibriumPoints_J2(secondary.MR,primary.J2,0,primary.R/rNorm,secondary.R_n);
% % L123 = [L123_noJ2n; L123_J2n]; % stacking
% 
% %%% Initial States (w/o J2)
% x0_n = [1-secondary.MR+4*secondary.R_n;0;0];
% xd0_n = [0; 0.01; 0.02];
% X0_n = [x0_n xd0_n];
% 
% %%% Setting time vector and normalizing 
% ti = 0; % sec
% dt = 10;
% tf = 3600*20;
% time0_n = [ti:dt:tf] ./ tNorm;
% 
% %%% Choosing ode45 tolerance and setting options
% tol = 1e-10;
% 
% % options = odeset('Events',@impactEvent_CR3Bn,'RelTol',tol,'AbsTol',tol);
% options = odeset('RelTol',tol,'AbsTol',tol);
% 
% %%% Propagating the States normally
% [time_n, X_BCR_n] = ode45(@Int_CR3Bn, time0_n, X0_n, options, secondary.MR, secondary.R_n);
%     
% %%% Propagating the States with J2
% [time_n, X_BCR_J2_n] = ode45(@Int_CR3Bn_J2, time0_n, X0_n, options, secondary.MR, primary.R/rNorm, secondary.R_n, primary.J2, 0);
% 
% %%% Propagating the states with J4
% [time_n, X_BCR_J4_n] = ode45(@Int_CR3Bn_perturb, time0_n, X0_n, options, secondary.MR,primary.R/rNorm,secondary.R_n,0, 0, 0, 0, primary.J4, 0);
% 
% %%% Propagating the States with J2 and J4
% [time_n, X_BCR_J24_n] = ode45(@Int_CR3Bn_perturb, time0_n, X0_n, options, secondary.MR,primary.R/rNorm,secondary.R_n,primary.J2, 0, 0, 0, primary.J4, 0);
% 
% 
% %%% Plotting
% figure; hold all
% % plotBody2(secondary.R_n,[1-secondary.MR,0,0],[1,1,1],[0,0,0],1.5,0.0)
% plotBody3(secondary.R_n,[1-secondary.MR,0,0],secondary.color,0.3)
% plot3(X_BCR_n(:,1),X_BCR_n(:,2),X_BCR_n(:,3),'k')
% plot3(X_BCR_J2_n(:,1),X_BCR_J2_n(:,2),X_BCR_J2_n(:,3),'k')
% axis equal
% 
% figure
% subplot(3,1,1); hold all
% plot(time_n.*tNorm./3600,(X_BCR_n(:,1) - X_BCR_J2_n(:,1)).*rNorm.*1000)
% title('Nom vs J2')
% PlotBoi2('','x, m',12)
% subplot(3,1,2); hold all
% plot(time_n.*tNorm./3600,(X_BCR_n(:,2) - X_BCR_J2_n(:,2)).*rNorm.*1000)
% PlotBoi2('','y, m',12)
% subplot(3,1,3); hold all
% plot(time_n.*tNorm./3600,(X_BCR_n(:,3) - X_BCR_J2_n(:,3)).*rNorm.*1000)
% PlotBoi2('Time, hr','z, m',12)
% 
% figure
% subplot(3,1,1); hold all
% plot(time_n.*tNorm./3600,(X_BCR_n(:,1) - X_BCR_J4_n(:,1)).*rNorm.*1000)
% title('Nom vs J4')
% PlotBoi2('','x, m',12)
% subplot(3,1,2); hold all
% plot(time_n.*tNorm./3600,(X_BCR_n(:,2) - X_BCR_J4_n(:,2)).*rNorm.*1000)
% PlotBoi2('','y, m',12)
% subplot(3,1,3); hold all
% plot(time_n.*tNorm./3600,(X_BCR_n(:,3) - X_BCR_J4_n(:,3)).*rNorm.*1000)
% PlotBoi2('Time, hr','z, m',12)
% 
% 
% figure
% subplot(3,1,1); hold all
% plot(time_n.*tNorm./3600,(X_BCR_J2_n(:,1) - X_BCR_J24_n(:,1)).*rNorm.*1000)
% title('J2 vs J24')
% PlotBoi2('','x, m',12)
% subplot(3,1,2); hold all
% plot(time_n.*tNorm./3600,(X_BCR_J2_n(:,2) - X_BCR_J24_n(:,2)).*rNorm.*1000)
% PlotBoi2('','y, m',12)
% subplot(3,1,3); hold all
% plot(time_n.*tNorm./3600,(X_BCR_J2_n(:,3) - X_BCR_J24_n(:,3)).*rNorm.*1000)
% PlotBoi2('Time, hr','z, m',12)











%%% Selecting bodies
primary   = bodies.sun;
secondary = bodies.mercury;

%%% Normalizing constants
rNorm = secondary.a;         % n <-> km
tNorm = 1/secondary.meanMot; % n <-> sec
vNorm = rNorm / tNorm;       % n <-> km/sec

%%% Acquire Lagrange points
L123_noJ2n = EquilibriumPoints(secondary.MR,1:3); % [3x3] of L1, L2, L3 normalized BCR coordinates
L123_J2n = EquilibriumPoints_J2(secondary.MR,primary.J2,0,primary.R/rNorm,secondary.R/rNorm);
% L123_J22n = EquilibriumPoints_J2(secondary.MR,primary.J2,secondary.J2,primary.R/rNorm,secondary.R/rNorm);
% L123 = [L123_noJ2n; L123_J2n]; % stacking

L123_noJ2n = L123_noJ2n.*rNorm;
L123_J2n   = L123_J2n.*rNorm;
% L123_J22n  = L123_J22n.*rNorm;








    