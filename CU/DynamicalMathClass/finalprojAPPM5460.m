clear
clc
close all
addpath(genpath('/Users/lukebury/Documents/MATLAB/mbin'))

% ========================================================================
%%% Importing Data
% ========================================================================
%%% General data on solar system bodies
bodies = getBodyData();

%%% Color options/schemes
colors = get_colors();

% ========================================================================
%%% Setup
% ========================================================================
primary = bodies.jupiter;
secondary = bodies.europa;

%%% Normalizing constants
rNorm = secondary.a;         % n <-> km
tNorm = 1/secondary.meanMot; % n <-> sec
vNorm = rNorm / tNorm;       % n <-> km/sec

%%% cirucular orbit velocity
rC = secondary.R + 300;
vC = sqrt(secondary.u/rC);
rC_n = rC./rNorm;
vC_n = vC./vNorm;

x0 = [rC,0,0];
xdot0 = [0, vC, 0];
X0 = [x0 xdot0]';

x0_n = [(1-secondary.MR) + rC_n,0,0];
xdot0_n = [0, vC_n, 0];
X0_n = [x0_n xdot0_n]';

% ========================================================================
%%% Integration
% ========================================================================
%%% Setting time vector and normalizing 
ti = 0; % sec
dt = 1; % sec
tf = 3600*24*0.2; % sec
time0 = [ti:dt:tf];
time0_n = time0 ./ tNorm;

secondary.MR = 2.5e-04;

%%% Choosing ode45 tolerance
tol = 1e-9;

%%% Setting integrator options
options = odeset('RelTol',tol,'AbsTol',tol);

%%% Propagating the States in 2B
[time_n, X_2BI] = ode45(@Int_2BI , time0, X0, options, secondary.u);
X_2BI_n = [X_2BI(:,1:3)./rNorm, X_2BI(:,4:6)./vNorm];

%%% Propagating the States in 3B
[time_n, X_BCR_n] = ode45(@Int_CR3Bn, time0_n, X0_n, options, secondary.MR, secondary.R_n);

%%% Plotting 2B
figure('pos',[10 10 700 350]); subplot(1,2,1); hold all
plotBody2(secondary.R_n,[1-secondary.MR,0,0],[1,1,1],[0,0,0],3)
plot3(nan,nan,nan,'color',[0,0,0],'linewidth',1.5);
plot3(X_2BI_n(:,1)+1-secondary.MR,X_2BI_n(:,2),X_2BI_n(:,3),'color',colors.std.purp,'linewidth',1.5);
plot3(X_2BI_n(1,1)+1-secondary.MR,X_2BI_n(1,2),X_2BI_n(1,3),'rx','linewidth',2,'markersize',12);
PlotBoi2('y$_{n}$','x$_{n}$',18,'Latex')
axis equal

%%% Plotting 2B
subplot(1,2,2); hold all
plotBody2(secondary.R_n,[1-secondary.MR,0,0],[1,1,1],[0,0,0],3)
p1 = plot3(nan,nan,nan,'color',[0,0,0],'linewidth',1.5);
p2 = plot3(X_BCR_n(:,1),X_BCR_n(:,2),X_BCR_n(:,3),'color',colors.std.purp,'linewidth',1.5);
p3 = plot3(X_BCR_n(1,1),X_BCR_n(1,2),X_BCR_n(1,3),'rx','linewidth',2,'markersize',12);
legend([p1,p2,p3],'Secondary Surface','Trajactory','Start Point')
PlotBoi2('y$_{n}$','x$_{n}$',18,'Latex')
axis equal

% figure; hold all
% plot3(X_BCR_n(:,1),X_BCR_n(:,2),X_BCR_n(:,3),'b','linewidth',2)                       % trajectory
% plot3(X_BCR_n(1,1),X_BCR_n(1,2),X_BCR_n(1,3),'bo','linewidth',2,'markersize',8)       % starting point
% plot3(X_BCR_n(end,1),X_BCR_n(end,2),X_BCR_n(end,3),'bx','linewidth',2,'markersize',8) % ending point
% plot3(L123_n(1,1),L123_n(1,2),L123_n(1,3),'rx','linewidth',2,'markersize',8)          % L1
% plot3(L123_n(2,1),L123_n(2,2),L123_n(2,3),'rx','linewidth',2,'markersize',8)          % L2
% plotBodyTexture3(secondary.R_n, [1-secondary.MR, 0, 0],secondary.img);
% % plotBody3(secondary.R_n, [1-secondary.MR, 0, 0], secondary.color)
% PlotBoi3('Normalized X Position','Normalized X Position','Z_n',12)
% axis equal
% view(0,90)
% grid on


