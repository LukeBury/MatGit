% ========================================================================
%%% Description
% ========================================================================
% For propagating manifolds of Cordelia L2 and Ophelia L1 Lyapunov
% manifolds in the bi-circular (4B) model

% Created: 10/15/19
% Author : Luke Bury, luke.bury@colorado.edu
% ========================================================================
%%% Initialization
% ========================================================================
clear
clc
close all
mbinPath = '~/CU_Google_Drive/Documents/MatGit/mbin';
manifoldPath = '~/CU_Google_Drive/Documents/MatGit/CU/Research/UranusRing/uranus_moon_manifolds/';
addpath(genpath(mbinPath))
ticWhole = tic;

% ========================================================================
%%% Importing Data
% ========================================================================
%%% General data on solar system bodies
bodies = getBodyData(mbinPath);

%%% Color options/schemes
colors = get_colors();

%%% Periodic orbit ICs
PO_ICs = get_PO_ICs();

% ========================================================================
%%% Run Switches
% ========================================================================


% ========================================================================
%%% Setup
% ========================================================================
% -------------------------------------------------
%%% Bodies and normalizations
% -------------------------------------------------
%%% Important bodies
uranus   = bodies.uranus;
cordelia = bodies.cordelia;
ophelia  = bodies.ophelia;

%%% Normalization factors
rNorm_cor = cordelia.a;            % n <-> km
tNorm_cor = 1 / cordelia.meanMot;  % n <-> sec
vNorm_cor = rNorm_cor / tNorm_cor; % n <-> km/sec

rNorm_oph = ophelia.a;             % n <-> km
tNorm_oph = 1 / ophelia.meanMot;   % n <-> sec
vNorm_oph = rNorm_oph / tNorm_oph; % n <-> km/sec

% -------------------------------------------------
%%% Equillibrium Points
% -------------------------------------------------
rLPs_cor_n = EquilibriumPoints(cordelia.MR);
rLPs_oph_n  = EquilibriumPoints(ophelia.MR);

% -------------------------------------------------
%%% Load manifold X0s
% -------------------------------------------------
X0s_cor_L2_stable_lowEnergy  = dlmread([manifoldPath,'X0s_cor_L2_Lyap_stable_lowEnergy.txt'], ',', 1, 0);
X0s_cor_L2_stable_highEnergy = dlmread([manifoldPath,'X0s_cor_L2_Lyap_stable_highEnergy.txt'], ',', 1, 0);

X0s_oph_L1_stable_lowEnergy  = dlmread([manifoldPath,'X0s_oph_L1_Lyap_stable_lowEnergy.txt'], ',', 1, 0);
X0s_oph_L1_stable_highEnergy = dlmread([manifoldPath,'X0s_oph_L1_Lyap_stable_highEnergy.txt'], ',', 1, 0);

% -------------------------------------------------
%%% Integrator options
% -------------------------------------------------
%%% Tolerance
tol = 1e-13;

%%% Options
options = odeset('RelTol',tol,'AbsTol',tol);


% ========================================================================
%%% Integrating  Manifolds
% ========================================================================
mu_sys_cor = cordelia.MR;
% mu3_cor    = ophelia.mass / (uranus.mass + cordelia.mass);
mu3_cor    = 0;
r13_cor    = ophelia.a / rNorm_cor;
theta0_cor = pi;
T2_cor     = cordelia.T_sec / tNorm_cor;
T3_cor     = ophelia.T_sec / tNorm_cor;

mu_sys_oph = ophelia.MR;
% mu3_oph    = cordelia.mass / (uranus.mass + ophelia.mass);
mu3_oph    = 0;
r13_oph = cordelia.a / rNorm_oph;
theta0_oph = pi;
T2_oph = ophelia.T_sec / tNorm_oph;
T3_oph = cordelia.T_sec / tNorm_oph;


% Tf_man_n    = 8000*2*pi;
% Tf_man_n    = 16000*2*pi;
Tf_man_n    = 400000*2*pi;

% timeVec_bkd = linspace(Tf_cor_n, 0, 513013);
timeVec_bkd_man = [Tf_man_n, 0];




n_manifolds = size(X0s_cor_L2_stable_lowEnergy,1);
% n_manifolds = 1;


trajs_cor = cell(n_manifolds,1);
trajs_oph = cell(n_manifolds,1);

parfor kk = 1:10
    dd = kk *5;
    [~, X_cor_4B_lowEnergy] = ode113(@xdot_CR4BP, timeVec_bkd_man, X0s_cor_L2_stable_lowEnergy(dd,:)', options, mu_sys_cor, mu3_cor, r13_cor, theta0_cor, T2_cor, T3_cor);
    [~, X_cor_4B_highEnergy] = ode113(@xdot_CR4BP, timeVec_bkd_man, X0s_cor_L2_stable_highEnergy(dd,:)', options, mu_sys_cor, mu3_cor, r13_cor, theta0_cor, T2_cor, T3_cor);
    
    [~, X_oph_4B_lowEnergy] = ode113(@xdot_CR4BP, timeVec_bkd_man, X0s_oph_L1_stable_lowEnergy(dd,:)', options, mu_sys_oph, mu3_oph, r13_oph, theta0_oph, T2_oph, T3_oph);
    [~, X_oph_4B_highEnergy] = ode113(@xdot_CR4BP, timeVec_bkd_man, X0s_oph_L1_stable_highEnergy(dd,:)', options, mu_sys_oph, mu3_oph, r13_oph, theta0_oph, T2_oph, T3_oph);
    
    trajs_cor{kk}.lowEnergy = X_cor_4B_lowEnergy(:,1:3);
    trajs_cor{kk}.highEnergy = X_cor_4B_highEnergy(:,1:3);
    
    trajs_oph{kk}.lowEnergy = X_oph_4B_lowEnergy(:,1:3);
    trajs_oph{kk}.highEnergy = X_oph_4B_highEnergy(:,1:3);
end

%%% Restructuring trajectory data from parfor - now one large
%%% struct array
trajs_cor_structArray = [trajs_cor{:}];
trajs_oph_structArray = [trajs_oph{:}];

Xs_cor_4B_lowEnergy_n  = [];
Xs_cor_4B_highEnergy_n = [];
Xs_oph_4B_lowEnergy_n  = [];
Xs_oph_4B_highEnergy_n = [];

for kk = 1:10
    Xs_cor_4B_lowEnergy_n  = [Xs_cor_4B_lowEnergy_n; trajs_cor_structArray(kk).lowEnergy; NaN(1,3)];
    Xs_cor_4B_highEnergy_n = [Xs_cor_4B_highEnergy_n; trajs_cor_structArray(kk).highEnergy; NaN(1,3)];
    
    Xs_oph_4B_lowEnergy_n  = [Xs_oph_4B_lowEnergy_n; trajs_oph_structArray(kk).lowEnergy; NaN(1,3)];
    Xs_oph_4B_highEnergy_n = [Xs_oph_4B_highEnergy_n; trajs_oph_structArray(kk).highEnergy; NaN(1,3)];
end

figure; hold all
plot3(Xs_cor_4B_lowEnergy_n(:,1).*rNorm_cor,Xs_cor_4B_lowEnergy_n(:,2).*rNorm_cor,Xs_cor_4B_lowEnergy_n(:,3).*rNorm_cor,'r')
plot3(Xs_oph_4B_lowEnergy_n(:,1).*rNorm_oph,Xs_oph_4B_lowEnergy_n(:,2).*rNorm_oph,Xs_oph_4B_lowEnergy_n(:,3).*rNorm_oph,'r')
PlotBoi3_CR3Bn(20)

figure; hold all
plot3(Xs_cor_4B_highEnergy_n(:,1).*rNorm_cor,Xs_cor_4B_highEnergy_n(:,2).*rNorm_cor,Xs_cor_4B_highEnergy_n(:,3).*rNorm_cor,'r')
plot3(Xs_oph_4B_highEnergy_n(:,1).*rNorm_oph,Xs_oph_4B_highEnergy_n(:,2).*rNorm_oph,Xs_oph_4B_highEnergy_n(:,3).*rNorm_oph,'r')
PlotBoi3_CR3B_km(20)











% -------------------------------------------------
%%% Data from paper
% -------------------------------------------------
% Mean distance: 51149.3 km
% velocity: 10.650 km/s
% Period: 8.3823 hr
% width: 20 to 96 km
% e: 0.00794
% i: 0
a_ring_km        = 51149.3; % km
e_ring           = 0.00794; % 
minWidth_ring_km = 26;      % km
maxWidth_ring_km = 96;      % km
% -------------------------------------------------
%%% Checking minimum and maximum radius of ring particles
% -------------------------------------------------
%%% Calculating periapsis and apoapsis distances
rp_ring_km = a_ring_km * (1 - e_ring); % km
ra_ring_km = a_ring_km * (1 + e_ring); % km

%%% Calculating minimum and maximum theoretical radius for ring
minRadius_ring_km = rp_ring_km - maxWidth_ring_km; % km
maxRadius_ring_km = ra_ring_km + maxWidth_ring_km; % km
r = minRadius_ring_km;
R = maxRadius_ring_km;
xf = 0;
Xf = 0;
yf = 0;
Yf = 0;
% Creates a annulus patch object and returns the handle.  Input arguments 
% are the inner radius, outer radius, inner x offset, outer x offset, inner
% y offset and outer Y offset.  Changes to the edgecolor and linestyle are
% allowed, and will preserve the correct look of the annulus
t = linspace(0,2*pi,200);
x = xf + r*cos(t);
y = yf + r*sin(t);
X = Xf + R*cos(t);
Y = Yf + R*sin(t);
P = patch([x X],[y Y],colors.std.grey,'linestyle','non','facealph',.5);
L(1) = line(x,y,'color',colors.std.white);
L(2) = line(X,Y,'color',colors.std.white);

% ========================================================================
%%% Formatting Structures
% ========================================================================
% -------------------------------------------------
%%% 
% -------------------------------------------------

% --------------------------
%%% 
% --------------------------

% ========================================================================
%%% Closeout
% ========================================================================
tocWhole = toc(ticWhole);
fprintf('\n(Elapsed time: %1.4f seconds)\n',tocWhole)
















