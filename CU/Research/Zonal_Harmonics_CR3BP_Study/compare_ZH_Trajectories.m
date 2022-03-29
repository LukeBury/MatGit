% ========================================================================
%%% Description
% ========================================================================
% For comparing CR3BP trajectories with and without zonal harmonics

% Created: 12/11/19
% Author : Luke Bury, luke.bury@colorado.edu
% ========================================================================
%%% Initialization
% ========================================================================
clear
clc
close all
mbinPath = '~/CU_Google_Drive/Documents/MatGit/mbin';
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
%%% Selecting system
% -------------------------------------------------
% sysName = 'Earth_Moon';
% sysName = 'Mars_Phobos';
% sysName = 'Jupiter_Io';
% sysName = 'Jupiter_Europa';
% sysName = 'Jupiter_Ganymede';
% sysName = 'Jupiter_Callisto';
sysName = 'Saturn_Enceladus';
% sysName = 'Saturn_Titan';
% sysName = 'Neptune_Triton';

% -------------------------------------------------
%%% Assigning parameters
% -------------------------------------------------
%%% Set primary & secondary
[primary, secondary] = assignPrimaryAndSecondary_CR3BP(sysName, bodies);

%%% Factor for normalizing distances
rNorm = secondary.a; % n <-> km
tNorm = 1/secondary.meanMot; % n <-> sec
vNorm = rNorm / tNorm;       % n <-> km/sec

%%% Setting each parameters structure
prms.u = secondary.MR;
prms.R1 = primary.R / rNorm;
prms.R2 = secondary.R_n;

prms_ZH     = prms;
prms_ZH.J2p = primary.J2;
prms_ZH.J4p = primary.J4;
prms_ZH.J6p = primary.J6;
prms_ZH.J2s = secondary.J2;

prms.n = 1;
prms_ZH.n = 1;

%%% Equillibrium Points
rLPs_n    = collinearEquilibriumPoints_ZH(prms);
rLPs_ZH_n = collinearEquilibriumPoints_ZH(prms_ZH);


% ========================================================================
%%% Integrating
% ========================================================================
% -------------------------------------------------
%%% Initial Conditions
% -------------------------------------------------
% r0_n = [rLPs_n(2,1); rLPs_n(2,2); rLPs_n(2,3)];
% v0_n = [0; 0; 0];
% X0_n = [r0_n; v0_n];
% Tf_n = 2*pi;


% X0_n = [1.0018504446716090,0.0030472627733265,0.0025600809355708,0.0031510274334065,-0.0017056019779290,0.0062339010071199]';
X0_n = [1.00185044467,0.00304726277,0.00256008094,0.00315102743,-0.00170560198,0.00623390101]';
% 
% 
% Tf_n = 12.0038099086530963;
%        12.003856214515384
%        12.003808188461973
% X02_n = [1.00185044467,0.00304726277,0.00256008094,0.00315102743,-0.00170560198,0.00623390101]';

Tf_n = 12.003856214515384;
       
% X0_n = [1.00185044467,0.00304726277,0.00256008094,0.00315102743,-0.00170560198,0.00623390101]';

%%% Integration options 
tol = 1e-13;
options = odeset('RelTol',tol,'AbsTol',tol);
options_impact = odeset('Event',@event_Impact_CR3Bn,'RelTol',tol,'AbsTol',tol);

% -------------------------------------------------
%%% Integration
% -------------------------------------------------
% [~, X1_n]    = ode113(@Int_CR3Bn,    [0 Tf_n*2], X02_n, options_impact, prms);

[~, X_n]    = ode113(@Int_CR3Bn,    [0 Tf_n], X0_n, options_impact, prms);
[~, X_ZH_n] = ode113(@Int_CR3Bn_ZH, [0 Tf_n], X0_n, options_impact, prms_ZH);


% ========================================================================
%%% Comparing
% ========================================================================
% -------------------------------------------------
%%% Plotting
% -------------------------------------------------
figure; hold all
p1 = plot3(X_n(:,1),X_n(:,2),X_n(:,3),'r','linewidth',1.5);
p2 = plot3(X_ZH_n(:,1),X_ZH_n(:,2),X_ZH_n(:,3),'b','linewidth',1.5);
PlotBoi3_CR3Bn(22)
plotSecondary(secondary)
axis equal
view(-30, 15)
legend([p1, p2], 'CR3BP', 'CR3BP w ZH')



% ========================================================================
%%% Closeout
% ========================================================================
tocWhole = toc(ticWhole);
fprintf('Elapsed time: %1.4f seconds\n',tocWhole)
















