% ========================================================================
%%% Description
% ========================================================================
% Propagate a trajectory with zonal harmonics and C22s

% Created: 3/27/22
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
%%% PO IC
% -------------------------------------------------

PO_IC = [0.9838431375478787;
 -0.0005086576720824;
 0.0000000000000000;
 0.0006611376703906;
 -0.0145367837458217;
 0.0000000000000000;
 7.8962068915452965];


% ========================================================================
%%% Integrate and plot the POs
% ========================================================================
% -------------------------------------------------
%%% Set up parameters
% -------------------------------------------------
%%% Set primary and secondary bodies
[primary, secondary] = assignPrimaryAndSecondary_CR3BP('Jupiter_Europa.CR3BP_J2pJ4pJ6pJ2s', bodies);

%%% Normalizing constants
[rNorm, tNorm, vNorm] = cr3bp_norms(primary, secondary, bodies.constants.G);
clear tNorm vNorm

%%% prms for integration
prms_ZH.u    = secondary.MR;
prms_ZH.R2 = secondary.R_n;
prms_ZH.J2p = primary.J2;
prms_ZH.J4p = primary.J4;
prms_ZH.J6p = primary.J6;
prms_ZH.J2s = secondary.J2;
prms_ZH.R1 = primary.R / rNorm;
prms_ZH.C22s = secondary.C22;

%%% Determine the mean motion via the ephemeris method
tN_ephemeris = sqrt((secondary.a^3) / (bodies.constants.G*(primary.mass+secondary.mass)));
vN_ephemeris = rNorm / tN_ephemeris;
prms_ZH.n = secondary.meanMot*tN_ephemeris;

%%% Collinear equillibrium points
rLPs_n_ZH = collinearEquilibriumPoints_ZH(prms_ZH);

%%% Set integrator handle
integratorHandle_ZH = @Int_CR3BnSTM_J2pJ4pJ6pJ2s;

%%% Integration options
tol = 1e-13;
options = odeset('RelTol',tol,'AbsTol',tol);
options_apsis = odeset('Events', @event_Apsis_CR3BP, 'RelTol',tol,'AbsTol',tol);

%%% Initialize STM as column vector to be added with state
stm0_colVec = reshape(eye(6),36,1);

[T_out, X_out] = ode113(@Int_CR3BnSTM_J2pJ4pJ6pJ2s_C22s, [0, PO_IC(end)], [PO_IC(1:6); stm0_colVec], options, prms_ZH);
figure; hold all
plot3(X_out(:,1),X_out(:,2),X_out(:,3),'linewidth', 2, 'color', colors.blue2)
PlotBoi3_CR3Bn(28)
plotSecondary(secondary)
axis equal









% ========================================================================
%%% Formatting Structures
% ========================================================================
% -------------------------------------------------
%%% 
% -------------------------------------------------

% --------------------------
% 
% --------------------------

% ========================================================================
%%% Closeout
% ========================================================================
tocWhole = toc(ticWhole);
fprintf('\nElapsed time: %1.4f seconds\n',tocWhole)
















