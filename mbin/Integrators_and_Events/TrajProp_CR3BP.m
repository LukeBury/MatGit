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

% PO_IC = [0.9996165286614656;
%  0.0000000000000000;
%  -0.0009890792364341;
%  0.0000000000000000;
%  -0.0161973279788486;
%  0.0000000000000000;
%  4.8871289984803123];

% % PO_IC = [0.9996165392231373;
% %  0.0000000000000000;
% %  -0.0009890876428918;
% %  0.0000000000000000;
% %  -0.0161901083758854;
% %  0.0000000000000000;
% %  4.8873805631699803];

PO_IC = [0.9842567339080418;
 0.0000000000000000;
 0.0000000000000000;
 -0.0000000000000240;
 -0.0162168583881011;
 0.0000000000000000;
 8.0742756181995201];

% ========================================================================
%%% Integrate and plot the POs
% ========================================================================
% -------------------------------------------------
%%% Set up parameters
% -------------------------------------------------
%%% Set primary and secondary bodies
[primary, secondary] = assignPrimaryAndSecondary_CR3BP('Jupiter_Europa.CR3BP', bodies)
% [primary, secondary] = assignPrimaryAndSecondary_CR3BP('Saturn_Enceladus.CR3BP', bodies)

% secondary.MR = 1.898884589251784e-07
% 989
% 989
% 989
% 989

%%% Normalizing constants
[rNorm, tNorm, vNorm] = cr3bp_norms(primary, secondary, bodies.constants.G);

%%% prms for integration
prms.u    = secondary.MR;
prms.R2 = secondary.R_n;
prms.R1 = primary.R / rNorm;

% prms.J2p  = primary.J2;
% prms.J4p  = primary.J4;
% prms.J6p  = primary.J6;
% prms.J2s  = secondary.J2;
% prms.C22s = secondary.C22;
prms.J2p  = primary.J2;
prms.J4p  = primary.J4;
prms.J6p  = primary.J6;
prms.J2s  = secondary.J2;
prms.C22s = secondary.C22;
% prms.C22s = 0;


%%% Determine the mean motion via the ephemeris method
tN_ephemeris = sqrt((secondary.a^3) / (bodies.constants.G*(primary.mass+secondary.mass)));
prms.n = secondary.meanMot*tN_ephemeris;
% prms.n = 1
%%% Set integrator handle
integratorHandle = @Int_CR3BnSTM_J2pJ4pJ6pJ2s_C22s;











%%% Collinear equillibrium points
% rLPs_n = EquilibriumPoints(prms.u, prms.n);
rLPs_n = collinearEquilibriumPoints_ZH_C22s(prms);

%%% Integration options
tol = 1e-13;
options = odeset('RelTol',tol,'AbsTol',tol);
options_apsis = odeset('Events', @event_Apsis_CR3BP, 'RelTol',tol,'AbsTol',tol);

%%% Initialize STM as column vector to be added with state
stm0_colVec = reshape(eye(6),36,1);

% [T_out, X_out, t_aps, x_aps, in_aps] = ode113(@Int_CR3BnSTM, [0, PO_IC(end)], [PO_IC(1:6); stm0_colVec], options_apsis, prms);
[T_out, X_out, t_aps, x_aps, in_aps] = ode113(integratorHandle, [0, PO_IC(end)], [PO_IC(1:6); stm0_colVec], options_apsis, prms);
figure; hold all
plot3(X_out(:,1),X_out(:,2),X_out(:,3),'linewidth', 2, 'color', colors.blue2)
PlotBoi3_CR3Bn(28)
% plotSecondary(secondary)
axis equal

stm_tf_t0                           = reshape(X_out(end,7:42),6,6);
monodromy                           = stm_tf_t0;
[eigenVectors_new, eigenValues_new] = eig(monodromy);
[S1, S2]                            = getStabilityIndices(diag(eigenValues_new));
% 
% 
alts = rowNorm(x_aps(:,1:3) - [1-prms.u, 0, 0]) - prms.R2;
min(alts)*rNorm




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
















