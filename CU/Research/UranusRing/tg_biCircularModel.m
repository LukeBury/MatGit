% ========================================================================
%%% Description
% ========================================================================
% For first tests of running bi-circular model code

% Created: 10/15/19
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
%%% 
% -------------------------------------------------

X0_n = PO_ICs.Jupiter_Europa.CR3BP.L2_Lyapunov;
Tp_n = X0_n(end);
X0_n = X0_n(1:6);

%%% Tolerance
tol = 1e-13;

%%% Options
options = odeset('RelTol',tol,'AbsTol',tol);

jupiter = bodies.jupiter;
europa = bodies.europa;
ganymede = bodies.ganymede;


rNorm = europa.a;            % n <-> km
tNorm = 1 / europa.meanMot;  % n <-> sec
vNorm = rNorm / tNorm; % n <-> km/sec

prms.u = bodies.europa.MR;

% [t, X] = ode113(@Int_CR3Bn, linspace(0, Tp_n, 100), X0_n, options, prms);




global mu mu3 r13 theta0 T2 T3
T2 = europa.Tp / tNorm;
T3 = ganymede.Tp / tNorm;
mu = europa.MR;
r13 = ganymede.a / rNorm;
mu3 = ganymede.mass / (jupiter.mass + europa.mass);
theta0 = pi;




[t_4B, X_4B] = ode113(@xdot_CR4BP, linspace(0, Tp_n, 100), X0_n, options);



figure; hold all
plot3(X_4B(:,1),X_4B(:,2),X_4B(:,3),'r','linewidth',2)
PlotBoi3_CR3Bn(20)
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
















