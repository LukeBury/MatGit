% ========================================================================
%%% Description
% ========================================================================
% 

% Created: 
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
%%% Set parameters
% -------------------------------------------------
mu = bodies.europa.MR;
G = bodies.constants.G;
a = bodies.europa.a;
rNorm = a;

m1 = bodies.jupiter.mass;
m2 = bodies.europa.mass;

J2p = bodies.jupiter.J2;
J2s = bodies.europa.J2;

A2p = J2p * (bodies.jupiter.R^2);
A2s = J2s * (bodies.europa.R^2);

A2p_n = J2p * ((bodies.jupiter.R/bodies.europa.a)^2);
A2s_n = J2s * (bodies.europa.R_n^2);


% [rNorm, tNorm, vNorm] = cr3bp_norms(primary, secondary, bodies.constants.G);

% ========================================================================
%%% Determining mean motion values
% ========================================================================

% -------------------------------------------------
%%% Classical CR3BP
% -------------------------------------------------
fprintf('================================ Classical ================================\n')

n_n_classical = 1

n_classical = sqrt(G*(m1+m2) / (a^3))

error_classical = (bodies.europa.meanMot - n_classical)*100 / bodies.europa.meanMot

clear prms
prms.u   = mu;
prms.n   = n_n_classical;
rLPs_n_classical   = collinearEquilibriumPoints_ZH(prms);
rLPs_classical = rLPs_n_classical.*rNorm;
clear prms


% -------------------------------------------------
%%% Theory method from paper
% -------------------------------------------------
fprintf('================================ Theory Method ================================\n')

n_theory = sqrt((G/(a^3))*((m1+m2) + 3*(m1*A2p+m2*A2s)/(2*a*a)))

n_n_theory = 1 + 1.5*((1-mu)*A2p_n + mu*A2s_n)

error_theory =  (bodies.europa.meanMot - n_theory)*100 / bodies.europa.meanMot

clear prms
prms.u   = mu;
prms.R1  = bodies.jupiter.R / bodies.europa.a;
prms.R2  = bodies.europa.R_n;
prms.J2p = J2p;
prms.J2s = J2s;
prms.n   = n_n_theory;
rLPs_n_theory   = collinearEquilibriumPoints_ZH(prms);
rLPs_theory = rLPs_n_theory.*rNorm;
deltaL123_theory = rLPs_theory - rLPs_classical
clear prms


% -------------------------------------------------
%%% "Lagrangian Method"
% -------------------------------------------------
fprintf('================================ Lagrangian Method ================================\n')

n_lagrangian = sqrt((G*(m1+m2)/(a^3)) * (1 + 3*(A2p+A2s)/(2*a*a)))

n_n_lagrangian = 1 + 1.5*(A2p_n + A2s_n)

error_lagrangian =  (bodies.europa.meanMot - n_lagrangian)*100 / bodies.europa.meanMot

clear prms
prms.u   = mu;
prms.R1  = bodies.jupiter.R / bodies.europa.a;
prms.R2  = bodies.europa.R_n;
prms.J2p = J2p;
prms.J2s = J2s;
prms.n   = n_n_lagrangian;
rLPs_n_lagrangian   = collinearEquilibriumPoints_ZH(prms);
rLPs_lagrangian = rLPs_n_lagrangian.*rNorm;
deltaL123_lagrangian = rLPs_lagrangian - rLPs_classical
clear prms

% -------------------------------------------------
%%% Ephemeris method from my paper
% -------------------------------------------------
fprintf('================================ Ephemeris Method ================================\n')

tN_ephemeris = sqrt(a*a*a / (G*(m1+m2)));
n_ephemeris = bodies.europa.meanMot
n_n_ephemeris = n_ephemeris*tN_ephemeris
error_ephemeris =  (bodies.europa.meanMot - n_ephemeris)*100 / bodies.europa.meanMot

clear prms
prms.u   = mu;
prms.R1  = bodies.jupiter.R / bodies.europa.a;
prms.R2  = bodies.europa.R_n;
prms.J2p = J2p;
prms.J2s = J2s;
prms.n   = n_n_ephemeris;
rLPs_n_ephemeris   = collinearEquilibriumPoints_ZH(prms);
rLPs_ephemeris = rLPs_n_ephemeris.*rNorm;
deltaL123_ephemeris = rLPs_ephemeris - rLPs_classical
clear prms
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
















