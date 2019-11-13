% ========================================================================
%%% Description
% ========================================================================
% For integrating trajectories in the normalized circular restricted
% three-body problem with zonal harmonics accounted for

% Created: 11/07/19
% Author : Luke Bury, luke.bury@colorado.edu
% ========================================================================
%%% Initialization
% ========================================================================
clear
clc
close all
ticWhole = tic;

% ========================================================================
%%% Importing Data
% ========================================================================
%%% General data on solar system bodies
currentFolder = pwd;
addpath(genpath(currentFolder))
bodies = getBodyData_exportVersion([currentFolder,'/Functions/SurfaceTextures']);

% ========================================================================
%%% Setup
% ========================================================================
% -------------------------------------------------
%%% Set initial conditions
% -------------------------------------------------
%%% Normalized state [6x1]
X0_n = [1.0029653092197379;
        -0.0000000000000000;
        0.0051782452251263;
        -0.0000000000000040;
        -0.0055255438466834;
        0.0000000000000147];

%%% Normalized final time
Tf_n = 2.4051204546426130;

% -------------------------------------------------
%%% Set up system
% -------------------------------------------------
%%% Choose primary and secondary bodies
primary   = bodies.saturn;
secondary = bodies.enceladus;

%%% Generate normalizing constants (for distance, time, and velocity) ...
%%% ex: "1/rNorm" is representative of 1 km in normalized units
rNorm = secondary.a;         % n <-> km
tNorm = 1/secondary.meanMot; % n <-> sec
vNorm = rNorm / tNorm;       % n <-> km/sec

%%% Set general parameter structure
%%% (prms.u = mass ratio. R1 and R2 are normalized radii
prms.u     = secondary.MR;
prms.R1    = primary.R / rNorm;
prms.R2    = secondary.R_n;

%%% Set desired zonal harmonics terms (any terms not set won't be included)
%%% *note: the 'p' and 's' denote terms for the 'primary' or 'secondary'
%%% body, so here we have J2, J4, and J6 set for the primary, and J2 set
%%% for the secondary
prms.J2p = primary.J2;
prms.J4p = primary.J4;
prms.J6p = primary.J6;
prms.J2s = secondary.J2;

% -------------------------------------------------
%%% Integrator options
% -------------------------------------------------
tol                      = 1e-13;
options                  = odeset('RelTol',tol,'AbsTol',tol);

% ========================================================================
%%% Integrate and Plot
% ========================================================================
% -------------------------------------------------
%%% Integrate trajectory
% -------------------------------------------------
[T2_n, X2_n] = ode113(@Int_CR3Bn_ZH, [0, Tf_n], X0_n, options, prms);

% -------------------------------------------------
%%% Plot trajectory
% -------------------------------------------------
figure; hold all
axis equal
xlabel('$x_n$','FontName','Times New Roman','Fontsize',20,'Interpreter','LaTex')
ylabel('$y_n$','FontName','Times New Roman','Fontsize',20,'Interpreter','LaTex')
zlabel('$z_n$','FontName','Times New Roman','Fontsize',20,'Interpreter','LaTex')
grid on
set(gcf,'color','white')

plot3(X2_n(:,1), X2_n(:,2), X2_n(:,3),'b','linewidth',1.5)

plotSecondary(secondary)

% ========================================================================
%%% Closeout
% ========================================================================
tocWhole = toc(ticWhole);
fprintf('\n(Elapsed time: %1.4f seconds)\n',tocWhole)
















