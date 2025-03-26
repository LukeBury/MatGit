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
mbinPath = '~/Documents/MatGit/mbin';
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




mu = 3.986e5; 

x0 = 6700;
y0 = 0;
vx0 = 0;
vy0 = 7.7;

X0 = [x0; y0; 0; vx0; vy0; 0];


%%% Integration options
tol = 1e-13;
options = odeset('RelTol',tol,'AbsTol',tol);


% [T_out, X_out, t_aps, x_aps, in_aps] = ode113(@Int_CR3BnSTM, [0, PO_IC(end)], [PO_IC(1:6); stm0_colVec], options_apsis, prms);
[T_out, X_out] = ode113(@Int_2BI, [0, 5600], X0, options, mu);

figure(80); hold all
plot3(X_out(:,1),X_out(:,2),X_out(:,3),'linewidth', 2, 'color', colors.mag)
PlotBoi2('x','y',16)
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
















