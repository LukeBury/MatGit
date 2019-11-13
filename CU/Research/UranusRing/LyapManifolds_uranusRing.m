% ========================================================================
%%% Description
% ========================================================================
% Looking at unstable manifolds of L2 Lyapunov orbits at Cordelia and L1
% Lyapunov orbits at Ophelia

% Created: 09/13/19
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
%%% Find energy of circular orbits near the ring in each system
% ========================================================================
% -------------------------------------------------
%%% 
% -------------------------------------------------



% ========================================================================
%%% Find the Lyapunov orbit for this energy level in each system
% ========================================================================
% -------------------------------------------------
%%% 
% -------------------------------------------------



% ========================================================================
%%% Propagate the manifolds of the Lyapunov orbit to check ring boundedness
% ========================================================================
% -------------------------------------------------
%%% 
% -------------------------------------------------






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
fprintf('Elapsed time: %1.4f seconds\n',tocWhole)















