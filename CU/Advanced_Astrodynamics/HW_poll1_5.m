clear
clc
close all
mbinPath = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/mbin';
addpath(genpath(mbinPath))

% ========================================================================
%%% Importing Data
% ========================================================================
%%% General data on solar system bodies
bodies = getBodyData(mbinPath);

%%% Color options/schemes
colors = get_colors();

% ========================================================================
%%% Run Switches
% ========================================================================


% ========================================================================
%%% Setup
% ========================================================================
% ------------------------------------
%%% 
% ------------------------------------
syms P1 P2 Q1 Q2 real

M = [cos(Q2), -(Q1^-1)*sin(Q2), P2*(Q1^-2)*sin(Q2), -P1*sin(Q2)-P2*(Q1^-1)*cos(Q2);...
     sin(Q2), (Q1^-1)*cos(Q2), -P2*(Q1^-2)*cos(Q2), P1*cos(Q2)-P2*(Q1^-1)*sin(Q2);...
     0, 0, cos(Q2), -Q1*sin(Q2);...
     0, 0, sin(Q2), Q1*cos(Q2)];

J = [zeros(2,2), eye(2,2);...
     -eye(2,2), zeros(2,2)];
 
eqn = M' * J * M;
simplify(eqn)







































