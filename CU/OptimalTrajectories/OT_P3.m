clear
clc
close all
addpath(genpath('/Users/CU_Google_Drive/lukebury/Documents/MatGit/mbin'))
tic
% ========================================================================
%%% Run/Plot Switches
% ========================================================================

% ========================================================================
%%% Importing Data
% ========================================================================
%%% General data on solar system bodies
bodies = getBodyData();

%%% Color options/schemes
colors = get_colors();

% ========================================================================
%%% Running
% ========================================================================


% % Hohmann cost
Jt_H = @(r) sqrt(2/(r*(1+r)))*(r-1) - (1/sqrt(r))*(sqrt(r)-1);



