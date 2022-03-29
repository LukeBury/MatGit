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
%%% 13-a
% ========================================================================
% syms e E real
% 
% eHat  = cos(E) - e*(cos(E)^2) - e + e*e*cos(E);
% ePerp = sqrt(1-e^2)*sin(E) - e*sqrt(1-e^2)*cos(E);
% 
% pretty(simplify(int(eHat,E)))
% pretty(simplify(int(ePerp,E)))
% 
% int(eHat,E,[0 2*pi])
% int(ePerp,E,[0 2*pi])

% ========================================================================
%%% 13-b
% ========================================================================
% clear
% syms e E real
% 
% eHat = sin(E);
% ePerp = sqrt(1-e^2)*cos(E);
% 
% pretty(simplify(int(eHat,E)))
% pretty(simplify(int(ePerp,E)))
% 
% int(eHat,E,[0 2*pi])
% int(ePerp,E,[0 2*pi])

% ========================================================================
%%% 13-d
% ========================================================================
clear
syms a e E real

expression = (1-e*cos(E))^2;

pretty(simplify(int(expression,E)))

int(expression,E,[0 2*pi])

int(expression,E,[0 2*pi]) * (a/(2*pi))































