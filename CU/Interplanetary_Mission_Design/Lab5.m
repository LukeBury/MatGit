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
%%% 
% ========================================================================

Earth = bodies.earth;
Sun = bodies.sun;
Mars = bodies.mars;

%%% Initial circular velocity
R0 = Earth.R + 300; % km
vc_Earth = sqrt(Earth.u/R0)

%%% Desired C3
c3_des = 12.8; % km^2/s^2

v0_des = sqrt(c3_des + 2*Earth.u/R0)

dV1 = v0_des - vc_Earth






VelPeriapsis_mars       =   5.6586530502471;
vc_Mars = sqrt(Mars.u/(Mars.R + 200))

VelPeriapsis_mars - vc_Mars



