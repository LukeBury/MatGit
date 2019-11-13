% ========================================================================
%%% Description
% ========================================================================
% This script is designed to plot the results generated in in
% getShallowImpactAngles.m. This script will generate a lat/lon plot of all
% the shallow impact-angle locations achieved by manifolds of a certain
% family of periodic orbits. Results can be color coded by energy

% Created: 11/13/19
% Author : Luke Bury, luke.bury@colorado.edu
% ========================================================================
%%% Initialization
% ========================================================================
clear
clc
close all
mbinPath    = '~/CU_Google_Drive/Documents/MatGit/mbin';
dataBinpath = '~/CU_Google_Drive/Documents/MatGit/MatlabOutputs';
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
%%% Choose data file
% -------------------------------------------------
% dataFile = 'shallowImpacts.M.Saturn_Enceladus.CR3BP.L2_SHalo.nodes1.txt';
dataFile = 'shallowImpacts.F.Saturn_Enceladus.CR3BP.L2_SHalo.nodes2000.txt';




%%% Create path to dataFile
dataFilePath = [dataBinpath, '/',dataFile];

% -------------------------------------------------
%%% Set up system
% -------------------------------------------------
%%% Set primary & secondary
[primary, secondary] = assignPrimaryAndSecondary_CR3BP(dataFile, bodies);

%%% Normalizing constants
rNorm = secondary.a;         % n <-> km
tNorm = 1/secondary.meanMot; % n <-> sec
vNorm = rNorm / tNorm;       % n <-> km/sec

%%% Setting parameter structure
prms.u     = secondary.MR;
prms.rNorm = rNorm;
prms.R1    = primary.R / rNorm;
prms.R2    = secondary.R_n;

if contains(dataFile,'.CR3BP_J2pJ4pJ6pJ2s.')
    prms.J2p = primary.J2; prms.J4p = primary.J4; prms.J6p = primary.J6; prms.J2s = secondary.J2;
end

% -------------------------------------------------
%%% Load PO family data
% -------------------------------------------------
%%% Load shallow angle impact data
shallowAngleImpactData = dlmread(dataFilePath,',',1,0);

%%% Column specifiers
c_x0_n            = 1;
c_y0_n            = 2;
c_z0_n            = 3;
c_xd0_n           = 4;
c_yd0_n           = 5;
c_zd0_n           = 6;
c_lat             = 7;
c_lon             = 8;
c_impactAngle_deg = 9;
c_impactHeading_x = 10;
c_impactHeading_y = 11;
c_JC              = 12;
c_landingVel_mps  = 13;

% -------------------------------------------------
%%% Integrator options
% -------------------------------------------------
tol                      = 1e-13;
options                  = odeset('RelTol',tol,'AbsTol',tol);
% options_impact           = odeset('Event',@event_Impact_CR3Bn,'RelTol',tol,'AbsTol',tol);



% ========================================================================
%%% Make plot of shallow-angle impacts
% ========================================================================
%%% Get number of impacts
n_shallowImpacts = size(shallowAngleImpactData,1);

%%% Generate figure
figure; hold all
xlim([-1, 1] .* 185)
ylim([-1, 1] .* 95)
PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 18, 'LaTex')
% 
% for impact_i = 1:n_shallowImpacts
%     plot(shallowAngleImpactData(impact_i, c_lon), shallowAngleImpactData(impact_i, c_lat), '.','markersize',17,...
%         'markeredgecolor',colors.std.black,'markerfacecolor',colors.std.grey)
%     
%     quiver(shallowAngleImpactData(impact_i, c_lon), shallowAngleImpactData(impact_i, c_lat),...
%         shallowAngleImpactData(impact_i, c_impactHeading_x), shallowAngleImpactData(impact_i, c_impactHeading_y),...
%         10,'color',colors.std.cyan,'linewidth',3)
% 
% end


quiver(shallowAngleImpactData(:, c_lon), shallowAngleImpactData(:, c_lat),...
    shallowAngleImpactData(:, c_impactHeading_x), shallowAngleImpactData(:, c_impactHeading_y),...
    1,'color',colors.std.ltred,'linewidth',1)

% plot(shallowAngleImpactData(:, c_lon), shallowAngleImpactData(:, c_lat), '.','markersize',17,...
%     'markeredgecolor',colors.std.red,'markerfacecolor',colors.std.red)
plot(shallowAngleImpactData(:, c_lon), shallowAngleImpactData(:, c_lat), '.','markersize',8,...
    'markeredgecolor',colors.std.red,'markerfacecolor',colors.std.red)
% ========================================================================
%%% Closeout
% ========================================================================
tocWhole = toc(ticWhole);
fprintf('Elapsed time: %1.4f seconds\n',tocWhole)
















