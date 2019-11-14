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
dataBinpath = '~/CU_Google_Drive/Documents/MatGit/MatlabOutputs/ShallowImpactDataSets';
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
plot_shallowImpactTrajectories = 1;

% ========================================================================
%%% Setup
% ========================================================================
% -------------------------------------------------
%%% Preferences
% -------------------------------------------------
headingLength = 6;

% -------------------------------------------------
%%% Choose data file
% -------------------------------------------------
% dataFile = 'shallowImpacts.M.Saturn_Enceladus.CR3BP.L2_Vertical.nodes1.txt';
% dataFile = 'shallowImpacts.M.Saturn_Enceladus.CR3BP.L2_SHalo.nodes1.txt';

% dataFile = 'shallowImpacts.F.Saturn_Enceladus.CR3BP.L2_SHalo.nodes2000.txt';
dataFile = 'shallowImpacts.F.Saturn_Enceladus.CR3BP.L2_Vertical.nodes2000.txt'; % n_shallowImpacts = 197;



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
c_x0_n              = 1;
c_y0_n              = 2;
c_z0_n              = 3;
c_xd0_n             = 4;
c_yd0_n             = 5;
c_zd0_n             = 6;
c_Tf_n              = 7;
c_lat               = 8;
c_lon               = 9;
c_impactAngle_deg   = 10;
c_impactHeading_lat = 11;
c_impactHeading_lon = 12;
c_JC                = 13;
c_landingVel_mps    = 14;

% -------------------------------------------------
%%% Integrator options
% -------------------------------------------------
tol                      = 1e-13;
options                  = odeset('RelTol',tol,'AbsTol',tol);

% ========================================================================
%%% Make plot of shallow-angle impacts
% ========================================================================
%%% Get number of impacts
n_shallowImpacts = size(shallowAngleImpactData,1);
989
n_shallowImpacts = 197;

%%% Generate figure
figure('position', [156 385 560 420]); hold all
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

for impact_i = 1:n_shallowImpacts
    % --------------------------
    % For plotting the heading directions, since the lat/lon plot isn't
    % equal, we need to artificially increase the magnitude as headings
    % become closer to purely horizontal (since '1' in the longitude-axis
    % is half as long as '1' in the latitude axis). To do this we calculate
    % the angle between the heading and [1,0], and use the cosine of that
    % angle to increase the heading's magnitude
    % --------------------------
    dx = shallowAngleImpactData(impact_i,c_impactHeading_lon);
    dy = shallowAngleImpactData(impact_i,c_impactHeading_lat);
    
    headingVec = [dx, dy]./norm([dx, dy]) .*headingLength;
    angle = acos((dot(headingVec, [1,0])./norm(headingVec)));
    scaling = 1 + abs(cos(angle));
    
    headingVec = headingVec.*scaling;
    
    % --------------------------
    % Plot the heading line
    % --------------------------
    xs = [shallowAngleImpactData(impact_i,c_lon), shallowAngleImpactData(impact_i,c_lon) + headingVec(1)];
    ys = [shallowAngleImpactData(impact_i,c_lat), shallowAngleImpactData(impact_i,c_lat) + headingVec(2)];
    plot(xs, ys, 'color', colors.std.ltred, 'linewidth', 1)
    
    
    % --------------------------
    % Plot the lat/lon impact point (iteratively)
    % --------------------------
    if n_shallowImpacts ~= size(shallowAngleImpactData,1)
        plot(shallowAngleImpactData(impact_i, c_lon), shallowAngleImpactData(impact_i, c_lat), '.','markersize',8,...
        'markeredgecolor',colors.std.red,'markerfacecolor',colors.std.red)
    end
end

% --------------------------
% Plot the lat/lon impact point (all at once)
% --------------------------
if n_shallowImpacts == size(shallowAngleImpactData,1)
    plot(shallowAngleImpactData(:, c_lon), shallowAngleImpactData(:, c_lat), '.','markersize',8,...
        'markeredgecolor',colors.std.red,'markerfacecolor',colors.std.red)
end

% ========================================================================
%%% Plot trajectories
% ========================================================================
if plot_shallowImpactTrajectories
    trajectories = cell(n_shallowImpacts,1);
    
    X0s_n = shallowAngleImpactData(:,c_x0_n:c_zd0_n);
    Tfs_n = shallowAngleImpactData(:,c_Tf_n);
    parfor impact_i = 1:n_shallowImpacts
        X0_n = X0s_n(impact_i,:)';
        Tf_n = Tfs_n(impact_i);
        
        [~, X_n] = ode113(@Int_CR3Bn, [0, Tf_n], X0_n, options, prms);
        
        trajectories{impact_i} = X_n;
    end
    
    figure('position', [717 385 560 420]); hold all
    PlotBoi3_CR3Bn(20)
    axis equal
    plotSecondary(secondary)
    for impact_i = 1:n_shallowImpacts
        plot3(trajectories{impact_i}(:,1), trajectories{impact_i}(:,2), trajectories{impact_i}(:,3), 'r', 'linewidth', 0.5)
    end
end





% ========================================================================
%%% Closeout
% ========================================================================
tocWhole = toc(ticWhole);
fprintf('Elapsed time: %1.4f seconds\n',tocWhole)
















