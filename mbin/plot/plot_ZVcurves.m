% ========================================================================
%%% Description
% ========================================================================
% For plotting zero-velocity curves in the CR3BP

% Created: 1/26/20
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
plot_secondary  = true;
plot_fullSystem = false;

% ========================================================================
%%% Setup
% ========================================================================

% -------------------------------------------------
%%% System
% -------------------------------------------------
% % primary = bodies.earth;     secondary = bodies.moon;
% primary = bodies.mars;      secondary = bodies.phobos;
primary = bodies.jupiter;   secondary = bodies.europa;
% % primary = bodies.jupiter;   secondary = bodies.ganymede;
% % primary = bodies.jupiter;   secondary = bodies.callisto;
% % primary = bodies.saturn;    secondary = bodies.enceladus;
% % primary = bodies.saturn;    secondary = bodies.titan;
% % primary = bodies.neptune;   secondary = bodies.triton;


%%% Normalizing constants
[rNorm, tNorm, vNorm] = cr3bp_norms(primary, secondary, bodies.constants.G);

%%% prms for integration
prms.u  = secondary.MR;
prms.n = 1;
prms.R2 = secondary.R_n;

% prms.u = 5e-3;
% primary.R_n = 0.075;
% secondary.R_n = 0.02;

%%% Equillibrium points
% rLPs_n = collinearEquilibriumPoints_ZH(prms);
rLPs_n = EquilibriumPoints(prms.u, prms.n,1:5);


% -------------------------------------------------
%%% System parameters
% -------------------------------------------------
%%% Select an energy level
JC_L1 = getJacobiConstant_ZH([rLPs_n(1,:), 0, 0, 0], prms);
JC_L2 = getJacobiConstant_ZH([rLPs_n(2,:), 0, 0, 0], prms);
JC_L3 = getJacobiConstant_ZH([rLPs_n(3,:), 0, 0, 0], prms);
JC_L45 = getJacobiConstant_ZH([rLPs_n(4,:), 0, 0, 0], prms);

% JC_des_vec = [2.9, 3.0, JC_L3, (JC_L3+JC_L2)/2, JC_L2, JC_L1, 3.2, 3.3];
% JC_des_vec = fliplr(JC_des_vec);
% JC_des_vec = [JC_L1 + (8.4270e-05), JC_L1 + (8.4270e-06)*2, linspace(JC_L1,JC_L2,5), 3.003600654998641, 3.003600654998641 - 8.4270e-06,  3.003600654998641 - (8.4270e-05)];

% JC_des_vec = [JC_L1+0.02, JC_L1, JC_L2, JC_L3+0.05, JC_L3, JC_L45];
JC_des_vec = L2FlyoverVelocity_2_JC(50, prms.u, rLPs_n(2,:), vNorm, 1);
% JC_des_vec = [JC_L1+0.06];

%%% Select the bounds of the illustrated region
% x_min = rLPs_n(3,1) - 0.5;
% x_max = rLPs_n(2,1) + 0.5;
% y_min = -1.3;
% y_max = 1.3;
% x_min = 0.995;
% x_max = 1.008;
% y_min = -7e-3;
% y_max = 7e-3;

if plot_secondary
    %%% Secondary Body
    x_min = 0.975;
    x_max = 1.025;
    y_min = -0.021;
    y_max = 0.021;
    
    %%% Resolution
    xRes = 750;
    yRes = 750;
end

if plot_fullSystem
    %%% Full system
    x_min = -1.3;
    x_max = 1.3;
    y_min = -1.3;
    y_max = 1.3;

    %%% Resolution
    xRes = 1500;
    yRes = 1500;
end


% ========================================================================
%%% Generate ZV Curves
% ========================================================================
% -------------------------------------------------
%%% Create grid
% -------------------------------------------------
x_vec = linspace(x_min, x_max, xRes);
y_vec = linspace(y_min, y_max, yRes);
[X_xy, Y_xy] = meshgrid(x_vec,y_vec);

% -------------------------------------------------
%%% Calculating JCs across x-y grid
% -------------------------------------------------
JCs_xy = zeros(size(X_xy));
for xk = 1:size(X_xy,1)
    for yk = 1:size(X_xy,2)
        %%% Zero-Velocity Curve        
        JCs_xy(xk,yk) = getJacobiConstant_ZH([X_xy(xk,yk), Y_xy(xk,yk), 0, 0, 0, 0], prms);
    end
end

% -------------------------------------------------
%%% Create contour
% -------------------------------------------------
for kk = 1:length(JC_des_vec)
    figure; hold all
    axis equal
    PlotBoi2('$x_n$','$y_n$',26,'LaTex')

    [xyContourPoints,href] = contour(X_xy,Y_xy,JCs_xy,[JC_des_vec(kk), JC_des_vec(kk)],...
        'color', colors.black, 'linewidth', 2, 'fill', 'on');
    colormap(colors.grn2)
%     [xyContourPoints,href] = contour(X_xy,Y_xy,JCs_xy,[JC_des_vec(kk), JC_des_vec(kk)],...
%         'color', colors.black, 'linewidth', 1.5);
    title(sprintf('JC = %1.12f', JC_des_vec(kk)))
    
    if plot_secondary
        plot(rLPs_n(1:2,1),rLPs_n(1:2,2),'^', 'markeredgecolor',colors.black,'markerfacecolor',colors.ltblue,'markersize',8)
    end
    
    if plot_fullSystem
        plot(rLPs_n(1:5,1),rLPs_n(1:5,2),'^', 'markeredgecolor',colors.black,'markerfacecolor',colors.ltblue,'markersize',8)
    end
    
    %%% plot bodies
%     plotBody2(primary.R/secondary.a, [-prms.u, 0, 0], colors.blue, colors.black, 1, 1)
%     plotBody2(primary.R_n, [-prms.u, 0, 0], colors.blue, colors.black, 1, 1)
%     plotBody2(secondary.R_n, [1-prms.u, 0, 0], colors.blue, colors.black, 1, 1)
    if plot_secondary
        plotBodyTexture3(secondary.R_n, [1-secondary.MR,0,0],secondary.img);
    end
    view(0,90)
%     xlim([0.75 1.25])
%     ylim([-1 1].*0.2)
end
% ========================================================================
%%% Closeout
% ========================================================================
tocWhole = toc(ticWhole);
fprintf('Elapsed time: %1.4f seconds\n',tocWhole)
















