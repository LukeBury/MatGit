% ========================================================================
%%% Description
% ========================================================================
% For analyzing the high-latitude-landing-problem data generated in Julia

% Created: 03/20/20
% Author : Luke Bury, luke.bury@colorado.edu
% ========================================================================
%%% Initialization
% ========================================================================
clear
clc
close all
mbinPath = '~/CU_Google_Drive/Documents/MatGit/mbin';
addpath(genpath(mbinPath))
dataPath = '/Volumes/LB_External_Drive/Research/High_Latitude_Landing_Study/Data/Full_Results/';
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
%%% System
% -------------------------------------------------
primary = bodies.jupiter;   secondary = bodies.europa;

% -------------------------------------------------
%%% Normalizing factors and equillibrium points
% -------------------------------------------------
%%% Normalizing factors
% rNorm = secondary.a;         % n <-> km
% tNorm = 1/secondary.meanMot; % n <-> sec
% vNorm = rNorm / tNorm;       % n <-> km/sec
[rNorm, tNorm, vNorm] = cr3bp_norms(primary, secondary, bodies.constants.G);


%%% Collinear equilibrium points
L123 = EquilibriumPoints(secondary.MR,1,1:3); % [3x3] of L1, L2, L3 normalized BCR coordinates

% -------------------------------------------------
%%% Load the data
% -------------------------------------------------
%%% Column specifiers
c_trajID    = 1;
c_x         = 2;
c_y         = 3;
c_z         = 4;
c_xd        = 5;
c_yd        = 6;
c_zd        = 7;
c_t         = 8;
c_JC        = 9;
c_R         = 10;
c_a         = 11;
c_e         = 12;
c_i         = 13;
c_E         = 14;
c_H         = 15;
c_hx        = 16;
c_hy        = 17;
c_hz        = 18;
c_lat       = 19;
c_lon       = 20;
c_isPer     = 21;
c_isApo     = 22;
c_isxyCross = 23;
c_L_del     = 24;
c_G_del     = 25;
c_H_del     = 26;
c_isL1Cross = 27;

% --------------------------
% Load .mat files of full data
% --------------------------
% load([dataPath,'lowTrajs_50mps_neckFull_L1_250km_22v0s_4pi.mat']); % data_L1_4pi
load([dataPath,'lowTrajs_50mps_neckFull_L1_250km_22v0s_1yr.mat']); % data_L1_1yr
% load([dataPath,'lowTrajs_50mps_neckFull_L2_250km_22v0s_4pi.mat']); % data_L2_4pi
% load([dataPath,'lowTrajs_50mps_neckImpact_L2_250km_22v0s_4pi.mat']); % data_L2_4pi_impact
load([dataPath,'lowTrajs_50mps_neckFull_L2_250km_22v0s_1yr.mat']); % data_L2_1yr
load([dataPath,'highTrajs_50mps_70lat_1pi.mat']); % data_high_1pi
data_title = '50 m/s';

% % --------------------------
% % Create near data (inside necks)
% % --------------------------
% %%% Set max distance from secondary ... (max distance to L1 or L2)
% Rmax = max([L123(2,1)-(1-secondary.MR), (1-secondary.MR)-L123(1,1)]);
% 
%%% Get distance to secondary for all neck data
% Rs_L1_1yr = rowNorm(data_L1_1yr(:,c_x:c_z) - [1-secondary.MR,0,0]);
% Rs_L2_1yr = rowNorm(data_L2_1yr(:,c_x:c_z) - [1-secondary.MR,0,0]);
% 
% %%% Find those that are within the neck
% nearIndices_L1   = Rs_L1_1yr <= Rmax;
% nearIndices_L2   = Rs_L2_1yr <= Rmax;
% 
% %%% Create the data
% nearData_L1_1yr = data_L1_1yr(nearIndices_L1,:);
% nearData_L2_1yr = data_L2_1yr(nearIndices_L2,:);
% 
% %%% Save the data
% save([dataPath,'lowTrajs_50mps_neckFull_L1_250km_22v0s_1yr_nearEur','.mat'], 'nearData_L1_1yr','-v7.3')
% save([dataPath,'lowTrajs_50mps_neckFull_L2_250km_22v0s_1yr_nearEur','.mat'], 'nearData_L2_1yr','-v7.3')

% --------------------------
% Load .mat files of near data (inside necks)
% --------------------------
load([dataPath,'lowTrajs_50mps_neckFull_L1_250km_22v0s_1yr_nearEur.mat']); % nearData_L1_1yr
load([dataPath,'lowTrajs_50mps_neckFull_L2_250km_22v0s_1yr_nearEur.mat']); % nearData_L2_1yr

% --------------------------
% Load .mat files of data with x-y crossings
% --------------------------
% Full data, but the only events tracked are x-y-plane crossings
load([dataPath,'lowTrajs_50mps_neckFull_L1_250km_22v0s_1yr_xyCrossings.mat']); % data_L1_1yr_xyCrossings
load([dataPath,'lowTrajs_50mps_neckFull_L2_250km_22v0s_1yr_xyCrossings.mat']); % data_L2_1yr_xyCrossings
load([dataPath,'lowTrajs_50mps_neckFull_L2_125km_22v0s_1yr_xyCrossings.mat']); % data_L2_1yr_xyCrossings_dense
load([dataPath,'highTrajs_50mps_70lat_1pi_xyCrossings.mat']); % data_high_1pi_xyCrossings

% % --------------------------
% % Create near data of trajs with x-y crossings tracked (inside necks)
% % --------------------------
% %%% Set max distance from secondary ... (max distance to L1 or L2)
% Rmax = max([L123(2,1)-(1-secondary.MR), (1-secondary.MR)-L123(1,1)]);
% 
% %%% Calculate distances
% Rs_L1_1yr_xy = rowNorm(data_L1_1yr_xyCrossings(:,c_x:c_z) - [1-secondary.MR,0,0]);
% Rs_L2_1yr_xy = rowNorm(data_L2_1yr_xyCrossings(:,c_x:c_z) - [1-secondary.MR,0,0]);
% Rs_L2_1yr_xy_dense = rowNorm(data_L2_1yr_xyCrossings_dense(:,c_x:c_z) - [1-secondary.MR,0,0]);
% 
% %%% Determine indicies data near Europa
% nearIndices_L1_xyCrossings = find(Rs_L1_1yr_xy <= Rmax);
% nearIndices_L2_xyCrossings = find(Rs_L2_1yr_xy <= Rmax);
% nearIndices_L2_xyCrossings_dense = find(Rs_L2_1yr_xy_dense <= Rmax);
% 
% %%% Create the data
% nearData_L1_1yr_xyCrossings = data_L1_1yr_xyCrossings(nearIndices_L1_xyCrossings,:);
% nearData_L2_1yr_xyCrossings = data_L2_1yr_xyCrossings(nearIndices_L2_xyCrossings,:);
% nearData_L2_1yr_xyCrossings_dense = data_L2_1yr_xyCrossings_dense(nearIndices_L2_xyCrossings_dense,:);
% 
% %%% Save the data
% save([dataPath,'lowTrajs_50mps_neckFull_L1_250km_22v0s_1yr_xyCrossings_nearEur','.mat'], 'nearData_L1_1yr_xyCrossings','-v7.3')
% save([dataPath,'lowTrajs_50mps_neckFull_L2_250km_22v0s_1yr_xyCrossings_nearEur','.mat'], 'nearData_L2_1yr_xyCrossings','-v7.3')
% save([dataPath,'lowTrajs_50mps_neckFull_L2_125km_22v0s_1yr_xyCrossings_nearEur','.mat'], 'nearData_L2_1yr_xyCrossings_dense','-v7.3')

% --------------------------
% Load .mat files of near data (inside necks) with x-y crossings tracked,
% and separate dataset of just those events
% --------------------------
%%% Loading data
load([dataPath,'lowTrajs_50mps_neckFull_L1_250km_22v0s_1yr_xyCrossings_nearEur.mat']); % nearData_L1_1yr_xyCrossings
load([dataPath,'lowTrajs_50mps_neckFull_L2_250km_22v0s_1yr_xyCrossings_nearEur.mat']); % nearData_L2_1yr_xyCrossings
load([dataPath,'lowTrajs_50mps_neckFull_L2_125km_22v0s_1yr_xyCrossings_nearEur.mat']); % nearData_L2_1yr_xyCrossings_dense

%%% Data sets for x-y crossings only
nearData_L1_1yr_xyCrossingsOnly = nearData_L1_1yr_xyCrossings(nearData_L1_1yr_xyCrossings(:,c_isxyCross)==1,:);
nearData_L2_1yr_xyCrossingsOnly = nearData_L2_1yr_xyCrossings(nearData_L2_1yr_xyCrossings(:,c_isxyCross)==1,:);
nearData_L2_1yr_xyCrossingsOnly_dense = nearData_L2_1yr_xyCrossings_dense(nearData_L2_1yr_xyCrossings_dense(:,c_isxyCross)==1,:);
data_high_1pi_xyCrossingsOnly   = data_high_1pi_xyCrossings(data_high_1pi_xyCrossings(:,c_isxyCross)==1,:);

% --------------------------
% Calculate distances to secondary for datasets
% --------------------------
Rs_L1_1yr        = rowNorm(data_L1_1yr(:,c_x:c_z) - [1-secondary.MR,0,0]);
Rs_L2_1yr        = rowNorm(data_L2_1yr(:,c_x:c_z) - [1-secondary.MR,0,0]);
Rs_high_1pi      = rowNorm(data_high_1pi(:,c_x:c_z) - [1-secondary.MR,0,0]);
Rs_L1_1yr_near   = rowNorm(nearData_L1_1yr(:,c_x:c_z) - [1-secondary.MR,0,0]);
Rs_L2_1yr_near   = rowNorm(nearData_L2_1yr(:,c_x:c_z) - [1-secondary.MR,0,0]);

Rs_L1_1yr_xyCrossings   = rowNorm(data_L1_1yr_xyCrossings(:,c_x:c_z) - [1-secondary.MR,0,0]);
Rs_L2_1yr_xyCrossings   = rowNorm(data_L2_1yr_xyCrossings(:,c_x:c_z) - [1-secondary.MR,0,0]);
Rs_L2_1yr_xyCrossings_dense   = rowNorm(data_L2_1yr_xyCrossings_dense(:,c_x:c_z) - [1-secondary.MR,0,0]);
Rs_high_1pi_xyCrossings = rowNorm(data_high_1pi_xyCrossings(:,c_x:c_z) - [1-secondary.MR,0,0]);
Rs_L1_1yr_near_xyCrossings = rowNorm(nearData_L1_1yr_xyCrossings(:,c_x:c_z) - [1-secondary.MR,0,0]);
Rs_L2_1yr_near_xyCrossings = rowNorm(nearData_L2_1yr_xyCrossings(:,c_x:c_z) - [1-secondary.MR,0,0]);
Rs_L2_1yr_near_xyCrossings_dense = rowNorm(nearData_L2_1yr_xyCrossings_dense(:,c_x:c_z) - [1-secondary.MR,0,0]);
Rs_L1_1yr_near_xyCrossingsOnly = rowNorm(nearData_L1_1yr_xyCrossingsOnly(:,c_x:c_z) - [1-secondary.MR,0,0]);
Rs_L2_1yr_near_xyCrossingsOnly = rowNorm(nearData_L2_1yr_xyCrossingsOnly(:,c_x:c_z) - [1-secondary.MR,0,0]);
Rs_L2_1yr_near_xyCrossingsOnly_dense = rowNorm(nearData_L2_1yr_xyCrossingsOnly_dense(:,c_x:c_z) - [1-secondary.MR,0,0]);
Rs_high_1pi_xyCrossingsOnly = rowNorm(data_high_1pi_xyCrossingsOnly(:,c_x:c_z) - [1-secondary.MR,0,0]);


% --------------------------
% labels and colors
% --------------------------
label_L1_4pi        = 'L1, 4pi';
label_L1_1yr        = 'L1, 1 year';
label_L2_4pi        = 'L2, 4pi';
label_L2_4pi_impact = 'L2, 4pi - impact';
label_L2_1yr        = 'L2, 1 year';
label_high_1pi      = 'Lat 70°+, 1pi';

color_L1_4pi        = colors.grn2;
% color_L1_1yr        = colors.purp;
color_L1_1yr        = colors.blue2;
color_L2_4pi        = colors.blue;
color_L2_4pi_impact = colors.brown;
% color_L2_1yr        = colors.cyan;
color_L2_1yr        = colors.red2;
color_high_1pi      = colors.red;


% -------------------------------------------------
%%% Description of data
% -------------------------------------------------
% --------------------------
% Family 1 of data
% --------------------------
% data_L1_1yr   - full data matrix of L1 X0s integrated for 1 year, xy
%                 events are not tracked
% data_L2_1yr   - full data matrix of L2 X0s integrated for 1 year, xy
%                 events are not tracked
% data_high_1pi - full data matrix of high lat (70+) X0s integrated for 1
%                 pi, xy events are not tracked
% 
% nearData_L1_1yr - Subset of data_L1_1yr matrix containing only data from
%                   states inside the neck, xy events are not tracked
% nearData_L2_1yr - Subset of data_L2_1yr matrix containing only data from
%                   states inside the neck, xy events are not tracked
% 
% --------------------------
% Family 2 of data
% --------------------------
% data_L1_1yr_xyCrossings - full data matrix of L1 X0s integrated for 1
%                           year, xy events are tracked
% data_L2_1yr_xyCrossings - full data matrix of L2 X0s integrated for 1
%                           year, xy events are tracked
% data_high_1pi_xyCrossings - full data matrix of high lat (70+) X0s 
%                             integrated for 1, pi, xy events are not 
%                             tracked
% 
% nearData_L1_1yr_xyCrossings - Subset of data_L1_1yr_xyCrossings matrix
%                               containing only data from states inside the
%                               neck, xy events are tracked
% nearData_L2_1yr_xyCrossings - Subset of data_L2_1yr_xyCrossings matrix
%                               containing only data from states inside the
%                               neck, xy events are tracked
% 
% nearData_L1_1yr_xyCrossingsOnly - Subset of nearData_L1_1yr_xyCrossings,
%                                   only contains xy-plane crossing events
% nearData_L2_1yr_xyCrossingsOnly - Subset of nearData_L2_1yr_xyCrossings,
%                                   only contains xy-plane crossing events
% data_high_1pi_xyCrossingsOnly - Subset of data_high_1pi_xyCrossings,
%                                   only contains xy-plane crossing events
%
% nearData_L2_1yr_xyCrossings_dense - Like nearData_L2_1yr_xyCrossings, but
%                                     with a more dense grid search
% nearData_L2_1yr_xyCrossingsOnly_dense - Like nearData_L2_1yr_xyCrossingsOnly,
%                                         but with a more dense grid search
% -------------------------------------------------
%%% Plotting datasets
% -------------------------------------------------
if 1+1==1
%     color_L1_1yr        = colors.purp;
%     color_L2_1yr        = colors.cyan;
%     color_high_1pi      = colors.red;
    
    figure; hold all
    axis equal
    xlim([L123(1,1)-secondary.R_n, L123(2,1)+secondary.R_n])
    PlotBoi3_CR3Bn(26)
    
%     p_high_1pi = plot3(data_high_1pi(:,c_x),data_high_1pi(:,c_y),data_high_1pi(:,c_z),'.','markersize',1, 'color', color_high_1pi);
    p_L1_1yr   = plot3(data_L1_1yr(:,c_x),data_L1_1yr(:,c_y),data_L1_1yr(:,c_z), '.','markersize',1, 'color', color_L1_1yr);
    p_L2_1yr   = plot3(data_L2_1yr(:,c_x),data_L2_1yr(:,c_y),data_L2_1yr(:,c_z), '.','markersize',1, 'color', color_L2_1yr);
    
%     [legh,objh] = legend([p_high_1pi, p_L1_1yr, p_L2_1yr], label_high_1pi  ,label_L1_1yr, label_L2_1yr,'FontSize',14);
    [legh,objh] = legend([ p_L1_1yr, p_L2_1yr]  ,label_L1_1yr, label_L2_1yr,'FontSize',14);
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);
end

% -------------------------------------------------
%%% Plotting lat/lon projections
% -------------------------------------------------
if 1+1==1
    figure; hold all
    PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$',26,'LaTex')
    h = image(xlim, -ylim, bodies.europa.img);
    xlim([-180 180])
    ylim([-90 90])
    
%     p_high_1pi = plot(data_high_1pi(1:200:end,c_lon),data_high_1pi(1:200:end,c_lat),'.','markersize',1, 'color', color_high_1pi);
    p_L1_1yr   = plot(data_L1_1yr(:,c_lon),data_L1_1yr(:,c_lat), '.','markersize',1, 'color', color_L1_1yr);
    p_L2_1yr   = plot(data_L2_1yr(:,c_lon),data_L2_1yr(:,c_lat), '.','markersize',1, 'color', color_L2_1yr);
    
%     [legh,objh] = legend([p_high_1pi, p_L1_1yr, p_L2_1yr], label_high_1pi  ,label_L1_1yr, label_L2_1yr,'FontSize',14);
    [legh,objh] = legend([p_L1_1yr, p_L2_1yr]  ,label_L1_1yr, label_L2_1yr,'FontSize',14);
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);

end

% -------------------------------------------------
%%% Plotting neck comparisons
% -------------------------------------------------clear

if 1+1==1
    figure; hold all
    JC = L2FlyoverVelocity_2_JC(50, secondary.MR, L123(2,:), vNorm, 1);
    prms.u = secondary.MR;
    prms.n = 1;
    prms.R2_n = secondary.R_n;
    plotCR3BP_YZNeck( JC, secondary.MR , 1, 0, prms, color_L1_1yr, 2);
    plotCR3BP_YZNeck( JC, secondary.MR , 2, 0, prms, color_L2_1yr, 2);
    plotBodyTexture3(secondary.R_n,[1-secondary.MR,0,0],secondary.img)
    axis equal
    PlotBoi3_CR3Bn(26)
    legend('L1 Neck','L2 Neck','FontSize',14)
    view(90,0)
end


% -------------------------------------------------
%%% Plotting initial conditions
% -------------------------------------------------
if 1+1==1
    L1_ICs = data_L1_1yr(data_L1_1yr(:,c_t)==0, :);
%     L2_ICs = data_L2_1yr(data_L2_1yr(:,c_t)==0, :);
    L2_ICs = test(:,2:4);
    
    figure; hold all
    JC = L2FlyoverVelocity_2_JC(50, secondary.MR, L123(2,:), vNorm, 1);
    prms.u = secondary.MR;
    prms.n = 1;
    prms.R2_n = secondary.R_n;
    plotCR3BP_YZNeck( JC, secondary.MR , 1, 0, prms, color_L1_1yr, 2);
    plotCR3BP_YZNeck( JC, secondary.MR , 2, 0, prms, color_L2_1yr, 2);
    plotBodyTexture3(secondary.R_n,[1-secondary.MR,0,0],secondary.img)
    
    plot3(L1_ICs(:,c_x), L1_ICs(:,c_y), L1_ICs(:,c_z),'.','color',colors.black)
%     plot3(L2_ICs(:,c_x), L2_ICs(:,c_y), L2_ICs(:,c_z),'.','color',colors.black)
    plot3(L2_ICs(:,1), L2_ICs(:,2), L2_ICs(:,3),'.','color',colors.black)
    
    axis equal
    PlotBoi3_CR3Bn(26)
    legend('L1 Neck','L2 Neck','FontSize',14)
    
    
end




% -------------------------------------------------
%%% Looking at trajectories with maximum inclincation at x-y crossing
% -------------------------------------------------
if 1+1 == 1
%%% Inclination Poincare results
figure; hold all
PlotBoi2('$r_2$','Inclination,$^\circ$',26,'LaTex')
p_high_1pi = plot(Rs_high_1pi_xyCrossingsOnly, data_high_1pi_xyCrossingsOnly(:,c_i),'.','color',color_high_1pi);
p_L1_1yr   = plot(Rs_L1_1yr_near_xyCrossingsOnly, nearData_L1_1yr_xyCrossingsOnly(:,c_i),'.','color',color_L1_1yr);
p_L2_1yr = plot(Rs_L2_1yr_near_xyCrossingsOnly, nearData_L2_1yr_xyCrossingsOnly(:,c_i),'.','color',color_L2_1yr);
p_surf   = plot([secondary.R_n secondary.R_n],[0 140],'k', 'linewidth',1.5);
[legh,objh] = legend([p_high_1pi, p_L1_1yr, p_L2_1yr, p_surf]  , label_high_1pi, label_L1_1yr, label_L2_1yr, 'Europa Surface','FontSize',14);
lineh = findobj(objh,'type','line');
set(lineh,'linestyle','-','linewidth',3);
ylim([0 140])
xlim([0 0.016])

%%% Max inc - L1
figure; hold all
PlotBoi2('$r_2$','Inclination,$^\circ$',26,'LaTex')
p_high_1pi = plot(Rs_high_1pi_xyCrossingsOnly, data_high_1pi_xyCrossingsOnly(:,c_i),'.','color',color_high_1pi);
p_L1_1yr   = plot(Rs_L1_1yr_near_xyCrossingsOnly, nearData_L1_1yr_xyCrossingsOnly(:,c_i),'.','color',color_L1_1yr);
p_surf   = plot([secondary.R_n secondary.R_n],[0 140],'k', 'linewidth',1.5);

% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
figure; hold all
PlotBoi2('$r_2$','Inclination,$^\circ$',26,'LaTex')
p_L1_1yr = plot(Rs_L1_1yr_near_xyCrossingsOnly, nearData_L1_1yr_xyCrossingsOnly(:,c_i),'.','color',colors.grey);
trajID_L1_inc_vs_R = nearData_L1_1yr_xyCrossingsOnly(nearData_L1_1yr_xyCrossingsOnly(:,c_i) == max(nearData_L1_1yr_xyCrossingsOnly(:,c_i)),c_trajID);
L1_maxInc_xycrossingData = nearData_L1_1yr_xyCrossingsOnly(nearData_L1_1yr_xyCrossingsOnly(:,c_trajID)==trajID_L1_inc_vs_R, :);
p_orbit = plot(rowNorm(L1_maxInc_xycrossingData(:,c_x:c_z)-[1-secondary.MR,0,0]), L1_maxInc_xycrossingData(:,c_i),'^','markersize',7,'markeredgecolor',colors.black,'markerfacecolor',colors.cyan);
p_surf   = plot([secondary.R_n secondary.R_n],[0 140],'k', 'linewidth',1.5);
ylim([0 140])
xlim([0 0.016])
[legh,objh] = legend([p_L1_1yr, p_orbit, p_surf],label_L1_1yr, 'Single Orbit', 'Europa Surface','FontSize',14);
% [legh,objh] = legend([p_L1_1yr, p_surf],label_L1_1yr, 'Europa Surface','FontSize',14);
lineh = findobj(objh,'type','line');
lineh = lineh([1,2,5,6]);
% lineh = lineh([1,2,3,4]);
set(lineh,'linestyle','-','linewidth',3);
% 3017 is the traj ID of that peak around inclination 37.5, r2 2.6e-3
% testTraj = nearData_L1_1yr_xyCrossings(nearData_L1_1yr_xyCrossings(:,c_trajID) == 3017, :);
% figure; hold all
% plot3(testTraj(:,c_x), testTraj(:,c_y), testTraj(:,c_z),'linewidth', 1, 'color', colors.blue2)
% PlotBoi3_CR3Bn(26)
% figure; hold all
% indices = 14000:14580;
% plot3(testTraj(indices,c_x), testTraj(indices,c_y), testTraj(indices,c_z),'linewidth', 2, 'color', colors.blue2)
% PlotBoi3_CR3Bn(26)

% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

figure; hold all
PlotBoi3_CR3Bn(26)
L1_maxInc_traj = nearData_L1_1yr_xyCrossings(nearData_L1_1yr_xyCrossings(:,c_trajID)==trajID_L1_inc_vs_R,:);
L1_maxInc_traj_full = data_L1_1yr_xyCrossings(data_L1_1yr_xyCrossings(:,c_trajID)==trajID_L1_inc_vs_R,:);
% plot3(L1_maxInc_traj(:,c_x),L1_maxInc_traj(:,c_y),L1_maxInc_traj(:,c_z),'color',color_L1_1yr)
plot3(L1_maxInc_traj(10:end,c_x),L1_maxInc_traj(10:end,c_y),L1_maxInc_traj(10:end,c_z),'color',colors.blue2)
% plot3(L1_maxInc_traj_full(:,c_x),L1_maxInc_traj_full(:,c_y),L1_maxInc_traj_full(:,c_z),'color',color_L1_1yr)
plot3(L123(1:2,1),[0,0],[0,0],'^','markeredgecolor',colors.black, 'markerfacecolor', colors.blue)

%%% Max inc - L2
figure; hold all
PlotBoi2('$r_2$','Inclination,$^\circ$',26,'LaTex')
p_L2_1yr = plot(Rs_L2_1yr_near_xyCrossingsOnly, nearData_L2_1yr_xyCrossingsOnly(:,c_i),'.','color',color_L2_1yr);
p_high_1pi = plot(Rs_high_1pi_xyCrossingsOnly, data_high_1pi_xyCrossingsOnly(:,c_i),'.','color',color_high_1pi);
p_surf   = plot([secondary.R_n secondary.R_n],[0 140],'k', 'linewidth',1.5);
    

L2_maxInc_trajID = nearData_L2_1yr_xyCrossingsOnly(nearData_L2_1yr_xyCrossingsOnly(:,c_i) == max(nearData_L2_1yr_xyCrossingsOnly(:,c_i)),c_trajID);
L2_maxInc_xycrossingData = nearData_L2_1yr_xyCrossingsOnly(nearData_L2_1yr_xyCrossingsOnly(:,c_trajID)==L2_maxInc_trajID, :);
plot(rowNorm(L2_maxInc_xycrossingData(:,c_x:c_z)-[1-secondary.MR,0,0]), L2_maxInc_xycrossingData(:,c_i),'.','color',colors.black);

[legh,objh] = legend([p_high_1pi, p_L1_1yr, p_L2_1yr, p_surf], label_high_1pi  ,label_L1_1yr, label_L2_1yr, 'Europa Surface','FontSize',14);
lineh = findobj(objh,'type','line');
set(lineh,'linestyle','-','linewidth',3);
    
figure; hold all
PlotBoi3_CR3Bn(26)
L2_maxInc_traj = nearData_L2_1yr_xyCrossings(nearData_L2_1yr_xyCrossings(:,c_trajID)==L2_maxInc_trajID,:);
% L2_maxInc_traj_full = data_L2_1yr_xyCrossings(data_L2_1yr_xyCrossings(:,c_trajID)==L2_maxInc_trajID,:);
plot3(L2_maxInc_traj(:,c_x),L2_maxInc_traj(:,c_y),L2_maxInc_traj(:,c_z),'color',colors.blue2)
% plot3(L2_maxInc_traj_full(:,c_x),L2_maxInc_traj_full(:,c_y),L2_maxInc_traj_full(:,c_z),'color',color_L2_1yr)
plot3(L123(1:2,1),[0,0],[0,0],'^','markeredgecolor',colors.black, 'markerfacecolor', colors.blue)


%%% Max inc - L2 dense
figure; hold all
PlotBoi2('$r_2$','Inclination,$^\circ$',26,'LaTex')
plot(Rs_high_1pi_xyCrossingsOnly, data_high_1pi_xyCrossingsOnly(:,c_i),'.','color',color_high_1pi);
plot(Rs_L2_1yr_near_xyCrossingsOnly_dense, nearData_L2_1yr_xyCrossingsOnly_dense(:,c_i),'.','color',color_L2_1yr);

L2_maxInc_trajID_dense = nearData_L2_1yr_xyCrossingsOnly_dense(nearData_L2_1yr_xyCrossingsOnly_dense(:,c_i) == max(nearData_L2_1yr_xyCrossingsOnly_dense(:,c_i)),c_trajID);
L2_maxInc_data_dense = nearData_L2_1yr_xyCrossingsOnly_dense(nearData_L2_1yr_xyCrossingsOnly_dense(:,c_trajID)==L2_maxInc_trajID_dense, :);
plot(rowNorm(L2_maxInc_data_dense(:,c_x:c_z)-[1-secondary.MR,0,0]), L2_maxInc_data_dense(:,c_i),'.','color',colors.black);

figure; hold all
PlotBoi3_CR3Bn(26)
L2_maxInc_traj_dense = nearData_L2_1yr_xyCrossings_dense(nearData_L2_1yr_xyCrossings_dense(:,c_trajID)==L2_maxInc_trajID_dense,:);
L2_maxInc_traj_dense_full = data_L2_1yr_xyCrossings_dense(data_L2_1yr_xyCrossings_dense(:,c_trajID)==L2_maxInc_trajID_dense,:);
% plot3(L2_maxInc_traj_dense(:,c_x),L2_maxInc_traj_dense(:,c_y),L2_maxInc_traj_dense(:,c_z),'color',colors.blue)
plot3(L2_maxInc_traj_dense_full(:,c_x),L2_maxInc_traj_dense_full(:,c_y),L2_maxInc_traj_dense_full(:,c_z),'color',colors.blue)
plot3(L123(1:2,1),[0,0],[0,0],'^','markeredgecolor',colors.black, 'markerfacecolor', colors.blue)

%3790:5000
figure; hold all; axis equal
% range = 3790:5000;
range = 3000:6000;
plot3(L2_maxInc_traj_dense(range,c_x),L2_maxInc_traj_dense(range,c_y),L2_maxInc_traj_dense(range,c_z))
axis normal


%%% trying to get a better maxInc orbit from L2-dense
newData = nearData_L2_1yr_xyCrossingsOnly_dense(nearData_L2_1yr_xyCrossingsOnly_dense(:,c_i)>25.87,:);
newData = newData(newData(:,c_i)<25.88,:);
newData = newData(newData(:,c_R)<0.00892, :);
newData = newData(newData(:,c_R)>0.00888, :);
newTrajID = newData(1,c_trajID);
newTraj =  nearData_L2_1yr_xyCrossings_dense(nearData_L2_1yr_xyCrossings_dense(:,c_trajID)==newTrajID,:);
figure; hold all
PlotBoi3_CR3Bn(26)
plot3(newTraj(:,c_x),newTraj(:,c_y),newTraj(:,c_z),'color',colors.blue)
newTraj_full = data_L2_1yr_xyCrossings_dense(data_L2_1yr_xyCrossings_dense(:,c_trajID)==newTrajID,:);
figure; hold all
PlotBoi3_CR3Bn(26)
plot3(newTraj_full(:,c_x),newTraj_full(:,c_y),newTraj_full(:,c_z),'color',colors.blue2)
plot3(L123(1:2,1),[0,0],[0,0],'^','markeredgecolor',colors.black, 'markerfacecolor', colors.blue2)

% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
figure; hold all
PlotBoi2('$r_2$','Inclination,$^\circ$',26,'LaTex')
p_L2_1yr = plot(Rs_L2_1yr_near_xyCrossingsOnly_dense, nearData_L2_1yr_xyCrossingsOnly_dense(:,c_i),'.','color',colors.grey);
newTraj_maxInc_data_dense = nearData_L2_1yr_xyCrossingsOnly_dense(nearData_L2_1yr_xyCrossingsOnly_dense(:,c_trajID)==newTrajID, :);
p_orbit = plot(rowNorm(newTraj_maxInc_data_dense(:,c_x:c_z)-[1-secondary.MR,0,0]), newTraj_maxInc_data_dense(:,c_i),'^','markersize',7,'markeredgecolor',colors.black,'markerfacecolor',colors.cyan);
p_surf   = plot([secondary.R_n secondary.R_n],[0 140],'k', 'linewidth',1.5);
ylim([0 140])
xlim([0 0.016])
[legh,objh] = legend([p_L2_1yr, p_orbit, p_surf],label_L2_1yr, 'Single Orbit', 'Europa Surface','FontSize',14);
% [legh,objh] = legend([p_L2_1yr, p_surf],label_L2_1yr, 'Europa Surface','FontSize',14);
lineh = findobj(objh,'type','line');
lineh = lineh([1,2,5,6]);
lineh = lineh([1,2,3,4]);
set(lineh,'linestyle','-','linewidth',3);
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

% --------------------------
%%% Trying to get this PO
% --------------------------
% % x0 = 1.018;
% x0 = L123(2,1)-0.002;
% % yd0 = -0.0005;
% yd0 = -0.0042411;
% X0 = [x0; 0; 0; 0; yd0; 0];
X0 = [1.018461802379925;
 -0.000000277654002;
 0.000000000000000;
 -0.000001832694740;
 -0.004241055848417;
 0.000000000000000;
 7.395200000000000];

prms.u  = secondary.MR;
prms.n  = 1;
tol     = 1e-13;
options = odeset('RelTol',tol,'AbsTol',tol);
[T, X] = ode113(@Int_CR3Bn, linspace(0, 7.3952,10000), X0(1:6), options, prms);
figure; hold all
plot3(X(:,1),X(:,2),X(:,3),'linewidth',2,'color',colors.blue)
axis normal
PlotBoi3_CR3Bn(26)

stm0        = eye(6);
stm0_colVec = reshape(stm0,36,1);
[~, X_new] = ode113(@Int_CR3BnSTM, [0, X0(end)], [X0(1:6); stm0_colVec], options, prms);
%%% Stability indices of new PO
stm_tf_t0                           = reshape(X_new(end,7:42),6,6);
monodromy                           = stm_tf_t0;
[eigenVectors_new, eigenValues_new] = eig(monodromy);
[S1, S2]                            = getStabilityIndices(diag(eigenValues_new));

end


% -------------------------------------------------
%%% Looking at trajectories with maximum zdot at x-y crossing
% -------------------------------------------------
if 1+1==1
%%% Poincaré section results
figure; hold all
p_high_1pi = plot(Rs_high_1pi_xyCrossingsOnly, abs(data_high_1pi_xyCrossingsOnly(:,c_zd)),'.','color',color_high_1pi);
p_L1_1yr = plot(Rs_L1_1yr_near_xyCrossingsOnly, abs(nearData_L1_1yr_xyCrossingsOnly(:,c_zd)),'.','color',color_L1_1yr);
p_L2_1yr = plot(Rs_L2_1yr_near_xyCrossingsOnly, nearData_L2_1yr_xyCrossingsOnly(:,c_zd),'.','color',color_L2_1yr);
p_surf   = plot([secondary.R_n secondary.R_n],[0 3],'k', 'linewidth',1.5);
[legh,objh] = legend([p_high_1pi, p_L1_1yr, p_L2_1yr, p_surf]  , label_high_1pi, label_L1_1yr, label_L2_1yr, 'Europa Surface','FontSize',14);
lineh = findobj(objh,'type','line');
set(lineh,'linestyle','-','linewidth',3);
ylim([0 0.09])
xlim([0 0.016])
PlotBoi2('$r_2$','$|\dot{z}|$',26,'LaTex')


%%% Max zdot - L1
figure; hold all
p_high_1pi = plot(Rs_high_1pi_xyCrossingsOnly, abs(data_high_1pi_xyCrossingsOnly(:,c_zd)),'.','color',color_high_1pi);
p_L1_1yr = plot(Rs_L1_1yr_near_xyCrossingsOnly, abs(nearData_L1_1yr_xyCrossingsOnly(:,c_zd)),'.','color',color_L1_1yr);


L1_maxZd_trajID = nearData_L1_1yr_xyCrossingsOnly(nearData_L1_1yr_xyCrossingsOnly(:,c_zd) == max(nearData_L1_1yr_xyCrossingsOnly(:,c_zd)),c_trajID);
L1_maxZd_xycrossingData = nearData_L1_1yr_xyCrossingsOnly(nearData_L1_1yr_xyCrossingsOnly(:,c_trajID)==L1_maxZd_trajID, :);
plot(rowNorm(L1_maxZd_xycrossingData(:,c_x:c_z)-[1-secondary.MR,0,0]), abs(L1_maxZd_xycrossingData(:,c_zd)),'.','color',colors.black);

figure; hold all
L1_maxZd_traj = nearData_L1_1yr_xyCrossings(nearData_L1_1yr_xyCrossings(:,c_trajID)==L1_maxZd_trajID,:);
plot3(L1_maxZd_traj(:,c_x),L1_maxZd_traj(:,c_y),L1_maxZd_traj(:,c_z),'color',color_L1_1yr)
plot3(L123(1:2,1),[0,0],[0,0],'^','markeredgecolor',colors.black, 'markerfacecolor', colors.blue2)


%%% Max zdot - L2
figure; hold all
p_high_1pi = plot(Rs_high_1pi_xyCrossingsOnly, data_high_1pi_xyCrossingsOnly(:,c_zd),'.','color',color_high_1pi);
p_L2_1yr = plot(Rs_L2_1yr_near_xyCrossingsOnly, nearData_L2_1yr_xyCrossingsOnly(:,c_zd),'.','color',color_L2_1yr);

L2_maxZd_trajID = nearData_L2_1yr_xyCrossingsOnly(nearData_L2_1yr_xyCrossingsOnly(:,c_zd) == max(nearData_L2_1yr_xyCrossingsOnly(:,c_zd)),c_trajID);
L2_maxZd_xycrossingData = nearData_L2_1yr_xyCrossingsOnly(nearData_L2_1yr_xyCrossingsOnly(:,c_trajID)==L2_maxZd_trajID, :);
plot(rowNorm(L2_maxZd_xycrossingData(:,c_x:c_z)-[1-secondary.MR,0,0]), L2_maxZd_xycrossingData(:,c_zd),'.','color',colors.black);

figure; hold all
L2_maxZd_traj = nearData_L2_1yr_xyCrossings(nearData_L2_1yr_xyCrossings(:,c_trajID)==L2_maxZd_trajID,:);
plot3(L2_maxZd_traj(:,c_x),L2_maxZd_traj(:,c_y),L2_maxZd_traj(:,c_z),'color',color_L2_1yr)


%%% Max zdot - L2 dense
figure; hold all
plot(Rs_L2_1yr_near_xyCrossingsOnly_dense, nearData_L2_1yr_xyCrossingsOnly_dense(:,c_zd),'.','color',color_L2_1yr);
plot(Rs_high_1pi_xyCrossingsOnly, data_high_1pi_xyCrossingsOnly(:,c_zd),'.','color',color_high_1pi);

L2_maxZd_trajID_dense = nearData_L2_1yr_xyCrossingsOnly_dense(nearData_L2_1yr_xyCrossingsOnly_dense(:,c_zd) == max(nearData_L2_1yr_xyCrossingsOnly_dense(:,c_zd)),c_trajID);
L2_maxZd_data_dense = nearData_L2_1yr_xyCrossingsOnly_dense(nearData_L2_1yr_xyCrossingsOnly_dense(:,c_trajID)==L2_maxZd_trajID_dense, :);
plot(rowNorm(L2_maxZd_data_dense(:,c_x:c_z)-[1-secondary.MR,0,0]), L2_maxZd_data_dense(:,c_zd),'.','color',colors.black);

figure; hold all
L2_maxZd_traj_dense = nearData_L2_1yr_xyCrossings_dense(nearData_L2_1yr_xyCrossings_dense(:,c_trajID)==L2_maxZd_trajID_dense,:);
plot3(L2_maxZd_traj_dense(:,c_x),L2_maxZd_traj_dense(:,c_y),L2_maxZd_traj_dense(:,c_z),'color',color_L2_1yr)



figure; hold all
PlotBoi2('$r_2$','$|\dot{z}|$',26,'LaTex')
p_high_1pi = plot(Rs_high_1pi_xyCrossingsOnly, abs(data_high_1pi_xyCrossingsOnly(:,c_zd)),'.','color',color_high_1pi);
p_L1_1yr   = plot(Rs_L1_1yr_near_xyCrossingsOnly, abs(nearData_L1_1yr_xyCrossingsOnly(:,c_zd)),'.','color',color_L1_1yr);
p_L2_1yr   = plot(Rs_L2_1yr_near_xyCrossingsOnly, abs(nearData_L2_1yr_xyCrossingsOnly(:,c_zd)),'.','color',color_L2_1yr);
p_surf   = plot([secondary.R_n secondary.R_n],[0 3],'k', 'linewidth',1.5);
[legh,objh] = legend([p_high_1pi, p_L1_1yr, p_L2_1yr, p_surf], label_high_1pi  ,label_L1_1yr, label_L2_1yr, 'Europa Surface','FontSize',14);
lineh = findobj(objh,'type','line');
set(lineh,'linestyle','-','linewidth',3);
xlim([0 0.016])
PlotBoi2('$r_2$','$|\dot{z}|$',26,'LaTex')

% ==================================================================
%%% trying to get a better maxZd orbit from L1 data
% ==================================================================
newData = nearData_L1_1yr_xyCrossingsOnly(nearData_L1_1yr_xyCrossingsOnly(:,c_zd)>0.0208,:);
newData = newData(newData(:,c_zd)<0.0212,:);
newData = newData(newData(:,c_R)<0.00965, :);
newData = newData(newData(:,c_R)>0.0096, :);
newTrajID = newData(1,c_trajID);
newTraj =  nearData_L1_1yr_xyCrossings(nearData_L1_1yr_xyCrossings(:,c_trajID)==newTrajID,:);
figure; hold all
PlotBoi3_CR3Bn(26)
plot3(newTraj(:,c_x),newTraj(:,c_y),newTraj(:,c_z),'color',colors.blue)
newTraj_full = data_L1_1yr_xyCrossings(data_L1_1yr_xyCrossings(:,c_trajID)==newTrajID,:);
figure; hold all
PlotBoi3_CR3Bn(26)
plot3(newTraj_full(:,c_x),newTraj_full(:,c_y),newTraj_full(:,c_z),'color',colors.blue)


%%% This looks like the underlying PO
figure; hold all
PlotBoi3_CR3Bn(26)
indices = 14000:14580; %Tp from 14000-14580 is 9.966629008513991
plot3(newTraj_full(indices,c_x),newTraj_full(indices,c_y),newTraj_full(indices,c_z),'color',colors.blue)

% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

figure; hold all
PlotBoi2('$r_2$','$|\dot{z}|$',26,'LaTex')
newData_xyCrossingsOnly = nearData_L1_1yr_xyCrossingsOnly(nearData_L1_1yr_xyCrossingsOnly(:,c_trajID)==newTrajID,:);
% p_high_1pi = plot(Rs_high_1pi_xyCrossingsOnly, abs(data_high_1pi_xyCrossingsOnly(:,c_zd)),'.','color',color_high_1pi);
p_L1_1yr = plot(Rs_L1_1yr_near_xyCrossingsOnly, abs(nearData_L1_1yr_xyCrossingsOnly(:,c_zd)),'.','color',colors.grey);
p_orbit = plot(rowNorm(newData_xyCrossingsOnly(:,c_x:c_z)-[1-secondary.MR,0,0]), abs(newData_xyCrossingsOnly(:,c_zd)),'^','markersize',7,'markeredgecolor',colors.black,'markerfacecolor',colors.cyan);
p_surf   = plot([secondary.R_n secondary.R_n],[0 0.2],'k', 'linewidth',1.5);
[legh,objh] = legend([p_L1_1yr, p_orbit, p_surf]  ,label_L1_1yr, 'Single Orbit', 'Europa Surface','FontSize',14);
% [legh,objh] = legend([p_L1_1yr, p_surf]  ,label_L1_1yr, 'Europa Surface','FontSize',14);
lineh = findobj(objh,'type','line');
lineh = lineh([1,2,5,6]);
% lineh = lineh([1,2,3,4]);
set(lineh,'linestyle','-','linewidth',3);
ylim([0 0.09])
xlim([0 0.016])
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$





% ==================================================================
%%% trying to get a better maxZd orbit from L2 dense data
% ==================================================================
newData = nearData_L2_1yr_xyCrossingsOnly_dense(nearData_L2_1yr_xyCrossingsOnly_dense(:,c_zd)>0.012,:);
newData = newData(newData(:,c_zd)<0.013,:);
newData = newData(newData(:,c_R)<0.0098, :);
newData = newData(newData(:,c_R)>0.00975, :);
newTrajID = newData(1,c_trajID);
newTraj =  nearData_L2_1yr_xyCrossings_dense(nearData_L2_1yr_xyCrossings_dense(:,c_trajID)==newTrajID,:);
figure; hold all
PlotBoi3_CR3Bn(26)
plot3(newTraj(:,c_x),newTraj(:,c_y),newTraj(:,c_z),'color',colors.blue)
newTraj_full = data_L2_1yr_xyCrossings_dense(data_L2_1yr_xyCrossings_dense(:,c_trajID)==newTrajID,:);
figure; hold all
PlotBoi3_CR3Bn(26)
plot3(newTraj_full(:,c_x),newTraj_full(:,c_y),newTraj_full(:,c_z),'color',colors.blue2)
plot3(L123(1:2,1),[0,0],[0,0],'^','markeredgecolor',colors.black, 'markerfacecolor', colors.blue2)

% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
figure; hold all
PlotBoi2('$r_2$','$|\dot{z}|$',26,'LaTex')
newData_xyCrossingsOnly = nearData_L2_1yr_xyCrossingsOnly_dense(nearData_L2_1yr_xyCrossingsOnly_dense(:,c_trajID)==newTrajID,:);
% p_high_1pi = plot(Rs_high_1pi_xyCrossingsOnly, abs(data_high_1pi_xyCrossingsOnly(:,c_zd)),'.','color',color_high_1pi);
p_L2_1yr = plot(Rs_L2_1yr_near_xyCrossingsOnly_dense, abs(nearData_L2_1yr_xyCrossingsOnly_dense(:,c_zd)),'.','color',colors.grey);
p_orbit = plot(rowNorm(newData_xyCrossingsOnly(:,c_x:c_z)-[1-secondary.MR,0,0]), abs(newData_xyCrossingsOnly(:,c_zd)),'^','markersize',7,'markeredgecolor',colors.black,'markerfacecolor',colors.cyan);
p_surf   = plot([secondary.R_n secondary.R_n],[0 0.2],'k', 'linewidth',1.5);
[legh,objh] = legend([p_L2_1yr, p_orbit, p_surf]  ,label_L2_1yr, 'Single Orbit', 'Europa Surface','FontSize',14);
% [legh,objh] = legend([p_L2_1yr, p_surf]  ,label_L2_1yr, 'Europa Surface','FontSize',14);
lineh = findobj(objh,'type','line');
lineh = lineh([1,2,5,6]);
% lineh = lineh([1,2,3,4]);
set(lineh,'linestyle','-','linewidth',3);
ylim([0 0.09])
xlim([0 0.016])
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
% $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

end


% ==================================================================
%%% Comparing z and zd within the neck
% ==================================================================
if 1+1==1
    figure; hold all
    PlotBoi2('$r_2$','$z$',26,'LaTex')
    plot(nearData_L2_1yr_xyCrossings_dense(:,c_R),nearData_L2_1yr_xyCrossings_dense(:,c_z), '.', 'markersize', 10,'color',colors.grey)
    xlim([0 L123(2,1)-1+prms.u])
    
    figure; hold all
    PlotBoi2('$r_2$','$\dot{z}$',26,'LaTex')
    plot(nearData_L2_1yr_xyCrossings_dense(:,c_R),nearData_L2_1yr_xyCrossings_dense(:,c_zd), '.', 'markersize', 10,'color',colors.grey)
    xlim([0 L123(2,1)-1+prms.u])
    
end


% ========================================================================
%%% Testing
% ========================================================================
% warning('Need to address the fact that I haven''t integrated the high trajs for very long')
% 
% 
% data = data_high_1pi_xyCrossingsOnly(Rs_high_1pi_xyCrossingsOnly>5.5e-3,:);
% data = data((rowNorm(data(:,c_x:c_z) - [1-prms.u,0,0]))<5.6e-3,:);
% data = data(data(:,c_i)<36,:);
% data = data(data(:,c_i)>35.5,:);
% 
% high_minInc_trajID = 48846;
% high_minInc_data = data_high_1pi_xyCrossingsOnly(data_high_1pi_xyCrossingsOnly(:,c_trajID)==high_minInc_trajID,:);
% plot(rowNorm(high_minInc_data(:,c_x:c_z)-[1-secondary.MR,0,0]), high_minInc_data(:,c_i),'.','color',colors.black);
% 
% 
% high_minInc_traj = data_high_1pi_xyCrossings(data_high_1pi_xyCrossings(:,c_trajID)==high_minInc_trajID,:);
% 
% 
% 
% figure; hold all
% plot3(high_minInc_traj(:,c_x),high_minInc_traj(:,c_y),high_minInc_traj(:,c_z))
% 
% 
% options_event = odeset('Event',@event_zEqualsZero_nonTerminal, 'RelTol',tol,'AbsTol',tol);
% X0 = high_minInc_traj(1,c_x:c_zd)';
% [T, X, time_event, X_event, index_event] = ode113(@Int_CR3Bn, [0 20*pi], X0, options_event, prms);
% figure; hold all
% plot3(X(:,1),X(:,2),X(:,3),'r')
% axis normal
% PlotBoi3_CR3Bn(26)
% 
% 
%     
% [X_SCI_event] = X_BaCR2SCI(X_event, time_event, prms);
% a_e_i_raan_w_ta = zeros(size(X_SCI_event,1),6);
% 
% for kk = 1:size(X_SCI_event,1)
%     [a,e,i,raan,w,ta] = ECI2OE(X_SCI_event(kk,1:3), X_SCI_event(kk,4:6),prms.u);
% 
%     a_e_i_raan_w_ta(kk,:) = [a,e,i,raan,w,ta];
% end
% 
% plot(rowNorm(X_event(:,1:3) - [1-prms.u,0,0]), a_e_i_raan_w_ta(:,3)*180/pi,'m.')

if 1+1==1
    tol     = 1e-13;
    options_event = odeset('Event',@event_zEqualsZero_nonTerminal, 'RelTol',tol,'AbsTol',tol);
    X0 = L1_maxInc_traj(10,c_x:c_zd)';
    [T, X, time_event, X_event, index_event] = ode113(@Int_CR3Bn, [0 646], X0, options_event, prms);
    figure; hold all
    plot3(X(:,1),X(:,2),X(:,3),'r')
    axis normal
    PlotBoi3_CR3Bn(26)

    figure; hold all
    plot3(X_event(:,1),X_event(:,2),X_event(:,3),'k.','markersize',10)
    PlotBoi3_CR3Bn(26)
    axis normal
end


% Just looking at a trajectory from L2 dense data
if 1+1==1
    newData = nearData_L2_1yr_xyCrossings_dense(nearData_L2_1yr_xyCrossings_dense(:,c_z)>0.00358,:);
    newData = newData(newData(:,c_z)<0.0036,:);
    newData = newData(newData(:,c_R)<0.0164, :);
    newData = newData(newData(:,c_R)>0.0163, :);
    newTrajID = newData(1,c_trajID);
    newTraj =  nearData_L2_1yr_xyCrossings_dense(nearData_L2_1yr_xyCrossings_dense(:,c_trajID)==newTrajID,:);
    
    figure; hold all
    PlotBoi3_CR3Bn(26)
    axis equal
    plot3(newTraj(:,c_x),newTraj(:,c_y),newTraj(:,c_z),'b')
    plot3(newData(1,c_x),newData(1,c_y),newData(1,c_z),'o','markeredgecolor',colors.black,'markerfacecolor',colors.red)
end









% Looking at inclination in rotating frame
if 1+1==1
%     Rs_high_1pi_xyCrossingsOnly
%     data_high_1pi_xyCrossingsOnly
%     nearData_L1_1yr_xyCrossingsOnly

    is_SCR_high_xyCrossingsOnly = NaN(size(Rs_high_1pi_xyCrossingsOnly));
    is_SCR_L1_xyCrossingsOnly = NaN(size(Rs_L1_1yr_near_xyCrossingsOnly));
    is_SCR_L2_xyCrossingsOnly = NaN(size(Rs_L2_1yr_near_xyCrossingsOnly));
    
    for kk = 1:length(is_SCR_high_xyCrossingsOnly)
        [a,e,i,raan,w,ta] = ECI2OE(data_high_1pi_xyCrossingsOnly(kk,c_x:c_z)-[1-secondary.MR,0,0],data_high_1pi_xyCrossingsOnly(kk,c_xd:c_zd),secondary.MR);
        is_SCR_high_xyCrossingsOnly(kk) = i*180/pi;
    end
    
    for kk = 1:length(is_SCR_L1_xyCrossingsOnly)
        [a,e,i,raan,w,ta] = ECI2OE(nearData_L1_1yr_xyCrossingsOnly(kk,c_x:c_z)-[1-secondary.MR,0,0],nearData_L1_1yr_xyCrossingsOnly(kk,c_xd:c_zd),secondary.MR);
        is_SCR_L1_xyCrossingsOnly(kk) = i*180/pi;
    end
    
    for kk = 1:length(is_SCR_L2_xyCrossingsOnly)
        [a,e,i,raan,w,ta] = ECI2OE(nearData_L2_1yr_xyCrossingsOnly(kk,c_x:c_z)-[1-secondary.MR,0,0],nearData_L2_1yr_xyCrossingsOnly(kk,c_xd:c_zd),secondary.MR);
        is_SCR_L2_xyCrossingsOnly(kk) = i*180/pi;
    end
    
    figure; hold all
    PlotBoi2('$r_2$','Inclination (SCR), $^\circ$',26,'LaTex')
    p_high_1pi = plot(Rs_high_1pi_xyCrossingsOnly, is_SCR_high_xyCrossingsOnly,'.','color',color_high_1pi);
    p_L1_1yr   = plot(Rs_L1_1yr_near_xyCrossingsOnly, is_SCR_L1_xyCrossingsOnly,'.','color',color_L1_1yr);
    p_L2_1yr   = plot(Rs_L2_1yr_near_xyCrossingsOnly, is_SCR_L2_xyCrossingsOnly,'.','color',color_L2_1yr);
    p_surf   = plot([secondary.R_n secondary.R_n],[0 180],'k', 'linewidth',1.5);
    [legh,objh] = legend([p_high_1pi, p_L1_1yr, p_L2_1yr, p_surf], label_high_1pi  ,label_L1_1yr, label_L2_1yr, 'Europa Surface','FontSize',14);
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);
    xlim([0 0.016])
    
end





% -------------------------------------------------
%%% Looking at Latitude vs R (1/19/22)
% -------------------------------------------------
if 1+1 == 1
    % -------------------------------------------------
    %%% L1 Latitude vs R
    % -------------------------------------------------
    figure; hold all
    PlotBoi2('$r_2$','Latitude,$^\circ$',26,'LaTex')
    p_L1_1yr   = plot3(nearData_L1_1yr(:,c_R), nearData_L1_1yr(:,c_lat), nearData_L1_1yr(:,c_trajID),'.','color',colors.mdgrey);
    p_surf   = plot([secondary.R_n secondary.R_n],[-90 90],'k', 'linewidth',1.5);
%     testTraj_L1_1 = nearData_L1_1yr(nearData_L1_1yr(:,c_trajID) == 2596,:); 
    testTraj_L1_1 = data_L1_1yr(data_L1_1yr(:,c_trajID) == 682,:); 
    p_testTraj_L1_1 = plot3(testTraj_L1_1(:,c_R), testTraj_L1_1(:,c_lat), ones(size(testTraj_L1_1(:,c_R))).*4000, '.', 'markersize', 7, 'color', colors.blue2);
    [legh,objh] = legend([p_L1_1yr, p_testTraj_L1_1, p_surf], label_L1_1yr, 'Single Orbit', 'Europa Surface','FontSize',14);
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);
%     xlim([0 0.016])
    xlim([0 0.02])

    figure; hold all
    PlotBoi3_CR3Bn(26)
    plot3(testTraj_L1_1(:,c_x), testTraj_L1_1(:,c_y), testTraj_L1_1(:,c_z),'color',color_L1_1yr)
    plot3(L123(1:2,1),[0,0],[0,0],'^','markeredgecolor',colors.black, 'markerfacecolor', colors.blue)
    xlim([0.977 1.023])



figure; hold all
PlotBoi2('$r_2$','Latitude,$^\circ$',26,'LaTex')
p_L1_1yr   = plot3(nearData_L1_1yr(:,c_R), nearData_L1_1yr(:,c_lat), nearData_L1_1yr(:,c_trajID),'.','color',colors.mdgrey);
p_surf   = plot([secondary.R_n secondary.R_n],[-90 90],'k', 'linewidth',1.5);
testTraj_L1_1 = data_L1_1yr(data_L1_1yr(:,c_trajID) == 682,:); 
p_testTraj_L1_1 = plot3(testTraj_L1_1(:,c_R), testTraj_L1_1(:,c_lat), ones(size(testTraj_L1_1(:,c_R))).*4000, '.', 'markersize', 7, 'color', colors.blue2);
[legh,objh] = legend([p_L1_1yr, p_testTraj_L1_1, p_surf], label_L1_1yr, 'Single Orbit', 'Europa Surface','FontSize',14);
lineh = findobj(objh,'type','line');
set(lineh,'linestyle','-','linewidth',3);
xlim([0 0.02]) 
% index = 13000:14000;
% index = 16000:17000;
% index = 30000:31000;
index = 1:size(testTraj_L1_1,1);
figure; hold all
PlotBoi3_CR3Bn(26)
plot3(testTraj_L1_1(index,c_x), testTraj_L1_1(index,c_y), testTraj_L1_1(index,c_z),'color',color_L1_1yr)
plot3(L123(1:2,1),[0,0],[0,0],'^','markeredgecolor',colors.black, 'markerfacecolor', colors.blue)
xlim([0.977 1.023])





    figure; hold all
    PlotBoi2('$r_2$','Latitude,$^\circ$',26,'LaTex')
    p_L1_1yr   = plot3(nearData_L1_1yr(:,c_R), nearData_L1_1yr(:,c_lat), nearData_L1_1yr(:,c_trajID),'.','color',colors.mdgrey);
    p_surf   = plot([secondary.R_n secondary.R_n],[-90 90],'k', 'linewidth',1.5);
%     testTraj_L1_2 = nearData_L1_1yr(nearData_L1_1yr(:,c_trajID) == 2441,:); 
    testTraj_L1_2 = data_L1_1yr(data_L1_1yr(:,c_trajID) == 2441,:); 
    p_testTraj_L1_2 = plot3(testTraj_L1_2(:,c_R), testTraj_L1_2(:,c_lat), ones(size(testTraj_L1_2(:,c_R))).*4000, '.', 'markersize', 7, 'color', colors.blue2);
    [legh,objh] = legend([p_L1_1yr, p_testTraj_L1_2, p_surf], label_L1_1yr, 'Single Orbit', 'Europa Surface','FontSize',14);
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);
%     xlim([0 0.016])
    xlim([0 0.02])

    figure; hold all
    PlotBoi3_CR3Bn(26)
    plot3(testTraj_L1_2(:,c_x), testTraj_L1_2(:,c_y), testTraj_L1_2(:,c_z),'color',color_L1_1yr)
    plot3(L123(1:2,1),[0,0],[0,0],'^','markeredgecolor',colors.black, 'markerfacecolor', colors.blue)
    xlim([0.977 1.023])
    
    % -------------------------------------------------
    %%% L2 Latitude vs R
    % -------------------------------------------------
    figure; hold all
    PlotBoi2('$r_2$','Latitude,$^\circ$',26,'LaTex')
    p_L2_1yr = plot3(nearData_L2_1yr_xyCrossings_dense(:,c_R), nearData_L2_1yr_xyCrossings_dense(:,c_lat), nearData_L2_1yr_xyCrossings_dense(:,c_trajID),'.','color',colors.mdgrey);
    p_surf   = plot([secondary.R_n secondary.R_n],[-90 90],'k', 'linewidth',1.5);
%     testTraj_L2_1 = nearData_L2_1yr_xyCrossings_dense(nearData_L2_1yr_xyCrossings_dense(:,c_trajID) == 2266,:); 
    testTraj_L2_1 = data_L2_1yr_xyCrossings_dense(data_L2_1yr_xyCrossings_dense(:,c_trajID) == 2266,:); 
    p_testTraj_L2_1 = plot3(testTraj_L2_1(:,c_R), testTraj_L2_1(:,c_lat), ones(size(testTraj_L2_1(:,c_R))).*5000, '.', 'markersize', 7, 'color', colors.blue2);
    [legh,objh] = legend([p_L2_1yr, p_testTraj_L2_1, p_surf], label_L2_1yr, 'Single Orbit', 'Europa Surface','FontSize',14);
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);
%     xlim([0 0.016])
    xlim([0 0.02])

    figure; hold all
    PlotBoi3_CR3Bn(26)
    plot3(testTraj_L2_1(:,c_x), testTraj_L2_1(:,c_y), testTraj_L2_1(:,c_z),'color',color_L1_1yr)
    plot3(L123(1:2,1),[0,0],[0,0],'^','markeredgecolor',colors.black, 'markerfacecolor', colors.blue)
    xlim([0.977 1.023])



    figure; hold all
    PlotBoi2('$r_2$','Latitude,$^\circ$',26,'LaTex')
    p_L2_1yr = plot3(nearData_L2_1yr_xyCrossings_dense(:,c_R), nearData_L2_1yr_xyCrossings_dense(:,c_lat), nearData_L2_1yr_xyCrossings_dense(:,c_trajID),'.','color',colors.mdgrey);
    p_surf   = plot([secondary.R_n secondary.R_n],[-90 90],'k', 'linewidth',1.5);
%     testTraj_L2_2 = nearData_L2_1yr_xyCrossings_dense(nearData_L2_1yr_xyCrossings_dense(:,c_trajID) == 3652,:); 
    testTraj_L2_2 = data_L2_1yr_xyCrossings_dense(data_L2_1yr_xyCrossings_dense(:,c_trajID) == 3652,:); 
    p_testTraj_L2_2 = plot3(testTraj_L2_2(:,c_R), testTraj_L2_2(:,c_lat), ones(size(testTraj_L2_2(:,c_R))).*5000, '.', 'markersize', 7, 'color', colors.blue2);
    [legh,objh] = legend([p_L2_1yr, p_testTraj_L2_2, p_surf], label_L2_1yr, 'Single Orbit', 'Europa Surface','FontSize',14);
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);
%     xlim([0 0.016])
    xlim([0 0.02])

    figure; hold all
    PlotBoi3_CR3Bn(26)
    plot3(testTraj_L2_2(:,c_x), testTraj_L2_2(:,c_y), testTraj_L2_2(:,c_z),'color',color_L1_1yr)
    plot3(L123(1:2,1),[0,0],[0,0],'^','markeredgecolor',colors.black, 'markerfacecolor', colors.blue)
    xlim([0.977 1.023])


end


% ========================================================================
%%% Figures in paper
% ========================================================================
if 1+1==1
    
    % -------------------------------------------------
    %%% Plotting neck comparisons
    % -------------------------------------------------
    figure; hold all
    JC = L2FlyoverVelocity_2_JC(50, secondary.MR, L123(2,:), vNorm, 1);
    prms.u = secondary.MR;
    prms.n = 1;
    prms.R2_n = secondary.R_n;
    plotCR3BP_YZNeck( JC, secondary.MR , 1, 0, prms, color_L1_1yr, 2);
    plotCR3BP_YZNeck( JC, secondary.MR , 2, 0, prms, color_L2_1yr, 2);
    plotBodyTexture3(secondary.R_n,[1-secondary.MR,0,0],secondary.img)
    axis equal
    ax = gca;
    ax.FontSize;
    ax.FontSize = 18;
    PlotBoi3_CR3Bn(26)
    legend('L1 Neck','L2 Neck','FontSize',16)
    view(90,0)
    xlim([0.978 1.022])

    % -------------------------------------------------
    %%% Plotting all the states
    % -------------------------------------------------
    figure; hold all
    axis equal
    xlim([L123(1,1)-secondary.R_n, L123(2,1)+secondary.R_n])
    ax = gca;
    ax.FontSize;
    ax.FontSize = 14;
    PlotBoi3_CR3Bn(28)
    PlotBoi2('$x_n$', '$y_n$', 26, 'LaTex')
    
%     p_high_1pi = plot3(data_high_1pi(:,c_x),data_high_1pi(:,c_y),data_high_1pi(:,c_z),'.','markersize',1, 'color', color_high_1pi);
%     p_L1_1yr   = plot3(data_L1_1yr(:,c_x),data_L1_1yr(:,c_y),data_L1_1yr(:,c_z), '.','markersize',1, 'color', color_L1_1yr);
%     p_L2_1yr   = plot3(data_L2_1yr(:,c_x),data_L2_1yr(:,c_y),data_L2_1yr(:,c_z), '.','markersize',1, 'color', color_L2_1yr);
    p_L1_1yr   = plot(data_L1_1yr(:,c_x),data_L1_1yr(:,c_z), '.','markersize',1, 'color', color_L1_1yr);
    p_L2_1yr   = plot(data_L2_1yr(:,c_x),data_L2_1yr(:,c_z), '.','markersize',1, 'color', color_L2_1yr);
    
%     [legh,objh] = legend([p_high_1pi, p_L1_1yr, p_L2_1yr], label_high_1pi  ,label_L1_1yr, label_L2_1yr,'FontSize',14);
    [legh,objh] = legend([ p_L1_1yr, p_L2_1yr]  ,label_L1_1yr, label_L2_1yr,'FontSize',16);
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);


    % -------------------------------------------------
    %%% Latitude/Longitdue projection of all states
    % -------------------------------------------------
    figure; hold all
    ax = gca;
    ax.FontSize;
    ax.FontSize = 14;
    PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$',26,'LaTex')
    h = image(xlim, -ylim, bodies.europa.img);
    xlim([-180 180])
    ylim([-90 90])
    
%     p_high_1pi = plot(data_high_1pi(1:200:end,c_lon),data_high_1pi(1:200:end,c_lat),'.','markersize',1, 'color', color_high_1pi);
    p_L1_1yr   = plot(data_L1_1yr(:,c_lon),data_L1_1yr(:,c_lat), '.','markersize',1, 'color', color_L1_1yr);
    p_L2_1yr   = plot(data_L2_1yr(:,c_lon),data_L2_1yr(:,c_lat), '.','markersize',1, 'color', color_L2_1yr);
    
%     [legh,objh] = legend([p_high_1pi, p_L1_1yr, p_L2_1yr], label_high_1pi  ,label_L1_1yr, label_L2_1yr,'FontSize',14);
    [legh,objh] = legend([p_L1_1yr, p_L2_1yr]  ,label_L1_1yr, label_L2_1yr,'FontSize',16);
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);

    % -------------------------------------------------
    %%% Poincare Results: all Inc vs R
    % -------------------------------------------------
    figure; hold all
    ax = gca;
    ax.FontSize;
    ax.FontSize = 14;
    PlotBoi2('$r_2$','Inclination,$^\circ$',26,'LaTex')
    p_high_1pi = plot(Rs_high_1pi_xyCrossingsOnly, data_high_1pi_xyCrossingsOnly(:,c_i),'.','color',colors.ltbrown);
    p_L1_1yr   = plot(Rs_L1_1yr_near_xyCrossingsOnly, nearData_L1_1yr_xyCrossingsOnly(:,c_i),'.','color',color_L1_1yr);
    p_L2_1yr = plot(Rs_L2_1yr_near_xyCrossingsOnly, nearData_L2_1yr_xyCrossingsOnly(:,c_i),'.','color',color_L2_1yr);
    p_surf   = plot([secondary.R_n secondary.R_n],[0 140],'k', 'linewidth',1.5);
    [legh,objh] = legend([p_high_1pi, p_L1_1yr, p_L2_1yr, p_surf]  , label_high_1pi, label_L1_1yr, label_L2_1yr, 'Europa Surface','FontSize',16);
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);
    ylim([0 140])
%     xlim([0 0.016])
    xlim([0, 0.0205])


    % -------------------------------------------------
    %%% Data picked out from L1 Inc vs R results
    % -------------------------------------------------
    figure; hold all
    ax = gca;
    ax.FontSize;
    ax.FontSize = 14;
    PlotBoi2('$r_2$','Inclination,$^\circ$',26,'LaTex')
    p_L1_1yr = plot(Rs_L1_1yr_near_xyCrossingsOnly, nearData_L1_1yr_xyCrossingsOnly(:,c_i),'.','color',colors.grey);
    trajID_L1_inc_vs_R = nearData_L1_1yr_xyCrossingsOnly(nearData_L1_1yr_xyCrossingsOnly(:,c_i) == max(nearData_L1_1yr_xyCrossingsOnly(:,c_i)),c_trajID);
%     trajID_L1_inc_vs_R = 3017 % This is the spike
    L1_maxInc_xycrossingData = nearData_L1_1yr_xyCrossingsOnly(nearData_L1_1yr_xyCrossingsOnly(:,c_trajID)==trajID_L1_inc_vs_R, :);
    p_orbit = plot(rowNorm(L1_maxInc_xycrossingData(:,c_x:c_z)-[1-secondary.MR,0,0]), L1_maxInc_xycrossingData(:,c_i),'^','markersize',7,'markeredgecolor',colors.black,'markerfacecolor',colors.cyan);
    p_surf   = plot([secondary.R_n secondary.R_n],[0 140],'k', 'linewidth',1.5);
%     ylim([0 140])
%     xlim([0 0.016])
    xlim([0, 0.0205])
    ylim([0 60])
    [legh,objh] = legend([p_L1_1yr, p_orbit, p_surf],label_L1_1yr, 'Single Orbit', 'Europa Surface','FontSize',16);
    % [legh,objh] = legend([p_L1_1yr, p_surf],label_L1_1yr, 'Europa Surface','FontSize',14);
    lineh = findobj(objh,'type','line');
    lineh = lineh([1,2,5,6]);
    % lineh = lineh([1,2,3,4]);
    set(lineh,'linestyle','-','linewidth',3);
    
    % 3017 is the traj ID of that peak around inclination 37.5, r2 2.6e-3
    % testTraj = nearData_L1_1yr_xyCrossings(nearData_L1_1yr_xyCrossings(:,c_trajID) == 3017, :);
    % figure; hold all
    % plot3(testTraj(:,c_x), testTraj(:,c_y), testTraj(:,c_z),'linewidth', 1, 'color', colors.blue2)
    % PlotBoi3_CR3Bn(26)
    % figure; hold all
    % indices = 14000:14580;
    % plot3(testTraj(indices,c_x), testTraj(indices,c_y), testTraj(indices,c_z),'linewidth', 2, 'color', colors.blue2)
    % PlotBoi3_CR3Bn(26)
    
    
    % -------------------------------------------------
    %%% Trajectory picked from L1 inc vs R results
    % -------------------------------------------------
    figure; hold all
    ax = gca;
    ax.FontSize;
    ax.FontSize = 14;
    PlotBoi3_CR3Bn(28)
    L1_maxInc_traj = nearData_L1_1yr_xyCrossings(nearData_L1_1yr_xyCrossings(:,c_trajID)==trajID_L1_inc_vs_R,:);
    L1_maxInc_traj_full = data_L1_1yr_xyCrossings(data_L1_1yr_xyCrossings(:,c_trajID)==trajID_L1_inc_vs_R,:);
    % plot3(L1_maxInc_traj(:,c_x),L1_maxInc_traj(:,c_y),L1_maxInc_traj(:,c_z),'color',color_L1_1yr)
    plot3(L1_maxInc_traj(10:end,c_x),L1_maxInc_traj(10:end,c_y),L1_maxInc_traj(10:end,c_z),'color',colors.blue2)
    % plot3(L1_maxInc_traj_full(:,c_x),L1_maxInc_traj_full(:,c_y),L1_maxInc_traj_full(:,c_z),'color',color_L1_1yr)
    plot3(L123(1:2,1),[0,0],[0,0],'^','markeredgecolor',colors.black, 'markerfacecolor', colors.blue)
    xlim([0.976, 1.024])

    % -------------------------------------------------
    %%% Data picked out from L2 Inc vs R results
    % -------------------------------------------------
    figure; hold all
    ax = gca;
    ax.FontSize;
    ax.FontSize = 14;
    PlotBoi2('$r_2$','Inclination,$^\circ$',26,'LaTex')
    newData = nearData_L2_1yr_xyCrossingsOnly_dense(nearData_L2_1yr_xyCrossingsOnly_dense(:,c_i)>25.87,:);
    newData = newData(newData(:,c_i)<25.88,:);
    newData = newData(newData(:,c_R)<0.00892, :);
    newData = newData(newData(:,c_R)>0.00888, :);
    trajID_L2_inc_vs_R = newData(1,c_trajID);
    p_L2_1yr = plot(Rs_L2_1yr_near_xyCrossingsOnly_dense, nearData_L2_1yr_xyCrossingsOnly_dense(:,c_i),'.','color',colors.grey);
    newTraj_maxInc_data_dense = nearData_L2_1yr_xyCrossingsOnly_dense(nearData_L2_1yr_xyCrossingsOnly_dense(:,c_trajID)==trajID_L2_inc_vs_R, :);
    p_orbit = plot(rowNorm(newTraj_maxInc_data_dense(:,c_x:c_z)-[1-secondary.MR,0,0]), newTraj_maxInc_data_dense(:,c_i),'^','markersize',7,'markeredgecolor',colors.black,'markerfacecolor',colors.cyan);
    p_surf   = plot([secondary.R_n secondary.R_n],[0 140],'k', 'linewidth',1.5);
%     ylim([0 140])
%     xlim([0 0.016])
    xlim([0, 0.0205])
    ylim([0 60])
    [legh,objh] = legend([p_L2_1yr, p_orbit, p_surf],label_L2_1yr, 'Single Orbit', 'Europa Surface','FontSize',16);
    % [legh,objh] = legend([p_L2_1yr, p_surf],label_L2_1yr, 'Europa Surface','FontSize',14);
    lineh = findobj(objh,'type','line');
    lineh = lineh([1,2,5,6]);
    lineh = lineh([1,2,3,4]);
    set(lineh,'linestyle','-','linewidth',3);
    
    % -------------------------------------------------
    %%% Trajectory picked from L2 inc vs R results
    % -------------------------------------------------
%     newTraj =  nearData_L2_1yr_xyCrossings_dense(nearData_L2_1yr_xyCrossings_dense(:,c_trajID)==trajID_L2_inc_vs_R,:);
%     figure; hold all
%     PlotBoi3_CR3Bn(26)
%     plot3(newTraj(:,c_x),newTraj(:,c_y),newTraj(:,c_z),'color',colors.blue)

    newTraj_full = data_L2_1yr_xyCrossings_dense(data_L2_1yr_xyCrossings_dense(:,c_trajID)==trajID_L2_inc_vs_R,:);
    figure; hold all
    ax = gca;
    ax.FontSize;
    ax.FontSize = 14;
    PlotBoi3_CR3Bn(28)
    plot3(newTraj_full(:,c_x),newTraj_full(:,c_y),newTraj_full(:,c_z),'color',colors.blue2)
    plot3(L123(1:2,1),[0,0],[0,0],'^','markeredgecolor',colors.black, 'markerfacecolor', colors.blue)
    xlim([0.976, 1.024])
    ylim([-1 1].*0.015)
    % -------------------------------------------------
    %%% Poincare Results: all zdot vs R
    % -------------------------------------------------
    figure; hold all
    p_high_1pi = plot(Rs_high_1pi_xyCrossingsOnly, abs(data_high_1pi_xyCrossingsOnly(:,c_zd)),'.','color',colors.ltbrown);
    p_L1_1yr = plot(Rs_L1_1yr_near_xyCrossingsOnly, abs(nearData_L1_1yr_xyCrossingsOnly(:,c_zd)),'.','color',color_L1_1yr);
    p_L2_1yr = plot(Rs_L2_1yr_near_xyCrossingsOnly, nearData_L2_1yr_xyCrossingsOnly(:,c_zd),'.','color',color_L2_1yr);
    p_surf   = plot([secondary.R_n secondary.R_n],[0 3],'k', 'linewidth',1.5);
    [legh,objh] = legend([p_high_1pi, p_L1_1yr, p_L2_1yr, p_surf]  , label_high_1pi, label_L1_1yr, label_L2_1yr, 'Europa Surface','FontSize',16);
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);
    ylim([0 0.09])
%     xlim([0 0.016])
    xlim([0, 0.0205])
    ax = gca;
    ax.FontSize;
    ax.FontSize = 14;
    PlotBoi2('$r_2$','$|\dot{z}|$',26,'LaTex')
    
    
    % -------------------------------------------------
    %%% Data picked out from L1 zdot vs R results
    % -------------------------------------------------
    figure; hold all
    ax = gca;
    ax.FontSize;
    ax.FontSize = 14;
    PlotBoi2('$r_2$','$|\dot{z}|$',26,'LaTex')
    newData = nearData_L1_1yr_xyCrossingsOnly(nearData_L1_1yr_xyCrossingsOnly(:,c_zd)>0.0208,:);
    newData = newData(newData(:,c_zd)<0.0212,:);
    newData = newData(newData(:,c_R)<0.00965, :);
    newData = newData(newData(:,c_R)>0.0096, :);
    newTrajID = newData(1,c_trajID);
    trajID_L1_zdot_vs_R = newData(1,c_trajID);
    newData_xyCrossingsOnly = nearData_L1_1yr_xyCrossingsOnly(nearData_L1_1yr_xyCrossingsOnly(:,c_trajID)==trajID_L1_zdot_vs_R,:);
    % p_high_1pi = plot(Rs_high_1pi_xyCrossingsOnly, abs(data_high_1pi_xyCrossingsOnly(:,c_zd)),'.','color',color_high_1pi);
    p_L1_1yr = plot(Rs_L1_1yr_near_xyCrossingsOnly, abs(nearData_L1_1yr_xyCrossingsOnly(:,c_zd)),'.','color',colors.grey);
    p_orbit = plot(rowNorm(newData_xyCrossingsOnly(:,c_x:c_z)-[1-secondary.MR,0,0]), abs(newData_xyCrossingsOnly(:,c_zd)),'^','markersize',7,'markeredgecolor',colors.black,'markerfacecolor',colors.cyan);
    p_surf   = plot([secondary.R_n secondary.R_n],[0 0.2],'k', 'linewidth',1.5);
    [legh,objh] = legend([p_L1_1yr, p_orbit, p_surf]  ,label_L1_1yr, 'Single Orbit', 'Europa Surface','FontSize',16);
    % [legh,objh] = legend([p_L1_1yr, p_surf]  ,label_L1_1yr, 'Europa Surface','FontSize',14);
    lineh = findobj(objh,'type','line');
    lineh = lineh([1,2,5,6]);
    % lineh = lineh([1,2,3,4]);
    set(lineh,'linestyle','-','linewidth',3);
    ylim([0 0.09])
%     xlim([0 0.016])
    xlim([0, 0.0205])
    
    % -------------------------------------------------
    %%% Trajectory picked from L1 zdot vs R results
    % -------------------------------------------------
    newTraj =  nearData_L1_1yr_xyCrossings(nearData_L1_1yr_xyCrossings(:,c_trajID)==trajID_L1_zdot_vs_R,:);
    figure; hold all
    ax = gca;
    ax.FontSize;
    ax.FontSize = 14;
    PlotBoi3_CR3Bn(28)
    plot3(newTraj(:,c_x),newTraj(:,c_y),newTraj(:,c_z),'color',colors.blue2)
    plot3(L123(1:2,1),[0,0],[0,0],'^','markeredgecolor',colors.black, 'markerfacecolor', colors.blue)
    xlim([0.976, 1.024])
    % -------------------------------------------------
    %%% Data picked out from L2 zdot vs R results
    % -------------------------------------------------
    newData = nearData_L2_1yr_xyCrossingsOnly_dense(nearData_L2_1yr_xyCrossingsOnly_dense(:,c_zd)>0.012,:);
    newData = newData(newData(:,c_zd)<0.013,:);
    newData = newData(newData(:,c_R)<0.0098, :);
    newData = newData(newData(:,c_R)>0.00975, :);
    trajID_L2_zdot_vs_R = newData(1,c_trajID);
    figure; hold all
    ax = gca;
    ax.FontSize;
    ax.FontSize = 14;
    PlotBoi2('$r_2$','$|\dot{z}|$',26,'LaTex')
    newData_xyCrossingsOnly = nearData_L2_1yr_xyCrossingsOnly_dense(nearData_L2_1yr_xyCrossingsOnly_dense(:,c_trajID)==trajID_L2_zdot_vs_R,:);
    % p_high_1pi = plot(Rs_high_1pi_xyCrossingsOnly, abs(data_high_1pi_xyCrossingsOnly(:,c_zd)),'.','color',color_high_1pi);
    p_L2_1yr = plot(Rs_L2_1yr_near_xyCrossingsOnly_dense, abs(nearData_L2_1yr_xyCrossingsOnly_dense(:,c_zd)),'.','color',colors.grey);
    p_orbit = plot(rowNorm(newData_xyCrossingsOnly(:,c_x:c_z)-[1-secondary.MR,0,0]), abs(newData_xyCrossingsOnly(:,c_zd)),'^','markersize',7,'markeredgecolor',colors.black,'markerfacecolor',colors.cyan);
    p_surf   = plot([secondary.R_n secondary.R_n],[0 0.2],'k', 'linewidth',1.5);
    [legh,objh] = legend([p_L2_1yr, p_orbit, p_surf]  ,label_L2_1yr, 'Single Orbit', 'Europa Surface','FontSize',16);
    % [legh,objh] = legend([p_L2_1yr, p_surf]  ,label_L2_1yr, 'Europa Surface','FontSize',14);
    lineh = findobj(objh,'type','line');
    lineh = lineh([1,2,5,6]);
    % lineh = lineh([1,2,3,4]);
    set(lineh,'linestyle','-','linewidth',3);
    ylim([0 0.09])
%     xlim([0 0.016])
    xlim([0, 0.0205])
    
    % -------------------------------------------------
    %%% Trajectory picked from L2 zdot vs R results
    % -------------------------------------------------
%     newTraj =  nearData_L2_1yr_xyCrossings_dense(nearData_L2_1yr_xyCrossings_dense(:,c_trajID)==trajID_L2_zdot_vs_R,:);
%     figure; hold all
%     PlotBoi3_CR3Bn(26)
%     plot3(newTraj(:,c_x),newTraj(:,c_y),newTraj(:,c_z),'color',colors.blue)

    newTraj_full = data_L2_1yr_xyCrossings_dense(data_L2_1yr_xyCrossings_dense(:,c_trajID)==trajID_L2_zdot_vs_R,:);
    figure; hold all
    ax = gca;
    ax.FontSize;
    ax.FontSize = 14;
    PlotBoi3_CR3Bn(28)
    plot3(newTraj_full(:,c_x),newTraj_full(:,c_y),newTraj_full(:,c_z),'color',colors.blue2)
    plot3(L123(1:2,1),[0,0],[0,0],'^','markeredgecolor',colors.black, 'markerfacecolor', colors.blue)
    xlim([0.976, 1.024])
    ylim([-1 1].*0.015)
    
    % -------------------------------------------------
    %%% 
    % -------------------------------------------------
    
    
    % -------------------------------------------------
    %%% 
    % -------------------------------------------------
    
    
    % -------------------------------------------------
    %%% 
    % -------------------------------------------------
    
    
    % -------------------------------------------------
    %%% 
    % -------------------------------------------------
    


end
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













    