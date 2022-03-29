% ========================================================================
%%% Description
% ========================================================================
% For analyzing the high-latitude-landing-problem data generated in Julia

% Created: 08/13/20
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
rNorm = secondary.a;         % n <-> km
tNorm = 1/secondary.meanMot; % n <-> sec
vNorm = rNorm / tNorm;       % n <-> km/sec

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

%%% Load data
load([dataPath,'lowTrajs_50mps_neckFull_L1_250km_22v0s_4pi_nearEur.mat']); % nearData_L1_4pi
load([dataPath,'lowTrajs_50mps_neckFull_L1_250km_22v0s_1yr_nearEur.mat']); % nearData_L1_1yr
load([dataPath,'lowTrajs_50mps_neckFull_L2_250km_22v0s_4pi_nearEur.mat']); % nearData_L2_4pi
load([dataPath,'lowTrajs_50mps_neckFull_L2_250km_22v0s_1yr_nearEur.mat']); % nearData_L2_1yr
load([dataPath,'lowTrajs_50mps_neckFull_L2_250km_22v0s_1yr_nearEur_xyCrossings.mat']); % nearData_L2_1yr_xyCrossings
load([dataPath,'highTrajs_50mps_70lat_1pi.mat']); % data_high_1pi

%%% Calculate distances from Europa
Rs_L1_4pi = rowNorm(nearData_L1_4pi(:,c_x:c_z) - [1-secondary.MR,0,0]);
Rs_L1_1yr = rowNorm(nearData_L1_1yr(:,c_x:c_z) - [1-secondary.MR,0,0]);
Rs_L2_4pi = rowNorm(nearData_L2_4pi(:,c_x:c_z) - [1-secondary.MR,0,0]);
Rs_L2_1yr = rowNorm(nearData_L2_1yr(:,c_x:c_z) - [1-secondary.MR,0,0]);
Rs_high_1pi = rowNorm(data_high_1pi(:,c_x:c_z) - [1-secondary.MR,0,0]);
Rs_L2_1yr_xyCrossings = rowNorm(nearData_L2_1yr_xyCrossings(:,c_x:c_z) - [1-secondary.MR,0,0]);
Rs_L2_1yr_xyCrossings2 = rowNorm(nearData_L2_1yr_xyCrossings2(:,c_x:c_z) - [1-secondary.MR,0,0]);

%%% Calculate Jacobi constant of dataset
prms.u = secondary.MR;
prms.R2_n = secondary.R_n;
prms.n = 1;
JC_data = getJacobiConstant_ZH(nearData_L1_4pi(1,c_x:c_zd), prms);

%%% Data tags and colors
label_L1_4pi        = 'L1, 4pi';
label_L1_1yr        = 'L1, 1 year';
label_L2_4pi        = 'L2, 4pi';
label_L2_4pi_impact = 'L2, 4pi - impact';
label_L2_1yr        = 'L2, 1 year';
label_high_1pi      = 'Lat 70°+, 1pi';

color_L1_4pi        = colors.grn2;
color_L1_1yr        = colors.ltgrn;
color_L2_4pi        = colors.blue;
color_L2_4pi_impact = colors.brown;
color_L2_1yr        = colors.cyan;
color_high_1pi      = colors.red;

%%% Alternate ranges
smallRange_high_1pi = 1:200:size(data_high_1pi,1);
smallRange_L1_1yr   = 1:100:size(nearData_L1_1yr,1);
smallRange_L2_1yr   = 1:100:size(nearData_L2_1yr,1);

%%% Code used to create the near dataset
if 1+1 == 1
    load([dataPath,'lowTrajs_50mps_neckFull_L1_250km_22v0s_4pi.mat']); % data_L1_4pi
    load([dataPath,'lowTrajs_50mps_neckFull_L1_250km_22v0s_1yr.mat']); % data_L1_1yr
    load([dataPath,'lowTrajs_50mps_neckFull_L2_250km_22v0s_4pi.mat']); % data_L2_4pi
%     load([dataPath,'lowTrajs_50mps_neckImpact_L2_250km_22v0s_4pi.mat']); % data_L2_4pi_impact
    load([dataPath,'lowTrajs_50mps_neckFull_L2_250km_22v0s_1yr.mat']); % data_L2_1yr
    load([dataPath,'highTrajs_50mps_70lat_1pi.mat']); % data_high_1pi
    
    Rs_L1_4pi = rowNorm(data_L1_4pi(:,c_x:c_z) - [1-secondary.MR,0,0]);
    Rs_L1_1yr = rowNorm(data_L1_1yr(:,c_x:c_z) - [1-secondary.MR,0,0]);
    Rs_L2_4pi = rowNorm(data_L2_4pi(:,c_x:c_z) - [1-secondary.MR,0,0]);
    Rs_L2_1yr = rowNorm(data_L2_1yr(:,c_x:c_z) - [1-secondary.MR,0,0]);
    Rs_high_1pi = rowNorm(data_high_1pi(:,c_x:c_z) - [1-secondary.MR,0,0]);
    
    Rmax = max(Rs_high_1pi);
    
    nearData_L1_4pi = data_L1_4pi( Rs_L1_4pi<Rmax,:);
    nearData_L1_1yr = data_L1_1yr( Rs_L1_1yr<Rmax,:);
    nearData_L2_4pi = data_L2_4pi( Rs_L2_4pi<Rmax,:);
    nearData_L2_1yr = data_L2_1yr( Rs_L2_1yr<Rmax,:);
    
    save([dataPath,'lowTrajs_50mps_neckFull_L1_250km_22v0s_4pi_nearEur','.mat'], 'nearData_L1_4pi','-v7.3')
    save([dataPath,'lowTrajs_50mps_neckFull_L1_250km_22v0s_1yr_nearEur','.mat'], 'nearData_L1_1yr','-v7.3')
    save([dataPath,'lowTrajs_50mps_neckFull_L2_250km_22v0s_4pi_nearEur','.mat'], 'nearData_L2_4pi','-v7.3')
    save([dataPath,'lowTrajs_50mps_neckFull_L2_250km_22v0s_1yr_nearEur','.mat'], 'nearData_L2_1yr','-v7.3')
    
    data_high_1pi_xyCrossings = dlmread([dataPath, 'highTrajs_50mps_70lat_1pi_xyCrossings.txt'],',',1,0);
    save([dataPath,'highTrajs_50mps_70lat_1pi_xyCrossings','.mat'], 'data_high_1pi_xyCrossings','-v7.3')
%     load([dataPath,'highTrajs_50mps_70lat_1pi_xyCrossings.mat']); % data_L2_4pi
    
    data_L2_1yr_xyCrossings2 = dlmread([dataPath, 'lowTrajs_50mps_neckFull_L2_250km_22v0s_1yr_xyCrossings2.txt'],',',1,0);
%     
    Rs_L2_1yr_xyCrossings2 = rowNorm(data_L2_1yr_xyCrossings2(:,c_x:c_z) - [1-secondary.MR,0,0]);
    nearData_L2_1yr_xyCrossings2 = data_L2_1yr_xyCrossings2( Rs_L2_1yr_xyCrossings2<Rmax,:);
%     save([dataPath,'lowTrajs_50mps_neckFull_L2_250km_22v0s_1yr_nearEur_xyCrossings','.mat'], 'nearData_L2_1yr_xyCrossings','-v7.3')

end

% traj1 = data_L2_1yr(data_L2_1yr(:,1)==400, :);
% traj2 = data_L2_1yr_xyCrossings(data_L2_1yr_xyCrossings(:,1)==400, :);



% ========================================================================
%%% Plots
% ========================================================================
% -------------------------------------------------
%%% Plot the trajectories
% -------------------------------------------------
if 1+1==1
    figure; hold all
    PlotBoi3_CR3Bn(23)
    axis equal
%     xlim([L123(1,1)-secondary.R_n, L123(2,1)+secondary.R_n])
    
%     p_high_1pi = plot3(data_high_1pi(:,c_x),data_high_1pi(:,c_y),data_high_1pi(:,c_z), '.','markersize',1, 'color', color_high_1pi);
    p_L1_1yr   = plot3(nearData_L1_1yr(:,c_x),nearData_L1_1yr(:,c_y),nearData_L1_1yr(:,c_z), '.','markersize',1, 'color', color_L1_1yr);
    p_L2_1yr   = plot3(nearData_L2_1yr(:,c_x),nearData_L2_1yr(:,c_y),nearData_L2_1yr(:,c_z), '.','markersize',1, 'color', color_L2_1yr);

    [legh,objh] = legend([p_L1_1yr, p_L2_1yr]  ,label_L1_1yr, label_L2_1yr);
%     [legh,objh] = legend([p_high_1pi, p_L1_1yr, p_L2_1yr], label_high_1pi  ,label_L1_1yr, label_L2_1yr);
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);
    

    plotCR3BP_Neck(prms, L123, JC_data, 600, 200, colors.black, 1.5)
    plotCR3BP_YZNeck(JC_data, prms.u , 1, 0, prms, colors.black, 1.5)
    plotCR3BP_YZNeck(JC_data, prms.u , 2, 0, prms, colors.black, 1.5)

end

% -------------------------------------------------
%%% Latitude vs Longitude
% -------------------------------------------------
if 1+1==1
    figure; hold all
    PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$',23,'LaTex')
    ylim([-90 90])
    xlim([-180 180])
    
    p_high_1pi = plot(data_high_1pi(1:200:size(data_high_1pi,1),c_lon),data_high_1pi(1:200:size(data_high_1pi,1),c_lat), '.','markersize',1, 'color', color_high_1pi);
    p_L1_1yr   = plot(nearData_L1_1yr(:,c_lon),nearData_L1_1yr(:,c_lat), '.','markersize',1, 'color', color_L1_1yr);
    p_L2_1yr   = plot(nearData_L2_1yr(:,c_lon),nearData_L2_1yr(:,c_lat), '.','markersize',1, 'color', color_L2_1yr);
%     
%     [legh,objh] = legend([p_L1_1yr, p_L2_1yr]  ,label_L1_1yr, label_L2_1yr);
    [legh,objh] = legend([p_high_1pi, p_L1_1yr, p_L2_1yr], label_high_1pi  ,label_L1_1yr, label_L2_1yr);
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);
end

% -------------------------------------------------
%%% Look at zdd_secondary
% -------------------------------------------------
if 1+1==1
    zdds_secondary_L2_1yr = zeros(size(nearData_L2_1yr,1),1);
    for kk = 1:length(zdds_secondary_L2_1yr)
        zdds_secondary_L2_1yr(kk) = -(prms.u * nearData_L2_1yr(kk,c_z)) / (Rs_L2_1yr(kk)^3);
    end
    
    %%% zdd_secondary vs R_n
    figure; hold all
    PlotBoi2('$R_n$','$-\frac{\mu z}{R_n^3}$',23,'LaTex')
    plot(Rs_L2_1yr, zdds_secondary_L2_1yr, '.', 'color', color_L2_1yr);
    
    %%% Latitude vs zdd_secondary
    figure; hold all
    PlotBoi2('$-\frac{\mu z}{R_n^3}$', 'Latitude, $^\circ$', 23, 'LaTex')
    plot(zdds_secondary_L2_1yr, nearData_L2_1yr(:,c_lat), '.', 'color', color_L2_1yr);
    
    %%% Latitude vs R_n
    figure; hold all
    PlotBoi2('$R_n$','Latitude, $^\circ$',23,'LaTex')
    plot(Rs_L2_1yr, nearData_L2_1yr(:,c_lat), '.', 'color', color_L2_1yr);
end

% -------------------------------------------------
%%% Inclination SCI
% -------------------------------------------------
if 1+1==1
    %%% inclination vs R
    figure; hold all
    PlotBoi2('$R_n$','$i_{SCI}$',23,'LaTex')
    p_high_1pi = plot(Rs_high_1pi, data_high_1pi(:,c_i), '.','color', color_high_1pi);
%     p_L1_1yr   = plot(Rs_L1_1yr, nearData_L1_1yr(:,c_i), '.','color', color_L1_1yr);
    p_L2_1yr   = plot(Rs_L2_1yr, nearData_L2_1yr(:,c_i), '.','color', color_L2_1yr);
%     [legh,objh] = legend([p_L1_1yr, p_L2_1yr]  ,label_L1_1yr, label_L2_1yr);
%     [legh,objh] = legend([p_high_1pi, p_L1_1yr, p_L2_1yr], label_high_1pi  ,label_L1_1yr, label_L2_1yr);
    [legh,objh] = legend([p_high_1pi, p_L2_1yr], label_high_1pi , label_L2_1yr);
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);
    
    
    %%% inclination vs latitude
    figure; hold all
    PlotBoi2('Latitude, $^\circ$','$i_{SCI}$',23,'LaTex')
    p_high_1pi = plot(data_high_1pi(:,c_lat), data_high_1pi(:,c_i), '.','color', color_high_1pi);
%     p_L1_1yr   = plot(nearData_L1_1yr(:,c_lat), nearData_L1_1yr(:,c_i), '.','color', color_L1_1yr);
    p_L2_1yr   = plot(nearData_L2_1yr(:,c_lat), nearData_L2_1yr(:,c_i), '.','color', color_L2_1yr);
%     [legh,objh] = legend([p_L1_1yr, p_L2_1yr]  ,label_L1_1yr, label_L2_1yr);
%     [legh,objh] = legend([p_high_1pi, p_L1_1yr, p_L2_1yr], label_high_1pi  ,label_L1_1yr, label_L2_1yr);
    [legh,objh] = legend([p_high_1pi, p_L2_1yr], label_high_1pi , label_L2_1yr);
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);
    
    
    %%% inclination vs latitude vs longitude
    figure; hold all
    PlotBoi3('Longitude, $^\circ$', 'Latitude, $^\circ$','$i_{SCI}$',23,'LaTex')
    p_high_1pi = plot3(data_high_1pi(smallRange_high_1pi,c_lon), data_high_1pi(smallRange_high_1pi,c_lat), data_high_1pi(smallRange_high_1pi,c_i), '.','color', color_high_1pi);
%     p_L1_1yr   = plot3(nearData_L1_1yr(:,c_lon), nearData_L1_1yr(:,c_lat), nearData_L1_1yr(:,c_i), '.','color', color_L1_1yr);
    p_L2_1yr   = plot3(nearData_L2_1yr(:,c_lon), nearData_L2_1yr(:,c_lat), nearData_L2_1yr(:,c_i), '.','color', color_L2_1yr);
%     [legh,objh] = legend([p_L1_1yr, p_L2_1yr]  ,label_L1_1yr, label_L2_1yr);
%     [legh,objh] = legend([p_high_1pi, p_L1_1yr, p_L2_1yr], label_high_1pi  ,label_L1_1yr, label_L2_1yr);
    [legh,objh] = legend([p_high_1pi, p_L2_1yr], label_high_1pi , label_L2_1yr);
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);
end

% -------------------------------------------------
%%% R vs z vs latitude
% -------------------------------------------------
if 1+1 == 1
    figure; hold all
    PlotBoi3('$R_n$', '$z$','Latitude, $^\circ$',23,'LaTex')

    p_high = plot3(Rs_high_1pi(smallRange_high_1pi), data_high_1pi(smallRange_high_1pi,c_z),   data_high_1pi(smallRange_high_1pi,c_lat),'.','color',color_high_1pi);
%     p_L1   = plot3(Rs_L1_1yr,   nearData_L1_1yr(:,c_z), nearData_L1_1yr(:,c_lat),'.','color',color_L1_1yr);
    p_L2   = plot3(Rs_L2_1yr,   nearData_L2_1yr(:,c_z), nearData_L2_1yr(:,c_lat),'.','color',color_L2_1yr);

    [legh,objh] = legend([p_high, p_L2], label_high_1pi, label_L2_1yr);
%     [legh,objh] = legend([p_high, p_L1, p_L2], label_high_1pi, label_L1_1yr, label_L2_1yr);
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);
end


% -------------------------------------------------
%%% inclination SCR - takes 10 minutes to compute
% -------------------------------------------------
if 1+1 == 1
    is_SCR_L1_yr = NaN(size(nearData_L1_1yr,1),1);
    for kk = 1:length(is_SCR_L1_yr)
        [a,e,i,raan,w,ta] = ECI2OE(nearData_L1_1yr(kk,c_x:c_z)-[1-secondary.MR,0,0],nearData_L1_1yr(kk,c_xd:c_zd),secondary.MR);
        is_SCR_L1_yr(kk) = i*180/pi;
    end
    is_SCR_L2_yr = NaN(size(nearData_L2_1yr,1),1);
    for kk = 1:length(is_SCR_L2_yr)
        [a,e,i,raan,w,ta] = ECI2OE(nearData_L2_1yr(kk,c_x:c_z)-[1-secondary.MR,0,0],nearData_L2_1yr(kk,c_xd:c_zd),secondary.MR);
        is_SCR_L2_yr(kk) = i*180/pi;
    end
    is_SCR_high_1pi = NaN(size(data_high_1pi,1),1);
    for kk = 1:length(is_SCR_high_1pi)
        [a,e,i,raan,w,ta] = ECI2OE(data_high_1pi(kk,c_x:c_z)-[1-secondary.MR,0,0],data_high_1pi(kk,c_xd:c_zd),secondary.MR);
        is_SCR_high_1pi(kk) = i*180/pi;
    end
    
    %%% inclination vs latitude vs longitude
    figure; hold all
    PlotBoi3('Longitude, $^\circ$', 'Latitude, $^\circ$','$i_{SCR}$',23,'LaTex')
    p_high_1pi = plot3(data_high_1pi(smallRange_high_1pi,c_lon), data_high_1pi(smallRange_high_1pi,c_lat), is_SCR_high_1pi(smallRange_high_1pi), '.','color', color_high_1pi);
%     p_L1_1yr   = plot3(nearData_L1_1yr(:,c_lon), nearData_L1_1yr(:,c_lat), is_SCR_L1_yr, '.','color', color_L1_1yr);
    p_L2_1yr   = plot3(nearData_L2_1yr(:,c_lon), nearData_L2_1yr(:,c_lat), is_SCR_L2_yr, '.','color', color_L2_1yr);
%     [legh,objh] = legend([p_L1_1yr, p_L2_1yr]  ,label_L1_1yr, label_L2_1yr);
%     [legh,objh] = legend([p_high_1pi, p_L1_1yr, p_L2_1yr], label_high_1pi  ,label_L1_1yr, label_L2_1yr);
    [legh,objh] = legend([p_high_1pi, p_L2_1yr], label_high_1pi , label_L2_1yr);
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);
    
    
    %%% a vs latitude vs longitude
    figure; hold all
    PlotBoi3('Longitude, $^\circ$', 'Latitude, $^\circ$','$a$',23,'LaTex')
    p_high_1pi = plot3(data_high_1pi(smallRange_high_1pi,c_lon), data_high_1pi(smallRange_high_1pi,c_lat), data_high_1pi(smallRange_high_1pi, c_a), '.','color', color_high_1pi);
%     p_L1_1yr   = plot3(nearData_L1_1yr(:,c_lon), nearData_L1_1yr(:,c_lat), nearData_L1_1yr(:, c_ia, '.','color', color_L1_1yr);
    p_L2_1yr   = plot3(nearData_L2_1yr(:,c_lon), nearData_L2_1yr(:,c_lat), nearData_L2_1yr(:,c_a), '.','color', color_L2_1yr);
%     [legh,objh] = legend([p_L1_1yr, p_L2_1yr]  ,label_L1_1yr, label_L2_1yr);
%     [legh,objh] = legend([p_high_1pi, p_L1_1yr, p_L2_1yr], label_high_1pi  ,label_L1_1yr, label_L2_1yr);
    [legh,objh] = legend([p_high_1pi, p_L2_1yr], label_high_1pi , label_L2_1yr);
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);
end


if 1+1==1
    figure; hold all
    PlotBoi2('$R_n$','$\dot{z}$',23,'LaTex')
    p_high_1pi = plot(Rs_high_1pi(data_high_1pi(:,c_isxyCross)==1), data_high_1pi(data_high_1pi(:,c_isxyCross)==1, c_zd), '.','color', color_high_1pi);
    p_L2_1yr   = plot(Rs_L2_1yr(nearData_L2_1yr(:,c_isxyCross)==1), nearData_L2_1yr(nearData_L2_1yr(:,c_isxyCross)==1, c_zd), '.','color', color_L2_1yr);
    [legh,objh] = legend([p_high_1pi, p_L2_1yr], label_high_1pi , label_L2_1yr);
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);
end

% -------------------------------------------------
%%% Hxy
% -------------------------------------------------
if 1+1==1
    %%% Calculate HxHx data
    HxHy_high_1pi = rowNorm(data_high_1pi(:,c_hx:c_hy));
    HxHy_L1_1yr   = rowNorm(nearData_L1_1yr(:,c_hx:c_hy));
    HxHy_L2_1yr   = rowNorm(nearData_L2_1yr(:,c_hx:c_hy));
    
    %%% Hxy vs H
    figure; hold all
    PlotBoi2('$|H|$','$|H_xH_y|$',23,'LaTex')
%     p_L1_1yr   = plot(nearData_L1_1yr(:,c_H),HxHy_L1_1yr, '.','markersize',1, 'color', color_L1_1yr);
    p_L2_1yr   = plot(nearData_L2_1yr(:,c_H),HxHy_L2_1yr, '.','markersize',1, 'color', color_L2_1yr);
    p_high_1pi = plot(data_high_1pi(smallRange_high_1pi,c_H),HxHy_high_1pi(smallRange_high_1pi), '.','markersize',1, 'color', color_high_1pi);
%     [legh,objh] = legend([p_high_1pi, p_L1_1yr, p_L2_1yr], label_high_1pi  ,label_L1_1yr, label_L2_1yr);
    [legh,objh] = legend([p_high_1pi, p_L2_1yr], label_high_1pi , label_L2_1yr);
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);
    
    
    %%% Hxy/H vs R
    figure; hold all
    PlotBoi2('$R_n$','$\frac{H_{xy}}{H}$',23,'LaTex')
%     p_L1_1yr   = plot(Rs_L1_1yr,HxHy_L1_1yr./nearData_L1_1yr(:,c_H), '.','markersize',1, 'color', color_L1_1yr);
    p_L2_1yr   = plot(Rs_L2_1yr,HxHy_L2_1yr./nearData_L2_1yr(:,c_H), '.','markersize',1, 'color', color_L2_1yr);
    p_high_1pi = plot(Rs_high_1pi(smallRange_high_1pi),HxHy_high_1pi(smallRange_high_1pi)./data_high_1pi(smallRange_high_1pi,c_H), '.','markersize',1, 'color', color_high_1pi);
%     [legh,objh] = legend([p_high_1pi, p_L1_1yr, p_L2_1yr], label_high_1pi  ,label_L1_1yr, label_L2_1yr);
    [legh,objh] = legend([p_high_1pi, p_L2_1yr], label_high_1pi , label_L2_1yr);
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);
    
    %%% z vs R
    figure; hold all
    PlotBoi2('Latitude, $^\circ$','$H_z$',23,'LaTex')
%     p_L1_1yr   = plot(Rs_L1_1yr,nearData_L1_1yr(:,c_hz), '.','markersize',1, 'color', color_L1_1yr);
    p_L2_1yr   = plot(nearData_L2_1yr(:,c_z),nearData_L2_1yr(:,c_hz), '.','markersize',1, 'color', color_L2_1yr);
    p_high_1pi = plot(data_high_1pi(smallRange_high_1pi, c_z),data_high_1pi(smallRange_high_1pi,c_hz), '.','markersize',1, 'color', color_high_1pi);
%     [legh,objh] = legend([p_high_1pi, p_L1_1yr, p_L2_1yr], label_high_1pi  ,label_L1_1yr, label_L2_1yr);
    [legh,objh] = legend([p_high_1pi, p_L2_1yr], label_high_1pi , label_L2_1yr);
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);
end








% -------------------------------------------------
%%% R vs zd for xyCrossing
% -------------------------------------------------
if 1+1==1
    high_xyCrossingIndicies = data_high_1pi(:,c_isxyCross)==1;
    L2_xyCrossingIndicies = nearData_L2_1yr_xyCrossings(:,c_isxyCross)==1;
    L2_xyCrossingIndicies2 = nearData_L2_1yr_xyCrossings2(:,c_isxyCross)==1;
    
    figure; hold all
    PlotBoi2('$R_n$','$\dot{z}$',23,'LaTex')
%     p_high_1pi = plot(Rs_high_1pi(high_xyCrossingIndicies), data_high_1pi(high_xyCrossingIndicies,c_zd),'.','color',color_high_1pi);
%     p_L2_1yr   = plot(Rs_L2_1yr_xyCrossings(L2_xyCrossingIndicies), nearData_L2_1yr_xyCrossings(L2_xyCrossingIndicies,c_zd),'.','color',color_L2_1yr);
    p_high_1pi = plot(Rs_high_1pi(high_xyCrossingIndicies), abs(data_high_1pi(high_xyCrossingIndicies,c_zd)),'.','color',color_high_1pi);
%     p_L2_1yr   = plot(Rs_L2_1yr_xyCrossings(L2_xyCrossingIndicies), abs(nearData_L2_1yr_xyCrossings(L2_xyCrossingIndicies,c_zd)),'.','color',color_L2_1yr);    [legh,objh] = legend([p_high_1pi, p_L2_1yr], label_high_1pi , label_L2_1yr);
    p_L2_1yr   = plot(Rs_L2_1yr_xyCrossings2(L2_xyCrossingIndicies2), abs(nearData_L2_1yr_xyCrossings2(L2_xyCrossingIndicies2,c_zd)),'.','color',color_L2_1yr);    
    [legh,objh] = legend([p_high_1pi, p_L2_1yr], label_high_1pi , label_L2_1yr);
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);
%     setLogPlot('y')


    figure; hold all
%     PlotBoi2('$R_n$','$\dot{z}$',23,'LaTex')
    PlotBoi3('$R_n$','$\dot{z}$', 'Longitude, $^\circ$',23,'LaTex')
%     p_high_1pi = plot(Rs_high_1pi(high_xyCrossingIndicies), data_high_1pi(high_xyCrossingIndicies,c_zd),'.','color',color_high_1pi);
%     p_L2_1yr   = plot(Rs_L2_1yr_xyCrossings(L2_xyCrossingIndicies), nearData_L2_1yr_xyCrossings(L2_xyCrossingIndicies,c_zd),'.','color',color_L2_1yr);
    p_high_1pi = plot3(Rs_high_1pi(high_xyCrossingIndicies), abs(data_high_1pi(high_xyCrossingIndicies,c_zd)), data_high_1pi(high_xyCrossingIndicies,c_lon),'.','color',color_high_1pi);
%     p_L2_1yr   = plot3(Rs_L2_1yr_xyCrossings(L2_xyCrossingIndicies), abs(nearData_L2_1yr_xyCrossings(L2_xyCrossingIndicies,c_zd)), nearData_L2_1yr_xyCrossings(L2_xyCrossingIndicies,c_lon),'.','color',color_L2_1yr);    [legh,objh] = legend([p_high_1pi, p_L2_1yr], label_high_1pi , label_L2_1yr);
    p_L2_1yr   = plot3(Rs_L2_1yr_xyCrossings2(L2_xyCrossingIndicies2), abs(nearData_L2_1yr_xyCrossings2(L2_xyCrossingIndicies2,c_zd)), nearData_L2_1yr_xyCrossings2(L2_xyCrossingIndicies2,c_lon),'.','color',color_L2_1yr);    [legh,objh] = legend([p_high_1pi, p_L2_1yr], label_high_1pi , label_L2_1yr);
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);
    
    figure; hold all
    PlotBoi2('$R_n$','$-\frac{\mu z}{R_n^3}$',23,'LaTex')
    plot(Rs_L2_1yr, zdds_secondary_L2_1yr,'.','color',color_L2_1yr)
    ylim([-1 1].*0.3)
end



if 2+2==2
    data_L2_xyCrossingsOnly = data_L2_1yr_xyCrossings2(data_L2_1yr_xyCrossings2(:,c_isxyCross)==1,:);
    Rs_L2_xyCrossingsOnly = rowNorm(data_L2_xyCrossingsOnly(:,c_x:c_z) - [1-secondary.MR,0,0]);
    
    data_high_xyCrossingsOnly = data_high_1pi(data_high_1pi(:,c_isxyCross)==1,:);
    Rs_high_xyCrossingsOnly = rowNorm(data_high_xyCrossingsOnly(:,c_x:c_z) - [1-secondary.MR,0,0]);
    
    
    figure; hold all
    PlotBoi2('$R_n$','$\dot{z}$',23,'LaTex')
    xlim([0 0.006])
    ylim([-1 1].*0.5)
    plot(Rs_L2_xyCrossingsOnly, abs(data_L2_xyCrossingsOnly(:,c_zd)),'.','color',color_L2_1yr);
    plot(Rs_high_xyCrossingsOnly, abs(data_high_xyCrossingsOnly(:,c_zd)),'.','color',color_high_1pi)
    plot(Rs_L2_xyCrossingsOnly(data_L2_xyCrossingsOnly(:,c_trajID)==989), abs(data_L2_xyCrossingsOnly(data_L2_xyCrossingsOnly(:,c_trajID)==989,c_zd)), 'x','color',colors.mag);
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






