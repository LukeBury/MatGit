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
load_fullData       = false;
load_nearData       = false;
load_xyCrossingData = true;

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
%%% Columb specifiers
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
if load_fullData
    load([dataPath,'lowTrajs_50mps_neckFull_L1_250km_22v0s_4pi.mat']); % data_L1_4pi
    load([dataPath,'lowTrajs_50mps_neckFull_L1_250km_22v0s_1yr.mat']); % data_L1_1yr
    load([dataPath,'lowTrajs_50mps_neckFull_L2_250km_22v0s_4pi.mat']); % data_L2_4pi
    load([dataPath,'lowTrajs_50mps_neckImpact_L2_250km_22v0s_4pi.mat']); % data_L2_4pi_impact
    load([dataPath,'lowTrajs_50mps_neckFull_L2_250km_22v0s_1yr.mat']); % data_L2_1yr
    load([dataPath,'highTrajs_50mps_70lat_1pi.mat']); % data_high_1pi
    
    %%% Calculate distances
    Rs_L1_4pi        = rowNorm(data_L1_4pi(:,c_x:c_z) - [1-secondary.MR,0,0]);
    Rs_L1_1yr        = rowNorm(data_L1_1yr(:,c_x:c_z) - [1-secondary.MR,0,0]);
    Rs_L2_4pi        = rowNorm(data_L2_4pi(:,c_x:c_z) - [1-secondary.MR,0,0]);
    Rs_L2_4pi_impact = rowNorm(data_L2_4pi_impact(:,c_x:c_z) - [1-secondary.MR,0,0]);
    Rs_L2_1yr        = rowNorm(data_L2_1yr(:,c_x:c_z) - [1-secondary.MR,0,0]);
    Rs_high_1pi      = rowNorm(data_high_1pi(:,c_x:c_z) - [1-secondary.MR,0,0]);
    
    
    %%% Indicies for near data
    Rmax = max(Rs_high_1pi);
    indicies_L1_near   = find(Rs_L1_1yr < Rmax);
    indicies_L2_near   = find(Rs_L2_1yr < Rmax);
    
    %%% Calculate Jacobi constant of dataset
    JC_data = getJacobiConstant_ZH(data_high_1pi(1,c_x:c_zd), prms);
end

if load_nearData
    %%% The distance of each state from Europa is less than the max
    %%% distance of the high set
    load([dataPath,'lowTrajs_50mps_neckFull_L1_250km_22v0s_4pi_nearEur.mat']); % nearData_L1_4pi
    load([dataPath,'lowTrajs_50mps_neckFull_L1_250km_22v0s_1yr_nearEur.mat']); % nearData_L1_1yr
    load([dataPath,'lowTrajs_50mps_neckFull_L2_250km_22v0s_4pi_nearEur.mat']); % nearData_L2_4pi
    load([dataPath,'lowTrajs_50mps_neckFull_L2_250km_22v0s_1yr_nearEur.mat']); % nearData_L2_1yr
    load([dataPath,'highTrajs_50mps_70lat_1pi.mat']); % data_high_1pi
    
    %%% Calculate distances
    Rs_L1_4pi = rowNorm(nearData_L1_4pi(:,c_x:c_z) - [1-secondary.MR,0,0]);
    Rs_L1_1yr = rowNorm(nearData_L1_1yr(:,c_x:c_z) - [1-secondary.MR,0,0]);
    Rs_L2_4pi = rowNorm(nearData_L2_4pi(:,c_x:c_z) - [1-secondary.MR,0,0]);
    Rs_L2_1yr = rowNorm(nearData_L2_1yr(:,c_x:c_z) - [1-secondary.MR,0,0]);
    Rs_high_1pi = rowNorm(data_high_1pi(:,c_x:c_z) - [1-secondary.MR,0,0]);
    
    %%% Calculate Jacobi constant of dataset
    JC_data = getJacobiConstant_ZH(data_high_1pi(1,c_x:c_zd), prms);
end

if load_xyCrossingData
    % Full data, but the only events tracked are x-y-plane crossings
    load([dataPath,'lowTrajs_50mps_neckFull_L1_250km_22v0s_1yr_xyCrossings.mat']); % data_L1_1yr_xyCrossings
    load([dataPath,'lowTrajs_50mps_neckFull_L2_250km_22v0s_1yr_xyCrossings.mat']); % data_L2_1yr_xyCrossings
%     load([dataPath,'lowTrajs_50mps_neckFull_L2_125km_22v0s_1yr_xyCrossings.mat']); % data_L2_1yr_xyCrossings
    load([dataPath,'highTrajs_50mps_70lat_1pi_xyCrossings.mat']); % data_high_1pi_xyCrossings
    
    %%% Calculate distances
    Rs_L1_1yr = rowNorm(data_L1_1yr_xyCrossings(:,c_x:c_z) - [1-secondary.MR,0,0]);
    Rs_L2_1yr = rowNorm(data_L2_1yr_xyCrossings(:,c_x:c_z) - [1-secondary.MR,0,0]);
    Rs_high_1pi = rowNorm(data_high_1pi_xyCrossings(:,c_x:c_z) - [1-secondary.MR,0,0]);
    
    %%% Determine indicies for x-y-plane crossing events
    indicies_L1_xyCrossingEvents = find(data_L1_1yr_xyCrossings(:,c_isxyCross)==1);
    indicies_L2_xyCrossingEvents = find(data_L2_1yr_xyCrossings(:,c_isxyCross)==1);
    indicies_high_xyCrossingEvents = find(data_high_1pi_xyCrossings(:,c_isxyCross)==1);
    
    %%% Determine indicies data near Europa
    indicies_L1_xyCrossings_near = find(Rs_L1_1yr <= 0.02);
    indicies_L2_xyCrossings_near = find(Rs_L2_1yr <= 0.02);
    
    %%% Data sets for x-y crossings only
    data_L1_1yr_xyCrossingsOnly   = data_L1_1yr_xyCrossings(data_L1_1yr_xyCrossings(:,c_isxyCross)==1,:);
    data_L2_1yr_xyCrossingsOnly   = data_L2_1yr_xyCrossings(data_L2_1yr_xyCrossings(:,c_isxyCross)==1,:);
    data_high_1pi_xyCrossingsOnly = data_high_1pi_xyCrossings(data_high_1pi_xyCrossings(:,c_isxyCross)==1,:);
    
    %%% Alternate ranges
    smallRange_high_1pi = 1:200:size(data_high_1pi_xyCrossings,1);
    smallRange_L1_1yr   = 1:100:size(data_L1_1yr_xyCrossings,1);
    smallRange_L2_1yr   = 1:100:size(data_L2_1yr_xyCrossings,1);
    
    %%% Calculate Jacobi constant of dataset
    prms.u = secondary.MR;
    prms.R2_n = secondary.R_n;
    prms.n = 1;
    JC_data = getJacobiConstant_ZH(data_high_1pi_xyCrossings(1,c_x:c_zd), prms);

end

if 1+1==1
    % load([dataPath,'lowTrajs_50mps_neckFull_L2_250km_22v0s_1yr_nearEur_xyCrossings.mat']); % nearData_L2_1yr_xyCrossings
end

%%% Set prms
prms.u = secondary.MR;
prms.R2_n = secondary.R_n;
prms.n = 1;

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



% ========================================================================
%%% Plotting
% ========================================================================

%%% For load_xyCrossingData
if 1+1==1
    %%% zd vs R
    figure; hold all
    PlotBoi2('$R_n$','$\dot{z}_n$',23,'LaTex')
    xlim([0 (L123(2,1) - (1-prms.u))]);
    ylim([-1 1].*0.2)
    
    plot(Rs_high_1pi(indicies_high_xyCrossingEvents), data_high_1pi_xyCrossingsOnly(:,c_zd), '.', 'color', color_high_1pi);
    plot(Rs_L1_1yr(indicies_L1_xyCrossingEvents), data_L1_1yr_xyCrossingsOnly(:,c_zd), '.', 'color', color_L1_1yr);
    plot(Rs_L2_1yr(indicies_L2_xyCrossingEvents), data_L2_1yr_xyCrossingsOnly(:,c_zd), '.', 'color', color_L2_1yr);
    
%     plot(rowNorm(maxZTraj_L1(maxZTraj_L1(:,c_isxyCross)==1, c_x:c_z) - [1-prms.u,0,0]), maxZTraj_L1(maxZTraj_L1(:,c_isxyCross)==1, c_zd),'k.')
    plot(rowNorm(X_event(:, 1:3) - [1-prms.u,0,0]), X_event(:, 6),'m.')
    
    
    
    
    %%% i vs R
    figure; hold all
    PlotBoi2('$R_n$','$i_{SCI}$, $^\circ$',23,'LaTex')
    xlim([0 (L123(2,1) - (1-prms.u))]);
    
    plot(Rs_L1_1yr(indicies_L1_xyCrossingEvents), data_L1_1yr_xyCrossingsOnly(:,c_i), '.', 'color', color_L1_1yr);
    
    X_SCI_event = X_BaCR2SCI(X_event, time_event, prms);
    a_e_i_raan_w_ta_event = zeros(size(X_SCI_event,1),6);
    for kk = 1:size(X_SCI_event,1)
        [a,e,i,raan,w,ta] = ECI2OE(X_SCI_event(kk,1:3), X_SCI_event(kk,4:6),prms.u);
        a_e_i_raan_w_ta_event(kk,:) = [a,e,i,raan,w,ta];
    end
    
    plot(rowNorm(X_event(:, 1:3) - [1-prms.u,0,0]), a_e_i_raan_w_ta_event(:, 3).*180/pi,'k.')
end

%%% Playing with zdot vs R
if 1+1==1
    zdRatio_L1 = zeros(size(data_L1_1yr_xyCrossingsOnly,1),1);
    zdRatio_L2 = zeros(size(data_L2_1yr_xyCrossingsOnly,1),1);
    zdRatio_high = zeros(size(data_high_1pi_xyCrossingsOnly,1),1);
    
    for kk = 1:length(zdRatio_L1)
        zdRatio_L1(kk) = data_L1_1yr_xyCrossingsOnly(kk,c_zd) / norm(data_L1_1yr_xyCrossingsOnly(kk,c_xd:c_zd));
    end
    for kk = 1:length(zdRatio_L2)
        zdRatio_L2(kk) = data_L2_1yr_xyCrossingsOnly(kk,c_zd) / norm(data_L2_1yr_xyCrossingsOnly(kk,c_xd:c_zd));
    end
    for kk = 1:length(zdRatio_high)
        zdRatio_high(kk) = data_high_1pi_xyCrossingsOnly(kk,c_zd) / norm(data_high_1pi_xyCrossingsOnly(kk,c_xd:c_zd));
    end
    
    %%% zd/v vs R
    figure; hold all
    PlotBoi2('$R_n$','$\frac{\dot{z}_n}{v_n}$',23,'LaTex')
    xlim([0 (L123(2,1) - (1-prms.u))]);
%     ylim([-1 1].*0.2)
    
    plot(Rs_L1_1yr(indicies_L1_xyCrossingEvents), zdRatio_L1, '.', 'color', color_L1_1yr);
    plot(Rs_L2_1yr(indicies_L2_xyCrossingEvents), zdRatio_L2, '.', 'color', color_L2_1yr);
    plot(Rs_high_1pi(indicies_high_xyCrossingEvents), zdRatio_high, '.', 'color', color_high_1pi);
    
end

%%% Looking at extreme orbits
if 1+1==1
    %%% max latitude traj from L1
    maxLatIndex_L1 = find(abs(data_L1_1yr_xyCrossings(:,c_lat)) == max(abs(data_L1_1yr_xyCrossings(:,c_lat))));
    maxLatTrajID_L1 = data_L1_1yr_xyCrossings(maxLatIndex_L1, c_trajID);
    maxLatTraj_L1 = data_L1_1yr_xyCrossings(data_L1_1yr_xyCrossings(:,c_trajID)==maxLatTrajID_L1(1), :);
    
    %%% max latitude traj from L2
    maxLatIndex_L2 = find(abs(data_L2_1yr_xyCrossings(:,c_lat)) == max(abs(data_L2_1yr_xyCrossings(:,c_lat))));
    maxLatTrajID_L2 = data_L2_1yr_xyCrossings(maxLatIndex_L2, c_trajID);
    maxLatTraj_L2 = data_L2_1yr_xyCrossings(data_L2_1yr_xyCrossings(:,c_trajID)==maxLatTrajID_L2(1), :);
    
    %%% max near z traj from L1
    maxZIndex_L1 = find(abs(data_L1_1yr_xyCrossings(:,c_z)) == max(abs(data_L1_1yr_xyCrossings(indicies_L1_xyCrossings_near,c_z))));
    maxZTrajID_L1 = data_L1_1yr_xyCrossings(maxZIndex_L1, c_trajID);
    maxZTraj_L1   = data_L1_1yr_xyCrossings(data_L1_1yr_xyCrossings(:,c_trajID)==maxZTrajID_L1(1), :);
%     newTraj1 = data_L1_1yr_xyCrossings(data_L1_1yr_xyCrossings(:,c_trajID)==222, :);
%     newTraj2 = data_L1_1yr_xyCrossings(data_L1_1yr_xyCrossings(:,c_trajID)==2587, :); % Weaker manifold
%     newTraj3 = data_L1_1yr_xyCrossings(data_L1_1yr_xyCrossings(:,c_trajID)==1035, :); % smaller center radius
%     newTraj4 = data_L1_1yr_xyCrossings(data_L1_1yr_xyCrossings(:,c_trajID)==2596, :); % smaller center radius
%     newTraj5 = data_L1_1yr_xyCrossings(data_L1_1yr_xyCrossings(:,c_trajID)==155, :);
%     newTraj6 = data_L1_1yr_xyCrossings(data_L1_1yr_xyCrossings(:,c_trajID)==1099, :);

%     testData = data_L1_1yr_xyCrossingsOnly((rowNorm(data_L1_1yr_xyCrossingsOnly(:,c_x:c_z) - [1-prms.u,0,0]) < (L123(2,1)-(1-prms.u))), :);
%     testData = testData(rowNorm(testData(:,c_x:c_z)-[1-prms.u,0,0])<(0.01738), :);
%     testData = testData(rowNorm(testData(:,c_x:c_z)-[1-prms.u,0,0])>(0.01737), :);
%     testData = testData(testData(:,c_zd)>0.0106, :);
%     testData = testData(testData(:,c_zd)<0.0107, :);

    %%% max near z traj from L2
    maxZIndex_L2 = find(abs(data_L2_1yr_xyCrossings(:,c_z)) == max(abs(data_L2_1yr_xyCrossings(indicies_L2_xyCrossings_near,c_z))));
    maxZTrajID_L2 = data_L2_1yr_xyCrossings(maxZIndex_L2, c_trajID);
    maxZTraj_L2 = data_L2_1yr_xyCrossings(data_L2_1yr_xyCrossings(:,c_trajID)==maxZTrajID_L2(1), :);
    
    
    
    %%% Plotting maxLatL1
    figure; hold all
    title('Max lat - L1')
    PlotBoi3_CR3Bn(23)
    axis equal
    plot3(maxLatTraj_L1(:,c_x),maxLatTraj_L1(:,c_y),maxLatTraj_L1(:,c_z));
    
    %%% Plotting maxLatL2
    figure; hold all
    title('Max lat - L2')
    PlotBoi3_CR3Bn(23)
    axis equal
    plot3(maxLatTraj_L2(:,c_x),maxLatTraj_L2(:,c_y),maxLatTraj_L2(:,c_z));
    
    %%% Plotting maxNearZL1
    figure; hold all
    title('Max z near Europa - L1')
    PlotBoi3_CR3Bn(23)
    axis equal
    plot3(maxZTraj_L1(:,c_x),maxZTraj_L1(:,c_y),maxZTraj_L1(:,c_z));
    
    %%% Plotting maxNearZL2
    figure; hold all
    title('Max z near Europa - L2')
    PlotBoi3_CR3Bn(23)
    axis equal
    plot3(maxZTraj_L2(:,c_x),maxZTraj_L2(:,c_y),maxZTraj_L2(:,c_z));
    
    
    
    
    %%% Looking at z=0 section of these orbits
    %%% Plotting maxLatL1
    figure; hold all
    title('Max lat - L1')
    PlotBoi3_CR3Bn(23)
    axis equal
    plot3(maxLatTraj_L1(maxLatTraj_L1(:,c_isxyCross)==1,c_x),maxLatTraj_L1(maxLatTraj_L1(:,c_isxyCross)==1,c_y),maxLatTraj_L1(maxLatTraj_L1(:,c_isxyCross)==1,c_z),'.');
    
     %%% Plotting maxNearZL1
    figure; hold all
    title('Max z near Europa - L1')
    PlotBoi3_CR3Bn(23)
    axis equal
    plot3(maxZTraj_L1(maxZTraj_L1(:,c_isxyCross)==1,c_x),maxZTraj_L1(maxZTraj_L1(:,c_isxyCross)==1,c_y),maxZTraj_L1(maxZTraj_L1(:,c_isxyCross)==1,c_z),'.');
    
    
    
    
    
    %%% Looking more at maxZtraj_L1
    %%% Important
%     t = maxZTraj_L1(19261, c_t) - maxZTraj_L1(10000,c_t); % 1.540075234322559e+02
    t = 1.540075234322559e+02;
    X0 = maxZTraj_L1(10000,c_x:c_zd)';
%     X0 = maxLatTraj_L1(10000,c_x:c_zd)';
    tol = 1e-13;
    options = odeset('RelTol',tol,'AbsTol',tol);
    % [T_test, X_test] = ode45(@Int_CR3Bn, [0 t*8], X0, options, prms);
    % figure; hold all
    % PlotBoi3_CR3Bn(23)
    % axis equal
    % plot3(X_test(:,1),X_test(:,2),X_test(:,3))

%     options_event = odeset('Event',@event_yEqualsZero_nonTerminal, 'RelTol',tol,'AbsTol',tol);
%     options_event = odeset('Event',@event_xEqualsSecondary_nonTerminal, 'RelTol',tol,'AbsTol',tol);
    options_event = odeset('Event',@event_zEqualsZero_nonTerminal, 'RelTol',tol,'AbsTol',tol);
    % [T_test, X_test, time_event, X_event, index_event] = ode45(@Int_CR3Bn, [0 t*20], X0, options_event, prms);
    [T_test, X_test, time_event, X_event, index_event] = ode45(@Int_CR3Bn, [0 t*49], X0, options_event, prms); % t*49 max for ode45 for maxZ_L1, t*3.5 for maxLat_L1
%     [T_test, X_test, time_event, X_event, index_event] = ode113(@Int_CR3Bn, [0 t*38], X0, options_event, prms); % t*38 max for ode113 for maxZ_L1, t*10 for maxLat_L1
%     [T_test, X_test, time_event, X_event, index_event] = ode45(@Int_CR3Bn, [0 t*3.5], X0, options_event, prms); % t*49 max for ode45 for maxZ_L1, t*3.5 for maxLat_L1
%     [T_test, X_test, time_event, X_event, index_event] = ode113(@Int_CR3Bn, [0 t*10], X0, options_event, prms); % t*38 max for ode113 for maxZ_L1, t*10 for maxLat_L1
    % [T_test, X_test, time_event, X_event, index_event] = ode45(@Int_CR3Bn, [0 t*2], X0, options_event, prms);
    % [T_test, X_test, time_event, X_event, index_event] = ode45(@Int_CR3Bn, [t*6 0], X0, options_event, prms);
    figure; hold all
    PlotBoi3_CR3Bn(23)
    axis equal
    plot3(X_test(:,1),X_test(:,2),X_test(:,3))

    figure; hold all
    PlotBoi2('1 thing','other thing',23,'LaTex')
    axis equal
    plot(X_event(:,1),X_event(:,2),'.')
    
%     [X_SCI] = X_BaCR2SCI(X_test, T_test, prms);
%     a_e_i_raan_w_ta = zeros(size(X_SCI,1),6);
%     ecosw = zeros(size(X_SCI,1),1);
%     esinw = zeros(size(X_SCI,1),1);
%     for kk = 1:size(X_SCI,1)
%         [a,e,i,raan,w,ta] = ECI2OE(X_SCI(kk,1:3), X_SCI(kk,4:6),prms.u);
% 
%         a_e_i_raan_w_ta(kk,:) = [a,e,i,raan,w,ta];
%         ecosw(kk) = e * cos(w);
%         esinw(kk) = e * sin(w);
%     end
%     figure; hold all
%     plot(a_e_i_raan_w_ta(:,1), a_e_i_raan_w_ta(:,2))





    indicies=1:30;
    figure; hold all
    PlotBoi3_CR3Bn(23)
    plot3(maxZTraj_L1(indicies,c_x), maxZTraj_L1(indicies,c_y), maxZTraj_L1(indicies,c_z))
    plot3(newTraj1(indicies,c_x), newTraj1(indicies,c_y), newTraj1(indicies,c_z))
    plot3(newTraj3(indicies,c_x), newTraj3(indicies,c_y), newTraj3(indicies,c_z))
    plot3(newTraj4(indicies,c_x), newTraj4(indicies,c_y), newTraj4(indicies,c_z))
    plot3(newTraj5(indicies,c_x), newTraj5(indicies,c_y), newTraj5(indicies,c_z))
    plot3(newTraj6(indicies,c_x), newTraj6(indicies,c_y), newTraj6(indicies,c_z))
    
    
    
    %%% Looking more at maxZtraj_L2
    t = 145;
    X0 = maxZTraj_L2(25,c_x:c_zd)';
    tol = 1e-13;
    options = odeset('RelTol',tol,'AbsTol',tol);
%     options_event = odeset('Event',@event_yEqualsZero_nonTerminal, 'RelTol',tol,'AbsTol',tol);
%     options_event = odeset('Event',@event_xEqualsSecondary_nonTerminal, 'RelTol',tol,'AbsTol',tol);
    options_event = odeset('Event',@event_zEqualsZero_nonTerminal, 'RelTol',tol,'AbsTol',tol);
    [T_test_L2, X_test_L2, time_event_L2, X_event_L2, index_event_L2] = ode113(@Int_CR3Bn, [0 t], X0, options_event, prms); % 
    figure; hold all
    PlotBoi3_CR3Bn(23)
    axis equal
    plot3(X_test_L2(:,1),X_test_L2(:,2),X_test_L2(:,3))

    figure; hold all
    PlotBoi2('1 thing','other thing',23,'LaTex')
    axis equal
    plot(X_event_L2(:,1),X_event_L2(:,2),'.')
    
    
    
    %%% Looking at high traj at the lower bound of zd vs R at xy plane
    %%% crossing
    traj = data_high_1pi_xyCrossings(data_high_1pi_xyCrossings(:,c_trajID)==44254,:);
    t = 100*pi;
    X0 = traj(1,c_x:c_zd)';
%     X0 = maxLatTraj_L1(10000,c_x:c_zd)';
    tol = 1e-13;
    options_event = odeset('Event',@event_zEqualsZero_nonTerminal, 'RelTol',tol,'AbsTol',tol);
    % [T_test, X_test, time_event, X_event, index_event] = ode45(@Int_CR3Bn, [0 t*20], X0, options_event, prms);
    [T_test_high, X_test_high, time_event_high, X_event_high, index_event_high] = ode45(@Int_CR3Bn, [0 t], X0, options_event, prms); % t*49 max for ode45 for maxZ_L1, t*3.5 for maxLat_L1
    figure; hold all
    PlotBoi3_CR3Bn(23)
    axis equal
    plot3(X_test_high(:,1),X_test_high(:,2),X_test_high(:,3),'r')
   
    figure
    plot(rowNorm(X_event_high(:,1:3)-[1-prms.u,0,0]), X_event_high(:,6),'r.')
    
    
%     tempTraj = maxLatTraj_L1(5650:5890,:);
%     n_bins = 100;
%     vBins = linspace(min(abs(tempTraj(:,c_zd))),max(abs(tempTraj(:,c_zd))),n_bins+1);
%     colorGradient = colorScale([colors.blue2;colors.red2],n_bins);
%     bin_values = zeros(size(tempTraj,1),1);
%     for kk = 1:size(tempTraj,1)
%         bin_values(kk) = discretize(abs(tempTraj(kk,c_zd)),vBins);
%     end
%     binData = cell(n_bins,1);
%     for kk = 1:n_bins
%         binData{kk} = tempTraj(bin_values==kk,:);
%     end
%     figure; hold all
%     xlim([0 0.02])
%     for kk = 1:n_bins
%         plot(rowNorm(binData{kk}(:,c_x:c_y)-[1-prms.u,0]), binData{kk}(:,c_z),'.','color',colorGradient(kk,:),'markersize',15)
%     end
%     grid on
%     title('50, max lat L1 traj, ind 5650:5890')
%     PlotBoi2('$\sqrt{(x_n-1+\mu)^2 + y_n^2}$','$z_n$',23,'LaTex')
    

    %%% z vs Rxy for segment of trajectory - amplitude crunching
    tempTraj = maxLatTraj_L1(5650:5890,:);
    figure; hold all
    xlim([0 0.02])
    plot(rowNorm(tempTraj(:,c_x:c_y)-[1-prms.u,0]),   abs(tempTraj(:,c_z)),        'linewidth',1.5)
    p1 = plot(rowNorm(tempTraj(1,c_x:c_y)-[1-prms.u,0]),   tempTraj(1,c_z),   'bo', 'linewidth',1.5, 'markersize', 11);
    p2 = plot(rowNorm(tempTraj(end,c_x:c_y)-[1-prms.u,0]), tempTraj(end,c_z), 'bx', 'linewidth',1.5, 'markersize', 11);
    legend([p1 p2],'Beginning','End')
    grid on
    title('50, max lat L1 traj, ind 5650:5890')
    PlotBoi2('$\sqrt{(x_n-1+\mu)^2 + y_n^2}$','$z_n$',23,'LaTex')
    
    figure; hold all
    PlotBoi3_CR3Bn(23)
    axis equal
    plot3(tempTraj(:,c_x),tempTraj(:,c_y),tempTraj(:,c_z))
end

%%% Looking at x-y plane crossings
if 1+1==1
    data_L1_1yr_NearXyCrossingsOnly = data_L1_1yr_xyCrossingsOnly((rowNorm(data_L1_1yr_xyCrossingsOnly(:,c_x:c_z) - [1-prms.u,0,0]) < max(Rs_high_1pi)), :);
    data_L2_1yr_NearXyCrossingsOnly = data_L2_1yr_xyCrossingsOnly((rowNorm(data_L2_1yr_xyCrossingsOnly(:,c_x:c_z) - [1-prms.u,0,0]) < max(Rs_high_1pi)), :);
    
    figure; hold all
    plot(Rs_high_1pi(indicies_high_xyCrossingEvents), rowNorm(data_high_1pi_xyCrossingsOnly(:,c_hx:c_hy)), '.', 'color', color_high_1pi);
    plot(rowNorm(data_L2_1yr_NearXyCrossingsOnly(:,c_x:c_z)-[1-prms.u,0,0]), rowNorm(data_L2_1yr_NearXyCrossingsOnly(:,c_hx:c_hy)), '.', 'color', color_L2_1yr);

    figure; hold all
    plot(Rs_high_1pi(indicies_high_xyCrossingEvents), data_high_1pi_xyCrossingsOnly(:,c_i), '.', 'color', color_high_1pi);
    plot(rowNorm(data_L1_1yr_NearXyCrossingsOnly(:,c_x:c_z)-[1-prms.u,0,0]), data_L1_1yr_NearXyCrossingsOnly(:,c_i), 'x', 'color', color_L1_1yr);
    plot(rowNorm(data_L2_1yr_NearXyCrossingsOnly(:,c_x:c_z)-[1-prms.u,0,0]), data_L2_1yr_NearXyCrossingsOnly(:,c_i), 'x', 'color', color_L2_1yr);
    
    
    figure; hold all
    PlotBoi2('$R_n$','Inclination, $^\circ$',23,'LaTex')
    p_high_1pi = plot(Rs_high_1pi(indicies_high_xyCrossingEvents), data_high_1pi_xyCrossingsOnly(:,c_i), '.', 'color', color_high_1pi);
    p_L1_1yr = plot(rowNorm(data_L1_1yr_NearXyCrossingsOnly(:,c_x:c_z)-[1-prms.u,0,0]), data_L1_1yr_NearXyCrossingsOnly(:,c_i), '.', 'color', color_L1_1yr);
    p_L2_1yr = plot(rowNorm(data_L2_1yr_NearXyCrossingsOnly(:,c_x:c_z)-[1-prms.u,0,0]), data_L2_1yr_NearXyCrossingsOnly(:,c_i), '.', 'color', color_L2_1yr);
    p_surf   = plot([secondary.R_n secondary.R_n],[0 140],'k', 'linewidth',1.5);
    [legh,objh] = legend([p_high_1pi, p_L1_1yr, p_L2_1yr, p_surf], label_high_1pi, label_L1_1yr, label_L2_1yr, 'Europa Surface');
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);
    
    
    figure; hold all
    PlotBoi2('Longitude, $^\circ$','Inclination, $^\circ$',23,'LaTex')
    xlim([-180 180])
    p_high_1pi = plot(data_high_1pi_xyCrossingsOnly(:,c_lon), data_high_1pi_xyCrossingsOnly(:,c_i), '.', 'color', color_high_1pi);
    p_L1_1yr = plot(data_L1_1yr_NearXyCrossingsOnly(:,c_lon), data_L1_1yr_NearXyCrossingsOnly(:,c_i), '.', 'color', color_L1_1yr);
    p_L2_1yr = plot(data_L2_1yr_NearXyCrossingsOnly(:,c_lon), data_L2_1yr_NearXyCrossingsOnly(:,c_i), '.', 'color', color_L2_1yr);
    [legh,objh] = legend([p_high_1pi, p_L1_1yr, p_L2_1yr], label_high_1pi, label_L1_1yr, label_L2_1yr);
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);
    
    
    figure; hold all
    PlotBoi2('$R_n$','$\dot{z}_n$',23,'LaTex')
    ylim([-1 1].*0.15)
    xlim([0 max(Rs_high_1pi)])
    p_high_1pi = plot(Rs_high_1pi(indicies_high_xyCrossingEvents), data_high_1pi_xyCrossingsOnly(:,c_zd), '.', 'color', color_high_1pi);
    p_L2_1yr = plot(rowNorm(data_L2_1yr_NearXyCrossingsOnly(:,c_x:c_z)-[1-prms.u,0,0]), data_L2_1yr_NearXyCrossingsOnly(:,c_zd), '.', 'color', color_L2_1yr);
    p_surf   = plot([secondary.R_n secondary.R_n],[-1 1].*.3,'k', 'linewidth',1.5);
    [legh,objh] = legend([p_high_1pi, p_L2_1yr, p_surf], label_high_1pi, label_L2_1yr, 'Europa Surface');
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);
    
    plot(rowNorm(maxLatTraj_L2(maxLatTraj_L2(:,c_isxyCross)==1,c_x:c_z)-[1-prms.u,0,0]), maxLatTraj_L2(maxLatTraj_L2(:,c_isxyCross)==1,c_zd),'k.')
   

    
    figure; hold all
    PlotBoi2('$R_n$','$\dot{z}_n$',23,'LaTex')
    ylim([-1 1].*0.15)
    xlim([0 max(Rs_high_1pi)])
    p_high_1pi = plot(Rs_high_1pi(indicies_high_xyCrossingEvents), data_high_1pi_xyCrossingsOnly(:,c_zd), '.', 'color', color_high_1pi);
    p_L1_1yr = plot(rowNorm(data_L1_1yr_NearXyCrossingsOnly(:,c_x:c_z)-[1-prms.u,0,0]), data_L1_1yr_NearXyCrossingsOnly(:,c_zd), '.', 'color', color_L1_1yr);
    p_surf   = plot([secondary.R_n secondary.R_n],[-1 1].*.3,'k', 'linewidth',1.5);
    [legh,objh] = legend([p_high_1pi, p_L1_1yr, p_surf], label_high_1pi, label_L1_1yr, 'Europa Surface');
    lineh = findobj(objh,'type','line');
    set(lineh,'linestyle','-','linewidth',3);
    
    plot(rowNorm(maxZTraj_L1(maxZTraj_L1(:,c_isxyCross)==1,c_x:c_z)-[1-prms.u,0,0]), maxZTraj_L1(maxZTraj_L1(:,c_isxyCross)==1,c_zd),'k.')
%     plot(rowNorm(maxLatTraj_L1(maxLatTraj_L1(:,c_isxyCross)==1,c_x:c_z)-[1-prms.u,0,0]), maxLatTraj_L1(maxLatTraj_L1(:,c_isxyCross)==1,c_zd),'k.')

    
%     figure; hold all
%     PlotBoi2('$R_n$','$\dot{z}_n$/$R_n$', 23, 'LaTex')
%     p_high_1pi = plot(Rs_high_1pi(indicies_high_xyCrossingEvents), data_high_1pi_xyCrossingsOnly(:,c_zd)./Rs_high_1pi(indicies_high_xyCrossingEvents), '.', 'color', color_high_1pi);
%     p_L2_1yr = plot(rowNorm(data_L2_1yr_NearXyCrossingsOnly(:,c_x:c_z)-[1-prms.u,0,0]), data_L2_1yr_NearXyCrossingsOnly(:,c_zd)./rowNorm(data_L2_1yr_NearXyCrossingsOnly(:,c_x:c_z)-[1-prms.u,0,0]), '.', 'color', color_L2_1yr);
%     p_L1_1yr = plot(rowNorm(data_L1_1yr_NearXyCrossingsOnly(:,c_x:c_z)-[1-prms.u,0,0]), data_L1_1yr_NearXyCrossingsOnly(:,c_zd)./rowNorm(data_L1_1yr_NearXyCrossingsOnly(:,c_x:c_z)-[1-prms.u,0,0]), '.', 'color', color_L1_1yr);
%     p_surf   = plot([secondary.R_n secondary.R_n],[-1 1].*.3,'k', 'linewidth',1.5);
%     [legh,objh] = legend([p_high_1pi, p_L1_1yr, p_surf], label_high_1pi, label_L1_1yr, 'Europa Surface');
%     lineh = findobj(objh,'type','line');
%     set(lineh,'linestyle','-','linewidth',3);
end




% testing = [];
% for kk = testIndicies(2):testIndicies(end)
%     if (maxZTraj_L1(kk,c_lon) < 0) && (maxZTraj_L1(kk-1,c_lon) > 0)
%         testing = [testing;NaN(1,27); maxZTraj_L1(kk,:)];
%     else
%         testing = [testing; maxZTraj_L1(kk,:)];
%     end
% end




% figure; hold all
% plot(maxZTraj_L1(testIndicies,c_lon), maxZTraj_L1(testIndicies,c_lat))
% 
% a_e_i_raan_w_ta = zeros(size(maxZTraj_L1,1),6);
% ecosw = zeros(size(maxZTraj_L1,1),1);
% esinw = zeros(size(maxZTraj_L1,1),1);
% for kk = 1:size(maxZTraj_L1,1)
%     [a,e,i,raan,w,ta] = ECI2OE(maxZTraj_L1(kk,c_x:c_z) - [1-prms.u,0,0],maxZTraj_L1(kk,c_xd:c_zd),prms.u);
%     
%     a_e_i_raan_w_ta(kk,:) = [a,e,i,raan,w,ta];
%     ecosw(kk) = e * cos(w);
%     esinw(kk) = e * sin(w);
% end
% 
% figure; plot(ecosw(20000:25000), esinw(20000:25000))


% X0 = maxZTraj_L1(1,c_x:c_zd)';
% tol = 1e-13;
% options = odeset('RelTol',tol,'AbsTol',tol);
% [T_test, X_test] = ode45(@Int_CR3Bn, [0 maxZTraj_L1(end,c_t)], X0, options, prms);
% figure; hold all
% PlotBoi3_CR3Bn(23)
% axis equal
% plot3(X_test(:,1),X_test(:,2),X_test(:,3))




% % 0.002325733869766
% x0 = 1-prms.u + 0.009;
% JC_X0 = getJacobiConstant_ZH([x0, 0, 0, 0, 0, 0],prms);
% yd0 = sqrt(JC_X0 - JC_data);
% X0 = [x0, 0, 0, 0, yd0, 0]';
% tol = 1e-13;
% options = odeset('RelTol',tol,'AbsTol',tol);
% [T_test, X_test] = ode113(@Int_CR3Bn, [0 1.4], X0, options, prms);
% figure; hold all
% PlotBoi3_CR3Bn(23)
% axis equal
% plot3(X_test(:,1),X_test(:,2),X_test(:,3))










% [X_SCI] = X_BaCR2SCI(X_test, T_test, prms);
% figure; plot3(X_SCI(:,1),X_SCI(:,2),X_SCI(:,3))

% %%% From Max z near Europa - L1, (45000:46170)
% X0_test = [1.001749837771900;
%              0.001002142717300;
%              0.000627496093000;
%              -0.006573771996300;
%              0.134193761362700;
%              0.044354471210300];
% tol = 1e-13;
% options = odeset('RelTol',tol,'AbsTol',tol);
% 
% 
% [T_test, X_test] = ode113(@Int_CR3Bn, [0 1*(86400*365.25)/tNorm], X0_test, options, prms);
% 
% figure; hold all
% title('Test, Max z near Europa - L1')
% PlotBoi3_CR3Bn(23)
% axis equal
% plot3(X_test(1:4000,1),X_test(1:4000,2),X_test(1:4000,3))
% 
% PO_test = [X0_test; 13.564585838486098];




% find(abs(data_L2_1yr(:, c_z)) == max(abs(data_L2_1yr((data_L2_1yr(:,c_R) < 0.02), c_z))))
% find(abs(data_L1_1yr(:, c_lat)) == max(abs(data_L1_1yr(:, c_lat))))
%%% From 100 L1, highest near-z trajID: 5334
%%% From 100 L2, highest near-z trajID: 1234
%%% From 100 L1, highest lat    trajID: 4955
%%% From 100 L2, highest lat    trajID: 3287
% 
% 
% figure; hold all
% axis equal
% PlotBoi3_CR3Bn(23)
% % indicies = data_L1_1yr(:,c_trajID)==4955;
% % plot3(data_L1_1yr(indicies,c_x),data_L1_1yr(indicies,c_y),data_L1_1yr(indicies,c_z));
% % title('100mps L1, max lat, trajID 4955')
% indicies = data_L2_1yr(:,c_trajID)==3287;
% plot3(data_L2_1yr(indicies,c_x),data_L2_1yr(indicies,c_y),data_L2_1yr(indicies,c_z));
% title('100mps L2, max lat, trajID 3287')





