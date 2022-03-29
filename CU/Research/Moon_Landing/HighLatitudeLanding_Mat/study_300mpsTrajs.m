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
load([dataPath,'lowTrajs_300mps_neckFull_L1_250km_22v0s_8pi.mat']); % data_L1_8pi
load([dataPath,'lowTrajs_300mps_neckFull_L2_250km_22v0s_8pi.mat']); % data_L2_8pi

data_title = '300 m/s';

% --------------------------
% Create near data (inside necks)
% --------------------------
% %%% Set max distance from secondary ... (max distance to L1 or L2)
% Rmax = max([L123(2,1)-(1-secondary.MR), (1-secondary.MR)-L123(1,1)]);
% 
% %%% Grab the near-Europa data
% nearData_L1_8pi = data_L1_8pi( (data_L1_8pi(:,c_R) <= Rmax) ,:);
% nearData_L2_8pi = data_L2_8pi( (data_L2_8pi(:,c_R) <= Rmax) ,:);
% 
% %%% Save the data
% save([dataPath,'lowTrajs_300mps_neckFull_L1_250km_22v0s_8pi_nearEur','.mat'], 'nearData_L1_8pi','-v7.3')
% save([dataPath,'lowTrajs_300mps_neckFull_L2_250km_22v0s_8pi_nearEur','.mat'], 'nearData_L2_8pi','-v7.3')

% --------------------------
% Load .mat files of near data (inside necks)
% --------------------------
load([dataPath,'lowTrajs_300mps_neckFull_L1_250km_22v0s_8pi_nearEur.mat']); % nearData_L1_8pi
load([dataPath,'lowTrajs_300mps_neckFull_L2_250km_22v0s_8pi_nearEur.mat']); % nearData_L2_8pi

% --------------------------
% Load .mat files of data with x-y crossings
% --------------------------
% neardata_L1_8pi_xyCrossingsOnly = nearData_L1_8pi((nearData_L1_8pi(:,c_z)==0), :);
% neardata_L2_8pi_xyCrossingsOnly = nearData_L1_8pi((nearData_L2_8pi(:,c_z)==0), :);
% 
% save([dataPath,'lowTrajs_300mps_neckFull_L1_250km_22v0s_8pi_nearEur_xyCrossingsOnly','.mat'], 'neardata_L1_8pi_xyCrossingsOnly','-v7.3')
% save([dataPath,'lowTrajs_300mps_neckFull_L2_250km_22v0s_8pi_nearEur_xyCrossingsOnly','.mat'], 'neardata_L2_8pi_xyCrossingsOnly','-v7.3')

load([dataPath,'lowTrajs_300mps_neckFull_L1_250km_22v0s_8pi_nearEur_xyCrossingsOnly.mat']); % neardata_L1_8pi_xyCrossingsOnly
load([dataPath,'lowTrajs_300mps_neckFull_L2_250km_22v0s_8pi_nearEur_xyCrossingsOnly.mat']); % neardata_L2_8pi_xyCrossingsOnly

% --------------------------
% labels and colors
% --------------------------
label_L1_8pi        = 'L1, 8pi';
label_L2_8pi        = 'L2, 8pi';

color_L1_8pi        = colors.blue2;
color_L2_8pi        = colors.red2;


% -------------------------------------------------
%%% Plotting datasets
% -------------------------------------------------
if 1+1==1
%     color_L1_8pi        = colors.purp;
%     color_L2_8pi        = colors.cyan;
% %     color_high_1pi      = colors.red;
%     
%     figure; hold all
%     axis equal
%     xlim([L123(1,1)-secondary.R_n, L123(2,1)+secondary.R_n])
%     PlotBoi3_CR3Bn(26)
%     
% %     p_high_1pi = plot3(data_high_1pi(:,c_x),data_high_1pi(:,c_y),data_high_1pi(:,c_z),'.','markersize',1, 'color', color_high_1pi);
%     p_L1_8pi   = plot3(data_L1_8pi(:,c_x),data_L1_8pi(:,c_y),data_L1_8pi(:,c_z), '.','markersize',1, 'color', color_L1_8pi);
%     p_L2_8pi   = plot3(data_L2_8pi(:,c_x),data_L2_8pi(:,c_y),data_L2_8pi(:,c_z), '.','markersize',1, 'color', color_L2_8pi);
%     
% %     [legh,objh] = legend([p_high_1pi, p_L1_8pi, p_L2_8pi], label_high_1pi  ,label_L1_8pi, label_L2_8pi,'FontSize',14);
%     [legh,objh] = legend([ p_L1_8pi, p_L2_8pi]  ,label_L1_8pi, label_L2_8pi,'FontSize',14);
%     lineh = findobj(objh,'type','line');
%     set(lineh,'linestyle','-','linewidth',3);
end

% -------------------------------------------------
%%% Plotting lat/lon projections
% -------------------------------------------------
if 1+1==1
    figure; hold all
    PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$',26,'LaTex')
    xlim([-180 180])
    ylim([-90 90])
    
    p_L1_8pi   = plot(data_L1_8pi(:,c_lon),data_L1_8pi(:,c_lat), '.','markersize',1, 'color', color_L1_8pi);
    p_L2_8pi   = plot(data_L2_8pi(:,c_lon),data_L2_8pi(:,c_lat), '.','markersize',1, 'color', color_L2_8pi);
    
    [legh,objh] = legend([p_L1_8pi, p_L2_8pi]  ,label_L1_8pi, label_L2_8pi,'FontSize',14);
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
    plotCR3BP_YZNeck( JC, secondary.MR , 1, 0, prms, color_L1_8pi, 2);
    plotCR3BP_YZNeck( JC, secondary.MR , 2, 0, prms, color_L2_8pi, 2);
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
    L1_ICs = data_L1_8pi(data_L1_8pi(:,c_t)==0, :);
    L2_ICs = data_L2_8pi(data_L2_8pi(:,c_t)==0, :);
    
    figure; hold all
    JC = L2FlyoverVelocity_2_JC(300, secondary.MR, L123(2,:), vNorm, 1);
    prms.u = secondary.MR;
    prms.n = 1;
    prms.R2_n = secondary.R_n;
    plotCR3BP_YZNeck( JC, secondary.MR , 1, 0, prms, color_L1_8pi, 2);
    plotCR3BP_YZNeck( JC, secondary.MR , 2, 0, prms, color_L2_8pi, 2);
    plotBodyTexture3(secondary.R_n,[1-secondary.MR,0,0],secondary.img)
    
    plot3(L1_ICs(:,c_x), L1_ICs(:,c_y), L1_ICs(:,c_z),'.','color',colors.black)
    plot3(L2_ICs(:,c_x), L2_ICs(:,c_y), L2_ICs(:,c_z),'.','color',colors.black)
    
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
p_L1_8pi   = plot(neardata_L1_8pi_xyCrossingsOnly(:,c_R), neardata_L1_8pi_xyCrossingsOnly(:,c_i),'.','color',color_L1_8pi);
p_L2_8pi = plot(neardata_L2_8pi_xyCrossingsOnly(:,c_R), neardata_L2_8pi_xyCrossingsOnly(:,c_i),'.','color',color_L2_8pi);
p_surf   = plot([secondary.R_n secondary.R_n],[0 140],'k', 'linewidth',1.5);
[legh,objh] = legend([p_L1_8pi, p_L2_8pi, p_surf], label_L1_8pi, label_L2_8pi, 'Europa Surface','FontSize',14);
lineh = findobj(objh,'type','line');
set(lineh,'linestyle','-','linewidth',3);
ylim([0 140])
xlim([0 0.016])
end


% -------------------------------------------------
%%% Looking at trajectories with maximum zdot at x-y crossing
% -------------------------------------------------
if 1+1==1
%%% Poincaré section results
figure; hold all
p_L1_8pi = plot(neardata_L1_8pi_xyCrossingsOnly(:,c_R), abs(neardata_L1_8pi_xyCrossingsOnly(:,c_zd)),'.','color',color_L1_8pi);
p_L2_8pi = plot(neardata_L2_8pi_xyCrossingsOnly(:,c_R), abs(neardata_L2_8pi_xyCrossingsOnly(:,c_zd)),'.','color',color_L2_8pi);
p_surf   = plot([secondary.R_n secondary.R_n],[0 3],'k', 'linewidth',1.5);
[legh,objh] = legend([p_L1_8pi, p_L2_8pi, p_surf], label_L1_8pi, label_L2_8pi, 'Europa Surface','FontSize',14);
lineh = findobj(objh,'type','line');
set(lineh,'linestyle','-','linewidth',3);
ylim([0 0.15])
xlim([0 0.016])
PlotBoi2('$r_2$','$|\dot{z}|$',26,'LaTex')


end


% -------------------------------------------------
%%% Looking at certain data
% -------------------------------------------------
if 1+1==1

    test_zd = neardata_L1_8pi_xyCrossingsOnly((neardata_L1_8pi_xyCrossingsOnly(:,c_R)>=(0.0129)), :); 
    test_zd = test_zd((test_zd(:,c_R) <= (0.01291)), :);
    test_zd = test_zd((abs(test_zd(:,c_zd)) >= (0.0332)), :);
    
    size(test_zd)
    test_zd(1,c_trajID)
    
    testTraj_zd = nearData_L1_8pi((nearData_L1_8pi(:,c_trajID) == test_zd(1,c_trajID)), :);
    
    figure; hold all
    plot3(testTraj_zd(:,c_x), testTraj_zd(:,c_y), testTraj_zd(:,c_z),'b')
    PlotBoi3_CR3Bn(26)
end


% % ==================================================================
% %%% Comparing z and zd within the neck
% % ==================================================================
% if 1+1==1
%     figure; hold all
%     PlotBoi2('$r_2$','$z$',26,'LaTex')
%     plot(neardata_L2_8pi_xyCrossings_dense(:,c_R),neardata_L2_8pi_xyCrossings_dense(:,c_z), '.', 'markersize', 10,'color',colors.grey)
%     xlim([0 L123(2,1)-1+prms.u])
%     
%     figure; hold all
%     PlotBoi2('$r_2$','$\dot{z}$',26,'LaTex')
%     plot(neardata_L2_8pi_xyCrossings_dense(:,c_R),neardata_L2_8pi_xyCrossings_dense(:,c_zd), '.', 'markersize', 10,'color',colors.grey)
%     xlim([0 L123(2,1)-1+prms.u])
%     
% end
% 
% 
% % ========================================================================
% %%% Testing
% % ========================================================================
% % warning('Need to address the fact that I haven''t integrated the high trajs for very long')
% % 
% % 
% % data = data_high_1pi_xyCrossingsOnly(Rs_high_1pi_xyCrossingsOnly>5.5e-3,:);
% % data = data((rowNorm(data(:,c_x:c_z) - [1-prms.u,0,0]))<5.6e-3,:);
% % data = data(data(:,c_i)<36,:);
% % data = data(data(:,c_i)>35.5,:);
% % 
% % high_minInc_trajID = 48846;
% % high_minInc_data = data_high_1pi_xyCrossingsOnly(data_high_1pi_xyCrossingsOnly(:,c_trajID)==high_minInc_trajID,:);
% % plot(rowNorm(high_minInc_data(:,c_x:c_z)-[1-secondary.MR,0,0]), high_minInc_data(:,c_i),'.','color',colors.black);
% % 
% % 
% % high_minInc_traj = data_high_1pi_xyCrossings(data_high_1pi_xyCrossings(:,c_trajID)==high_minInc_trajID,:);
% % 
% % 
% % 
% % figure; hold all
% % plot3(high_minInc_traj(:,c_x),high_minInc_traj(:,c_y),high_minInc_traj(:,c_z))
% % 
% % 
% % options_event = odeset('Event',@event_zEqualsZero_nonTerminal, 'RelTol',tol,'AbsTol',tol);
% % X0 = high_minInc_traj(1,c_x:c_zd)';
% % [T, X, time_event, X_event, index_event] = ode113(@Int_CR3Bn, [0 20*pi], X0, options_event, prms);
% % figure; hold all
% % plot3(X(:,1),X(:,2),X(:,3),'r')
% % axis normal
% % PlotBoi3_CR3Bn(26)
% % 
% % 
% %     
% % [X_SCI_event] = X_BaCR2SCI(X_event, time_event, prms);
% % a_e_i_raan_w_ta = zeros(size(X_SCI_event,1),6);
% % 
% % for kk = 1:size(X_SCI_event,1)
% %     [a,e,i,raan,w,ta] = ECI2OE(X_SCI_event(kk,1:3), X_SCI_event(kk,4:6),prms.u);
% % 
% %     a_e_i_raan_w_ta(kk,:) = [a,e,i,raan,w,ta];
% % end
% % 
% % plot(rowNorm(X_event(:,1:3) - [1-prms.u,0,0]), a_e_i_raan_w_ta(:,3)*180/pi,'m.')
% 
% if 1+1==1
%     tol     = 1e-13;
%     options_event = odeset('Event',@event_zEqualsZero_nonTerminal, 'RelTol',tol,'AbsTol',tol);
%     X0 = L1_maxInc_traj(10,c_x:c_zd)';
%     [T, X, time_event, X_event, index_event] = ode113(@Int_CR3Bn, [0 646], X0, options_event, prms);
%     figure; hold all
%     plot3(X(:,1),X(:,2),X(:,3),'r')
%     axis normal
%     PlotBoi3_CR3Bn(26)
% 
%     figure; hold all
%     plot3(X_event(:,1),X_event(:,2),X_event(:,3),'k.','markersize',10)
%     PlotBoi3_CR3Bn(26)
%     axis normal
% end
% 
% 
% % Just looking at a trajectory from L2 dense data
% if 1+1==1
%     newData = neardata_L2_8pi_xyCrossings_dense(neardata_L2_8pi_xyCrossings_dense(:,c_z)>0.00358,:);
%     newData = newData(newData(:,c_z)<0.0036,:);
%     newData = newData(newData(:,c_R)<0.0164, :);
%     newData = newData(newData(:,c_R)>0.0163, :);
%     newTrajID = newData(1,c_trajID);
%     newTraj =  neardata_L2_8pi_xyCrossings_dense(neardata_L2_8pi_xyCrossings_dense(:,c_trajID)==newTrajID,:);
%     
%     figure; hold all
%     PlotBoi3_CR3Bn(26)
%     axis equal
%     plot3(newTraj(:,c_x),newTraj(:,c_y),newTraj(:,c_z),'b')
%     plot3(newData(1,c_x),newData(1,c_y),newData(1,c_z),'o','markeredgecolor',colors.black,'markerfacecolor',colors.red)
% end
% 
% 
% 
% 
% 
% 
% 
% 
% 
% % Looking at inclination in rotating frame
% if 1+1==1
% %     Rs_high_1pi_xyCrossingsOnly
% %     data_high_1pi_xyCrossingsOnly
% %     neardata_L1_8pi_xyCrossingsOnly
% 
%     is_SCR_high_xyCrossingsOnly = NaN(size(Rs_high_1pi_xyCrossingsOnly));
%     is_SCR_L1_xyCrossingsOnly = NaN(size(Rs_L1_1yr_near_xyCrossingsOnly));
%     is_SCR_L2_xyCrossingsOnly = NaN(size(Rs_L2_1yr_near_xyCrossingsOnly));
%     
%     for kk = 1:length(is_SCR_high_xyCrossingsOnly)
%         [a,e,i,raan,w,ta] = ECI2OE(data_high_1pi_xyCrossingsOnly(kk,c_x:c_z)-[1-secondary.MR,0,0],data_high_1pi_xyCrossingsOnly(kk,c_xd:c_zd),secondary.MR);
%         is_SCR_high_xyCrossingsOnly(kk) = i*180/pi;
%     end
%     
%     for kk = 1:length(is_SCR_L1_xyCrossingsOnly)
%         [a,e,i,raan,w,ta] = ECI2OE(neardata_L1_8pi_xyCrossingsOnly(kk,c_x:c_z)-[1-secondary.MR,0,0],neardata_L1_8pi_xyCrossingsOnly(kk,c_xd:c_zd),secondary.MR);
%         is_SCR_L1_xyCrossingsOnly(kk) = i*180/pi;
%     end
%     
%     for kk = 1:length(is_SCR_L2_xyCrossingsOnly)
%         [a,e,i,raan,w,ta] = ECI2OE(neardata_L2_8pi_xyCrossingsOnly(kk,c_x:c_z)-[1-secondary.MR,0,0],neardata_L2_8pi_xyCrossingsOnly(kk,c_xd:c_zd),secondary.MR);
%         is_SCR_L2_xyCrossingsOnly(kk) = i*180/pi;
%     end
%     
%     figure; hold all
%     PlotBoi2('$r_2$','Inclination (SCR), $^\circ$',26,'LaTex')
%     p_high_1pi = plot(Rs_high_1pi_xyCrossingsOnly, is_SCR_high_xyCrossingsOnly,'.','color',color_high_1pi);
%     p_L1_8pi   = plot(Rs_L1_1yr_near_xyCrossingsOnly, is_SCR_L1_xyCrossingsOnly,'.','color',color_L1_8pi);
%     p_L2_8pi   = plot(Rs_L2_1yr_near_xyCrossingsOnly, is_SCR_L2_xyCrossingsOnly,'.','color',color_L2_8pi);
%     p_surf   = plot([secondary.R_n secondary.R_n],[0 180],'k', 'linewidth',1.5);
%     [legh,objh] = legend([p_high_1pi, p_L1_8pi, p_L2_8pi, p_surf], label_high_1pi  ,label_L1_8pi, label_L2_8pi, 'Europa Surface','FontSize',14);
%     lineh = findobj(objh,'type','line');
%     set(lineh,'linestyle','-','linewidth',3);
%     xlim([0 0.016])
%     
% end
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













    