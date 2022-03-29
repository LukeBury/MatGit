% ========================================================================
%%% Description
% ========================================================================
% For looking at abort-viable PO solutions that resulted from bifurcation
% study at Europa

% Created: 6/17/21
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

% -------------------------------------------------
%%% Load families that have abort solutions
% -------------------------------------------------
family_SHalo           = 'Jupiter_Europa.CR3BP.L2_SHalo.txt';
family_L2_L_1T_4P2_1P3 = 'Jupiter_Europa.CR3BP.L2_L_1T_4P2_1P3.txt'; % the end of this is the 2P2 bifurcation from LoPO_2P2_1T (so this is also LoPO_2P2_1T_2P2)
family_DRO             = 'Jupiter_Europa.CR3BP.DRO.txt'; 
family_DPO             = 'Jupiter_Europa.CR3BP.DPO.txt'; 
family_DPO_2P3         = 'Jupiter_Europa.CR3BP.DPO_2P3.txt'; 
family_DPO_2P4         = 'Jupiter_Europa.CR3BP.DPO_2P4.txt'; 
family_Hg2             = 'Jupiter_Europa.CR3BP.Hg2.txt'; 
family_Hg2_2P2         = 'Jupiter_Europa.CR3BP.Hg2_2P2.txt';
family_Hg2_2T          = 'Jupiter_Europa.CR3BP.Hg2_2T.txt'; 
family_Hg2_3T          = 'Jupiter_Europa.CR3BP.Hg2_3T.txt'; 
family_Se7             = 'Jupiter_Europa.CR3BP.Se7.txt'; 
family_Se7_6P2         = 'Jupiter_Europa.CR3BP.Se7_6P2.txt'; 
family_LoPO            = 'Jupiter_Europa.CR3BP.LoPO.txt'; 
family_LoPO_2P2_1T     = 'Jupiter_Europa.CR3BP.LoPO_2P2_1T.txt'; 
family_LoPO_2P4        = 'Jupiter_Europa.CR3BP.LoPO_2P4.txt'; 

% --------------------------
% Actually load families
% --------------------------
%%% Path from mbin to data
dataPathFromMBin = '/Data/InitialConditions/PO_Families/';

%%% PO data file
PO_data_SHalo           = dlmread([mbinPath, dataPathFromMBin, family_SHalo],           ',',1,0);
PO_data_L2_L_1T_4P2_1P3 = dlmread([mbinPath, dataPathFromMBin, family_L2_L_1T_4P2_1P3], ',',1,0);
PO_data_DRO             = dlmread([mbinPath, dataPathFromMBin, family_DRO],             ',',1,0);
PO_data_DPO             = dlmread([mbinPath, dataPathFromMBin, family_DPO],             ',',1,0);
PO_data_DPO_2P3         = dlmread([mbinPath, dataPathFromMBin, family_DPO_2P3],         ',',1,0);
PO_data_DPO_2P4         = dlmread([mbinPath, dataPathFromMBin, family_DPO_2P4],         ',',1,0);
PO_data_Hg2             = dlmread([mbinPath, dataPathFromMBin, family_Hg2],             ',',1,0);
PO_data_Hg2_2P2         = dlmread([mbinPath, dataPathFromMBin, family_Hg2_2P2],         ',',1,0);
PO_data_Hg2_2T          = dlmread([mbinPath, dataPathFromMBin, family_Hg2_2T],          ',',1,0);
PO_data_Hg2_3T          = dlmread([mbinPath, dataPathFromMBin, family_Hg2_3T],          ',',1,0);
PO_data_Se7             = dlmread([mbinPath, dataPathFromMBin, family_Se7],             ',',1,0);
PO_data_Se7_6P2         = dlmread([mbinPath, dataPathFromMBin, family_Se7_6P2],         ',',1,0);
PO_data_LoPO            = dlmread([mbinPath, dataPathFromMBin, family_LoPO],            ',',1,0);
PO_data_LoPO_2P2_1T     = dlmread([mbinPath, dataPathFromMBin, family_LoPO_2P2_1T],     ',',1,0);
PO_data_LoPO_2P4        = dlmread([mbinPath, dataPathFromMBin, family_LoPO_2P4],        ',',1,0);

%%% Grabbing specific POs
PO_SHalo__ns1           = PO_data_SHalo(493,:);
PO_L2_L_1T_4P2_1P3__ns1 = PO_data_L2_L_1T_4P2_1P3(1788,:);
PO_L2_L_1T_4P2_1P3__ns2 = PO_data_L2_L_1T_4P2_1P3(1858,:);
PO_DRO__s1              = PO_data_DRO(28,:);
PO_DPO__s1              = PO_data_DPO(245,:);
PO_DPO_2P3__ns1         = PO_data_DPO_2P3(42,:);
PO_DPO_2P4__ns1         = PO_data_DPO_2P4(106,:);
PO_Hg2__ns1             = PO_data_Hg2(927,:);
PO_Hg2_2P2__ns1         = PO_data_Hg2_2P2(253,:);
PO_Hg2_2T__ns1          = PO_data_Hg2_2T(39,:);
PO_Hg2_3T__ns1          = PO_data_Hg2_3T(8,:);
PO_Hg2_3T__ns2          = PO_data_Hg2_3T(1245,:);
PO_Se7__s1              = PO_data_Se7(305,:);
PO_Se7__ns1             = PO_data_Se7(468,:);
PO_Se7_6P2__ns1         = PO_data_Se7_6P2(119,:);
PO_LoPO__s1             = PO_data_LoPO(404,:);
PO_LoPO_2P2_1T__ns1     = PO_data_LoPO_2P2_1T(289,:);
PO_LoPO_2P4__ns1        = PO_data_LoPO_2P4(124,:);

%%% Column specifiers 
c_x0 = 1;   c_y0 = 2;   c_z0 = 3;
c_xd0 = 4;  x_yd0 = 5;  c_zd0 = 6;
c_Tp = 7;   c_JC = 8;   c_s1 = 9;   c_s2 = 10;
c_alpha = 11;   c_beta = 12;    c_impactFlag = 13;
c_error = 14;

return
% ========================================================================
%%% Integrate and plot the POs
% ========================================================================
% -------------------------------------------------
%%% Set up parameters
% -------------------------------------------------
%%% Set primary and secondary bodies
[primary, secondary] = assignPrimaryAndSecondary_CR3BP('Jupiter_Europa.CR3BP', bodies);

%%% Normalizing constants
[rNorm, tNorm, vNorm] = cr3bp_norms(primary, secondary, bodies.constants.G);

%%% prms for integration
prms.u    = secondary.MR;
prms.R2 = secondary.R_n;
prms.n    = 1;
rLPs_n = EquilibriumPoints(prms.u, prms.n);


%%% Integration options
tol = 1e-13;
options = odeset('RelTol',tol,'AbsTol',tol);
options_apsis = odeset('Events', @event_Apsis_CR3BP, 'RelTol',tol,'AbsTol',tol);


%%% Initialize STM as column vector to be added with state
stm0_colVec = reshape(eye(6),36,1);
        

% % -------------------------------------------------
% %%% Southern Halo, ns1 .... (L2_L_1T  ns1)
% % -------------------------------------------------
% [T_SHalo__ns1, X_SHalo__ns1, t_ev, x_ev, in_ev] = ode113(@Int_CR3BnSTM, [0, PO_SHalo__ns1(c_Tp)], [PO_SHalo__ns1(c_x0:c_zd0)'; stm0_colVec], options_apsis, prms);
% figure; hold all
% plot3(X_SHalo__ns1(:,1),X_SHalo__ns1(:,2),X_SHalo__ns1(:,3), 'linewidth', 2, 'color', colors.black)
% plotSecondary(secondary)
% PlotBoi3_CR3Bn(26)
% axis equal
% title('SHalo ... ns1','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)
% 
% apses_alts = rowNorm(x_ev(:,1:3) - [1-prms.u,0,0]) - prms.R2;
% plot3(x_ev((apses_alts < 5/rNorm),1), x_ev((apses_alts < 5/rNorm),2), x_ev((apses_alts < 5/rNorm),3), '.','markersize',20, 'color', colors.red2)
% 
% latLons_SHalo__ns1 = zeros(length(T_SHalo__ns1),2);
% for kk = 1:length(T_SHalo__ns1)
%     [lat_deg, lon_deg] = BCR2latlon(X_SHalo__ns1(kk,1:3)', 'secondary', prms.u);
%     latLons_SHalo__ns1(kk,:) = [lat_deg, lon_deg];
% end
% latLons_apses_SHalo__ns1 = zeros(length(t_ev),2);
% for kk = 1:length(t_ev)
%     [lat_deg, lon_deg] = BCR2latlon(x_ev(kk,1:3)', 'secondary', prms.u);
%     latLons_apses_SHalo__ns1(kk,:) = [lat_deg, lon_deg];
% end
% 
% figure; hold all
% xlim([0 360])
% ylim([-90 90])
% PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')
% h = image(xlim, -ylim, bodies.europa.img);
% [lons_new] = convert_lon180_to_lon360(latLons_SHalo__ns1(:,2));
% plot(lons_new,latLons_SHalo__ns1(:,1),'k.')
% [lons_new_impact] = convert_lon180_to_lon360(latLons_apses_SHalo__ns1((apses_alts < 5/rNorm),2));
% p_landing = plot(lons_new_impact, latLons_apses_SHalo__ns1((apses_alts < 5/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.drkred, 'markerfacecolor', colors.red2);
% legend([p_landing], 'Tangent Impact', 'fontsize', 20, 'fontname', 'Times New Roman', 'location', 'best')
% title('SHalo ... ns1','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)
% 
% 
% % -------------------------------------------------
% %%% L2_L_1T_4P2_1P3__ns1
% % -------------------------------------------------
% [T_L2_L_1T_4P2_1P3__ns1, X_L2_L_1T_4P2_1P3__ns1, t_ev, x_ev, in_ev] = ode113(@Int_CR3BnSTM, [0, PO_L2_L_1T_4P2_1P3__ns1(c_Tp)], [PO_L2_L_1T_4P2_1P3__ns1(c_x0:c_zd0)'; stm0_colVec], options_apsis, prms);
% figure; hold all
% plot3(X_L2_L_1T_4P2_1P3__ns1(:,1),X_L2_L_1T_4P2_1P3__ns1(:,2),X_L2_L_1T_4P2_1P3__ns1(:,3), 'linewidth', 2, 'color', colors.black)
% plotSecondary(secondary)
% PlotBoi3_CR3Bn(26)
% axis equal
% title('L2\_L\_1T\_4P2\_1P3 ... ns1','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)
% 
% apses_alts = rowNorm(x_ev(:,1:3) - [1-prms.u,0,0]) - prms.R2;
% plot3(x_ev((apses_alts < 5/rNorm),1), x_ev((apses_alts < 5/rNorm),2), x_ev((apses_alts < 5/rNorm),3), '.','markersize',20, 'color', colors.red2)
% 
% latLons_L2_L_1T_4P2_1P3__ns1 = zeros(length(T_L2_L_1T_4P2_1P3__ns1),2);
% for kk = 1:length(T_L2_L_1T_4P2_1P3__ns1)
%     [lat_deg, lon_deg] = BCR2latlon(X_L2_L_1T_4P2_1P3__ns1(kk,1:3)', 'secondary', prms.u);
%     latLons_L2_L_1T_4P2_1P3__ns1(kk,:) = [lat_deg, lon_deg];
% end
% latLons_apses_L2_L_1T_4P2_1P3__ns1 = zeros(length(t_ev),2);
% for kk = 1:length(t_ev)
%     [lat_deg, lon_deg] = BCR2latlon(x_ev(kk,1:3)', 'secondary', prms.u);
%     latLons_apses_L2_L_1T_4P2_1P3__ns1(kk,:) = [lat_deg, lon_deg];
% end
% 
% figure; hold all
% xlim([0 360])
% ylim([-90 90])
% PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')
% h = image(xlim, -ylim, bodies.europa.img);
% [lons_new] = convert_lon180_to_lon360(latLons_L2_L_1T_4P2_1P3__ns1(:,2));
% plot(lons_new,latLons_L2_L_1T_4P2_1P3__ns1(:,1),'k.')
% [lons_new_impact] = convert_lon180_to_lon360(latLons_apses_L2_L_1T_4P2_1P3__ns1((apses_alts < 5/rNorm),2));
% p_landing = plot(lons_new_impact, latLons_apses_L2_L_1T_4P2_1P3__ns1((apses_alts < 5/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.drkred, 'markerfacecolor', colors.red2);
% legend([p_landing], 'Tangent Impact', 'fontsize', 20, 'fontname', 'Times New Roman', 'location', 'best')
% title('L2\_L\_1T\_4P2\_1P3 ... ns1','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)
% 
% % -------------------------------------------------
% %%% L2_L_1T_4P2_1P3__ns2
% % -------------------------------------------------
% [T_L2_L_1T_4P2_1P3__ns2, X_L2_L_1T_4P2_1P3__ns2, t_ev, x_ev, in_ev] = ode113(@Int_CR3BnSTM, [0, PO_L2_L_1T_4P2_1P3__ns2(c_Tp)], [PO_L2_L_1T_4P2_1P3__ns2(c_x0:c_zd0)'; stm0_colVec], options_apsis, prms);
% figure; hold all
% plot3(X_L2_L_1T_4P2_1P3__ns2(:,1),X_L2_L_1T_4P2_1P3__ns2(:,2),X_L2_L_1T_4P2_1P3__ns2(:,3), 'linewidth', 2, 'color', colors.black)
% plotSecondary(secondary)
% PlotBoi3_CR3Bn(26)
% axis equal
% title('L2\_L\_1T\_4P2\_1P3 ... ns2','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)
% 
% apses_alts = rowNorm(x_ev(:,1:3) - [1-prms.u,0,0]) - prms.R2;
% plot3(x_ev((apses_alts < 5/rNorm),1), x_ev((apses_alts < 5/rNorm),2), x_ev((apses_alts < 5/rNorm),3), '.','markersize',20, 'color', colors.red2)
% 
% latLons_L2_L_1T_4P2_1P3__ns2 = zeros(length(T_L2_L_1T_4P2_1P3__ns2),2);
% for kk = 1:length(T_L2_L_1T_4P2_1P3__ns2)
%     [lat_deg, lon_deg] = BCR2latlon(X_L2_L_1T_4P2_1P3__ns2(kk,1:3)', 'secondary', prms.u);
%     latLons_L2_L_1T_4P2_1P3__ns2(kk,:) = [lat_deg, lon_deg];
% end
% latLons_apses_L2_L_1T_4P2_1P3__ns2 = zeros(length(t_ev),2);
% for kk = 1:length(t_ev)
%     [lat_deg, lon_deg] = BCR2latlon(x_ev(kk,1:3)', 'secondary', prms.u);
%     latLons_apses_L2_L_1T_4P2_1P3__ns2(kk,:) = [lat_deg, lon_deg];
% end
% 
% 
% figure; hold all
% xlim([0 360])
% ylim([-90 90])
% PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')
% h = image(xlim, -ylim, bodies.europa.img);
% [lons_new] = convert_lon180_to_lon360(latLons_L2_L_1T_4P2_1P3__ns2(:,2));
% plot(lons_new,latLons_L2_L_1T_4P2_1P3__ns2(:,1),'k.')
% [lons_new_impact] = convert_lon180_to_lon360(latLons_apses_L2_L_1T_4P2_1P3__ns2((apses_alts < 5/rNorm),2));
% p_landing = plot(lons_new_impact, latLons_apses_L2_L_1T_4P2_1P3__ns2((apses_alts < 5/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.drkred, 'markerfacecolor', colors.red2);
% legend([p_landing], 'Tangent Impact', 'fontsize', 20, 'fontname', 'Times New Roman', 'location', 'best')
% title('L2\_L\_1T\_4P2\_1P3 ... ns2','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)
% 
% 
% % -------------------------------------------------
% %%% DRO__s1
% % -------------------------------------------------
% [T_DRO__s1, X_DRO__s1, t_ev, x_ev, in_ev] = ode113(@Int_CR3BnSTM, [0, PO_DRO__s1(c_Tp)], [PO_DRO__s1(c_x0:c_zd0)'; stm0_colVec], options_apsis, prms);
% figure; hold all
% plot3(X_DRO__s1(:,1),X_DRO__s1(:,2),X_DRO__s1(:,3), 'linewidth', 2, 'color', colors.black)
% plotSecondary(secondary)
% PlotBoi3_CR3Bn(26)
% axis equal
% title('DRO ... s1','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)
% 
% apses_alts = rowNorm(x_ev(:,1:3) - [1-prms.u,0,0]) - prms.R2;
% plot3(x_ev((apses_alts < 1/rNorm),1), x_ev((apses_alts < 1/rNorm),2), x_ev((apses_alts < 1/rNorm),3), '.','markersize',20, 'color', colors.red2)
% 
% latLons_DRO__s1 = zeros(length(T_DRO__s1),2);
% for kk = 1:length(T_DRO__s1)
%     [lat_deg, lon_deg] = BCR2latlon(X_DRO__s1(kk,1:3)', 'secondary', prms.u);
%     latLons_DRO__s1(kk,:) = [lat_deg, lon_deg];
% end
% latLons_apses_DRO__s1 = zeros(length(t_ev),2);
% for kk = 1:length(t_ev)
%     [lat_deg, lon_deg] = BCR2latlon(x_ev(kk,1:3)', 'secondary', prms.u);
%     latLons_apses_DRO__s1(kk,:) = [lat_deg, lon_deg];
% end
% 
% figure; hold all
% xlim([0 360])
% ylim([-90 90])
% PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')
% h = image(xlim, -ylim, bodies.europa.img);
% [lons_new] = convert_lon180_to_lon360(latLons_DRO__s1(:,2));
% plot(lons_new,latLons_DRO__s1(:,1),'k.')
% [lons_new_impact] = convert_lon180_to_lon360(latLons_apses_DRO__s1((apses_alts < 1/rNorm),2));
% p_landing = plot(lons_new_impact, latLons_apses_DRO__s1((apses_alts < 1/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.drkred, 'markerfacecolor', colors.red2);
% legend([p_landing], 'Tangent Impact', 'fontsize', 20, 'fontname', 'Times New Roman', 'location', 'best')
% title('DRO ... s1','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)
% 
% % -------------------------------------------------
% %%% DPO__s1
% % -------------------------------------------------
% [T_DPO__s1, X_DPO__s1, t_ev, x_ev, in_ev] = ode113(@Int_CR3BnSTM, [0, PO_DPO__s1(c_Tp)], [PO_DPO__s1(c_x0:c_zd0)'; stm0_colVec], options_apsis, prms);
% figure; hold all
% plot3(X_DPO__s1(:,1),X_DPO__s1(:,2),X_DPO__s1(:,3), 'linewidth', 2, 'color', colors.black)
% plotSecondary(secondary)
% PlotBoi3_CR3Bn(26)
% axis equal
% title('DPO ... s1','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)
% 
% apses_alts = rowNorm(x_ev(:,1:3) - [1-prms.u,0,0]) - prms.R2;
% plot3(x_ev((apses_alts < 5/rNorm),1), x_ev((apses_alts < 5/rNorm),2), x_ev((apses_alts < 5/rNorm),3), '.','markersize',20, 'color', colors.red2)
% 
% latLons_DPO__s1 = zeros(length(T_DPO__s1),2);
% for kk = 1:length(T_DPO__s1)
%     [lat_deg, lon_deg] = BCR2latlon(X_DPO__s1(kk,1:3)', 'secondary', prms.u);
%     latLons_DPO__s1(kk,:) = [lat_deg, lon_deg];
% end
% latLons_apses_DPO__s1 = zeros(length(t_ev),2);
% for kk = 1:length(t_ev)
%     [lat_deg, lon_deg] = BCR2latlon(x_ev(kk,1:3)', 'secondary', prms.u);
%     latLons_apses_DPO__s1(kk,:) = [lat_deg, lon_deg];
% end
% 
% 
% figure; hold all
% xlim([0 360])
% ylim([-90 90])
% PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')
% h = image(xlim, -ylim, bodies.europa.img);
% [lons_new] = convert_lon180_to_lon360(latLons_DPO__s1(:,2));
% plot(lons_new,latLons_DPO__s1(:,1),'k.')
% [lons_new_impact] = convert_lon180_to_lon360(latLons_apses_DPO__s1((apses_alts < 5/rNorm),2));
% p_landing = plot(lons_new_impact, latLons_apses_DPO__s1((apses_alts < 5/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.drkred, 'markerfacecolor', colors.red2);
% legend([p_landing], 'Tangent Impact', 'fontsize', 20, 'fontname', 'Times New Roman', 'location', 'best')
% title('DPO ... s1','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)
% 
% % -------------------------------------------------
% %%% DPO_2P3__ns1
% % -------------------------------------------------
% [T_DPO_2P3__ns1, X_DPO_2P3__ns1, t_ev, x_ev, in_ev] = ode113(@Int_CR3BnSTM, [0, PO_DPO_2P3__ns1(c_Tp)], [PO_DPO_2P3__ns1(c_x0:c_zd0)'; stm0_colVec], options_apsis, prms);
% figure; hold all
% plot3(X_DPO_2P3__ns1(:,1),X_DPO_2P3__ns1(:,2),X_DPO_2P3__ns1(:,3), 'linewidth', 2, 'color', colors.black)
% plotSecondary(secondary)
% PlotBoi3_CR3Bn(26)
% axis equal
% title('DPO\_2P3 ... ns1','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)
% 
% apses_alts = rowNorm(x_ev(:,1:3) - [1-prms.u,0,0]) - prms.R2;
% plot3(x_ev((apses_alts < 5/rNorm),1), x_ev((apses_alts < 5/rNorm),2), x_ev((apses_alts < 5/rNorm),3), '.','markersize',20, 'color', colors.red2)
% 
% latLons_DPO_2P3__ns1 = zeros(length(T_DPO_2P3__ns1),2);
% for kk = 1:length(T_DPO_2P3__ns1)
%     [lat_deg, lon_deg] = BCR2latlon(X_DPO_2P3__ns1(kk,1:3)', 'secondary', prms.u);
%     latLons_DPO_2P3__ns1(kk,:) = [lat_deg, lon_deg];
% end
% latLons_apses_DPO_2P3__ns1 = zeros(length(t_ev),2);
% for kk = 1:length(t_ev)
%     [lat_deg, lon_deg] = BCR2latlon(x_ev(kk,1:3)', 'secondary', prms.u);
%     latLons_apses_DPO_2P3__ns1(kk,:) = [lat_deg, lon_deg];
% end
% 
% figure; hold all
% xlim([0 360])
% ylim([-90 90])
% PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')
% h = image(xlim, -ylim, bodies.europa.img);
% [lons_new] = convert_lon180_to_lon360(latLons_DPO_2P3__ns1(:,2));
% plot(lons_new,latLons_DPO_2P3__ns1(:,1),'k.')
% [lons_new_impact] = convert_lon180_to_lon360(latLons_apses_DPO_2P3__ns1((apses_alts < 5/rNorm),2));
% p_landing = plot(lons_new_impact, latLons_apses_DPO_2P3__ns1((apses_alts < 5/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.drkred, 'markerfacecolor', colors.red2);
% legend([p_landing], 'Tangent Impact', 'fontsize', 20, 'fontname', 'Times New Roman', 'location', 'best')
% title('DPO\_2P3 ... ns1','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)
% 
% % -------------------------------------------------
% %%% DPO_2P4__ns1
% % -------------------------------------------------
% [T_DPO_2P4__ns1, X_DPO_2P4__ns1, t_ev, x_ev, in_ev] = ode113(@Int_CR3BnSTM, [0, PO_DPO_2P4__ns1(c_Tp)], [PO_DPO_2P4__ns1(c_x0:c_zd0)'; stm0_colVec], options_apsis, prms);
% figure; hold all
% plot3(X_DPO_2P4__ns1(:,1),X_DPO_2P4__ns1(:,2),X_DPO_2P4__ns1(:,3), 'linewidth', 2, 'color', colors.black)
% plotSecondary(secondary)
% PlotBoi3_CR3Bn(26)
% axis equal
% title('DPO\_2P4 ... ns1','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)
% 
% apses_alts = rowNorm(x_ev(:,1:3) - [1-prms.u,0,0]) - prms.R2;
% plot3(x_ev((apses_alts < 5/rNorm),1), x_ev((apses_alts < 5/rNorm),2), x_ev((apses_alts < 5/rNorm),3), '.','markersize',20, 'color', colors.red2)
% 
% latLons_DPO_2P4__ns1 = zeros(length(T_DPO_2P4__ns1),2);
% for kk = 1:length(T_DPO_2P4__ns1)
%     [lat_deg, lon_deg] = BCR2latlon(X_DPO_2P4__ns1(kk,1:3)', 'secondary', prms.u);
%     latLons_DPO_2P4__ns1(kk,:) = [lat_deg, lon_deg];
% end
% latLons_apses_DPO_2P4__ns1 = zeros(length(t_ev),2);
% for kk = 1:length(t_ev)
%     [lat_deg, lon_deg] = BCR2latlon(x_ev(kk,1:3)', 'secondary', prms.u);
%     latLons_apses_DPO_2P4__ns1(kk,:) = [lat_deg, lon_deg];
% end
% 
% figure; hold all
% xlim([0 360])
% ylim([-90 90])
% PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')
% h = image(xlim, -ylim, bodies.europa.img);
% [lons_new] = convert_lon180_to_lon360(latLons_DPO_2P4__ns1(:,2));
% plot(lons_new,latLons_DPO_2P4__ns1(:,1),'k.')
% [lons_new_impact] = convert_lon180_to_lon360(latLons_apses_DPO_2P4__ns1((apses_alts < 5/rNorm),2));
% p_landing = plot(lons_new_impact, latLons_apses_DPO_2P4__ns1((apses_alts < 5/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.drkred, 'markerfacecolor', colors.red2);
% legend([p_landing], 'Tangent Impact', 'fontsize', 20, 'fontname', 'Times New Roman', 'location', 'best')
% title('DPO\_2P4 ... ns1','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)
% 
% % -------------------------------------------------
% %%% Hg2__ns1
% % -------------------------------------------------
% [T_Hg2__ns1, X_Hg2__ns1, t_ev, x_ev, in_ev] = ode113(@Int_CR3BnSTM, [0, PO_Hg2__ns1(c_Tp)], [PO_Hg2__ns1(c_x0:c_zd0)'; stm0_colVec], options_apsis, prms);
% figure; hold all
% plot3(X_Hg2__ns1(:,1),X_Hg2__ns1(:,2),X_Hg2__ns1(:,3), 'linewidth', 2, 'color', colors.black)
% plotSecondary(secondary)
% PlotBoi3_CR3Bn(26)
% axis equal
% title('Hg2 ... ns1','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)
% 
% apses_alts = rowNorm(x_ev(:,1:3) - [1-prms.u,0,0]) - prms.R2;
% plot3(x_ev((apses_alts < 5/rNorm),1), x_ev((apses_alts < 5/rNorm),2), x_ev((apses_alts < 5/rNorm),3), '.','markersize',20, 'color', colors.red2)
% 
% latLons_Hg2__ns1 = zeros(length(T_Hg2__ns1),2);
% for kk = 1:length(T_Hg2__ns1)
%     [lat_deg, lon_deg] = BCR2latlon(X_Hg2__ns1(kk,1:3)', 'secondary', prms.u);
%     latLons_Hg2__ns1(kk,:) = [lat_deg, lon_deg];
% end
% latLons_apses_Hg2__ns1 = zeros(length(t_ev),2);
% for kk = 1:length(t_ev)
%     [lat_deg, lon_deg] = BCR2latlon(x_ev(kk,1:3)', 'secondary', prms.u);
%     latLons_apses_Hg2__ns1(kk,:) = [lat_deg, lon_deg];
% end
% 
% figure; hold all
% xlim([0 360])
% ylim([-90 90])
% PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')
% h = image(xlim, -ylim, bodies.europa.img);
% [lons_new] = convert_lon180_to_lon360(latLons_Hg2__ns1(:,2));
% plot(lons_new,latLons_Hg2__ns1(:,1),'k.')
% [lons_new_impact] = convert_lon180_to_lon360(latLons_apses_Hg2__ns1((apses_alts < 5/rNorm),2));
% p_landing = plot(lons_new_impact, latLons_apses_Hg2__ns1((apses_alts < 5/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.drkred, 'markerfacecolor', colors.red2);
% legend([p_landing], 'Tangent Impact', 'fontsize', 20, 'fontname', 'Times New Roman', 'location', 'best')
% title('Hg2 ... ns1','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)
% 
% % -------------------------------------------------
% %%% Hg2_2P2__ns1
% % -------------------------------------------------
% [T_Hg2_2P2__ns1, X_Hg2_2P2__ns1, t_ev, x_ev, in_ev] = ode113(@Int_CR3BnSTM, [0, PO_Hg2_2P2__ns1(c_Tp)], [PO_Hg2_2P2__ns1(c_x0:c_zd0)'; stm0_colVec], options_apsis, prms);
% figure; hold all
% plot3(X_Hg2_2P2__ns1(:,1),X_Hg2_2P2__ns1(:,2),X_Hg2_2P2__ns1(:,3), 'linewidth', 2, 'color', colors.black)
% plotSecondary(secondary)
% PlotBoi3_CR3Bn(26)
% axis equal
% title('Hg2\_2P2 ... ns1','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)
% 
% apses_alts = rowNorm(x_ev(:,1:3) - [1-prms.u,0,0]) - prms.R2;
% 
% plot3(x_ev((apses_alts < 5/rNorm),1), x_ev((apses_alts < 5/rNorm),2), x_ev((apses_alts < 5/rNorm),3), '.','markersize',20, 'color', colors.red2)
% plot3(x_ev(((136/rNorm < apses_alts) & (apses_alts < 138/rNorm)),1), x_ev(((136/rNorm < apses_alts) & (apses_alts < 138/rNorm)),2), x_ev(((136/rNorm < apses_alts) & (apses_alts < 138/rNorm)),3), '.','markersize',20, 'color', colors.purp)
% 
% latLons_Hg2_2P2__ns1 = zeros(length(T_Hg2_2P2__ns1),2);
% for kk = 1:length(T_Hg2_2P2__ns1)
%     [lat_deg, lon_deg] = BCR2latlon(X_Hg2_2P2__ns1(kk,1:3)', 'secondary', prms.u);
%     latLons_Hg2_2P2__ns1(kk,:) = [lat_deg, lon_deg];
% end
% latLons_apses_Hg2_2P2__ns1 = zeros(length(t_ev),2);
% for kk = 1:length(t_ev)
%     [lat_deg, lon_deg] = BCR2latlon(x_ev(kk,1:3)', 'secondary', prms.u);
%     latLons_apses_Hg2_2P2__ns1(kk,:) = [lat_deg, lon_deg];
% end
% 
% 
% figure; hold all
% xlim([0 360])
% ylim([-90 90])
% PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')
% h = image(xlim, -ylim, bodies.europa.img);
% [lons_new] = convert_lon180_to_lon360(latLons_Hg2_2P2__ns1(:,2));
% plot(lons_new,latLons_Hg2_2P2__ns1(:,1),'k.')
% [lons_new_impact] = convert_lon180_to_lon360(latLons_apses_Hg2_2P2__ns1((apses_alts < 5/rNorm),2));
% p_landing = plot(lons_new_impact, latLons_apses_Hg2_2P2__ns1((apses_alts < 5/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.drkred, 'markerfacecolor', colors.red2);
% [lons_new_impact_close] = convert_lon180_to_lon360(latLons_apses_Hg2_2P2__ns1(((136/rNorm < apses_alts) & (apses_alts < 138/rNorm)),2));
% p_close = plot(lons_new_impact_close, latLons_apses_Hg2_2P2__ns1(((136/rNorm < apses_alts) & (apses_alts < 138/rNorm)),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.purp, 'markerfacecolor', colors.ltpurp);
% legend([p_landing, p_close], 'Tangent Impact', '135 km Periapsis', 'fontsize', 20, 'fontname', 'Times New Roman')
% title('Hg2\_2P2 ... ns1','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)
% 
% % -------------------------------------------------
% %%% Hg2_2T__ns1
% % -------------------------------------------------
[T_Hg2_2T__ns1, X_Hg2_2T__ns1, t_ev, x_ev, in_ev] = ode113(@Int_CR3BnSTM, [0, PO_Hg2_2T__ns1(c_Tp)], [PO_Hg2_2T__ns1(c_x0:c_zd0)'; stm0_colVec], options_apsis, prms);
figure; hold all
plot3(X_Hg2_2T__ns1(:,1),X_Hg2_2T__ns1(:,2),X_Hg2_2T__ns1(:,3), 'linewidth', 2, 'color', colors.red2)
plotTrajShadows(X_Hg2_2T__ns1, 2, colors.grey, 'x', 1.025, 'y', 0.025, 'z', -1.7e-2, 'bodyshadow', [1-prms.u, prms.R2])
plotSecondary(secondary)
PlotBoi3_CR3Bn(26)
axis equal
title('Hg2\_2T ... ns1','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)


apses_alts = rowNorm(x_ev(:,1:3) - [1-prms.u,0,0]) - prms.R2;
plot3(x_ev((apses_alts < 5/rNorm),1), x_ev((apses_alts < 5/rNorm),2), x_ev((apses_alts < 5/rNorm),3), '.','markersize',20, 'color', colors.black)
plot3(x_ev(((212/rNorm < apses_alts) & (apses_alts < 213/rNorm)),1), x_ev(((212/rNorm < apses_alts) & (apses_alts < 213/rNorm)),2), x_ev(((212/rNorm < apses_alts) & (apses_alts < 213/rNorm)),3), '.','markersize',20, 'color', colors.blue2)

latLons_Hg2_2T__ns1 = zeros(length(T_Hg2_2T__ns1),2);
for kk = 1:length(T_Hg2_2T__ns1)
    [lat_deg, lon_deg] = BCR2latlon(X_Hg2_2T__ns1(kk,1:3)', 'secondary', prms.u);
    latLons_Hg2_2T__ns1(kk,:) = [lat_deg, lon_deg];
end
latLons_apses_Hg2_2T__ns1 = zeros(length(t_ev),2);
for kk = 1:length(t_ev)
    [lat_deg, lon_deg] = BCR2latlon(x_ev(kk,1:3)', 'secondary', prms.u);
    latLons_apses_Hg2_2T__ns1(kk,:) = [lat_deg, lon_deg];
end

figure; hold all
xlim([0 360])
ylim([-90 90])
PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')
h = image(xlim, -ylim, bodies.europa.img);
[lons_new] = convert_lon180_to_lon360(latLons_Hg2_2T__ns1(:,2));
plot(lons_new,latLons_Hg2_2T__ns1(:,1),'.', 'color', colors.red2)
[lons_new_impact] = convert_lon180_to_lon360(latLons_apses_Hg2_2T__ns1((apses_alts < 5/rNorm),2));
p_landing = plot(lons_new_impact, latLons_apses_Hg2_2T__ns1((apses_alts < 5/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.black, 'markerfacecolor', colors.drkgrey);
[lons_new_impact_close] = convert_lon180_to_lon360(latLons_apses_Hg2_2T__ns1(((212/rNorm < apses_alts) & (apses_alts < 213/rNorm)),2));
p_close = plot(lons_new_impact_close, latLons_apses_Hg2_2T__ns1(((212/rNorm < apses_alts) & (apses_alts < 213/rNorm)),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.drkblue, 'markerfacecolor', colors.blue2);
legend([p_landing, p_close], 'Tangent Impact', '212 km Periapsis', 'fontsize', 20, 'fontname', 'Times New Roman')
title('Hg2\_2T ... ns1','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)

% 
% % -------------------------------------------------
% %%% Hg2_3T__ns1
% % -------------------------------------------------
[T_Hg2_3T__ns1, X_Hg2_3T__ns1, t_ev, x_ev, in_ev] = ode113(@Int_CR3BnSTM, [0, PO_Hg2_3T__ns1(c_Tp)], [PO_Hg2_3T__ns1(c_x0:c_zd0)'; stm0_colVec], options_apsis, prms);
figure; hold all
plot3(X_Hg2_3T__ns1(:,1),X_Hg2_3T__ns1(:,2),X_Hg2_3T__ns1(:,3), 'linewidth', 2, 'color', colors.red2)
plotTrajShadows(X_Hg2_3T__ns1, 2, colors.grey, 'x', 1.025, 'y', 0.018, 'z', -1.4e-2, 'bodyshadow', [1-prms.u, prms.R2])
plotSecondary(secondary)
PlotBoi3_CR3Bn(26)
axis equal
title('Hg2\_3T ... ns1','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)

apses_alts = rowNorm(x_ev(:,1:3) - [1-prms.u,0,0]) - prms.R2;
plot3(x_ev((apses_alts < 5/rNorm),1), x_ev((apses_alts < 5/rNorm),2), x_ev((apses_alts < 5/rNorm),3), '.','markersize',20, 'color', colors.black)

latLons_Hg2_3T__ns1 = zeros(length(T_Hg2_3T__ns1),2);
for kk = 1:length(T_Hg2_3T__ns1)
    [lat_deg, lon_deg] = BCR2latlon(X_Hg2_3T__ns1(kk,1:3)', 'secondary', prms.u);
    latLons_Hg2_3T__ns1(kk,:) = [lat_deg, lon_deg];
end
latLons_apses_Hg2_3T__ns1 = zeros(length(t_ev),2);
for kk = 1:length(t_ev)
    [lat_deg, lon_deg] = BCR2latlon(x_ev(kk,1:3)', 'secondary', prms.u);
    latLons_apses_Hg2_3T__ns1(kk,:) = [lat_deg, lon_deg];
end

figure; hold all
xlim([0 360])
ylim([-90 90])
PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')
h = image(xlim, -ylim, bodies.europa.img);
[lons_new] = convert_lon180_to_lon360(latLons_Hg2_3T__ns1(:,2));
plot(lons_new,latLons_Hg2_3T__ns1(:,1),'.', 'color', colors.red2)
[lons_new_impact] = convert_lon180_to_lon360(latLons_apses_Hg2_3T__ns1((apses_alts < 5/rNorm),2));
p_landing = plot(lons_new_impact, latLons_apses_Hg2_3T__ns1((apses_alts < 5/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.black, 'markerfacecolor', colors.drkgrey);
legend([p_landing], 'Tangent Impact', 'fontsize', 20, 'fontname', 'Times New Roman', 'location', 'best')
title('Hg2\_3T ... ns1','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)
% 
% % -------------------------------------------------
% %%% Hg2_3T__ns2
% % -------------------------------------------------
% [T_Hg2_3T__ns2, X_Hg2_3T__ns2, t_ev, x_ev, in_ev] = ode113(@Int_CR3BnSTM, [0, PO_Hg2_3T__ns2(c_Tp)], [PO_Hg2_3T__ns2(c_x0:c_zd0)'; stm0_colVec], options_apsis, prms);
% figure; hold all
% plot3(X_Hg2_3T__ns2(:,1),X_Hg2_3T__ns2(:,2),X_Hg2_3T__ns2(:,3), 'linewidth', 2, 'color', colors.black)
% plotSecondary(secondary)
% PlotBoi3_CR3Bn(26)
% axis equal
% title('Hg2\_3T ... ns2','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)
% 
% apses_alts = rowNorm(x_ev(:,1:3) - [1-prms.u,0,0]) - prms.R2;
% plot3(x_ev((apses_alts < 8/rNorm),1), x_ev((apses_alts < 8/rNorm),2), x_ev((apses_alts < 8/rNorm),3), '.','markersize',20, 'color', colors.red2)
% 
% latLons_Hg2_3T__ns2 = zeros(length(T_Hg2_3T__ns2),2);
% for kk = 1:length(T_Hg2_3T__ns2)
%     [lat_deg, lon_deg] = BCR2latlon(X_Hg2_3T__ns2(kk,1:3)', 'secondary', prms.u);
%     latLons_Hg2_3T__ns2(kk,:) = [lat_deg, lon_deg];
% end
% latLons_apses_Hg2_3T__ns2 = zeros(length(t_ev),2);
% for kk = 1:length(t_ev)
%     [lat_deg, lon_deg] = BCR2latlon(x_ev(kk,1:3)', 'secondary', prms.u);
%     latLons_apses_Hg2_3T__ns2(kk,:) = [lat_deg, lon_deg];
% end
% 
% figure; hold all
% xlim([0 360])
% ylim([-90 90])
% PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')
% h = image(xlim, -ylim, bodies.europa.img);
% [lons_new] = convert_lon180_to_lon360(latLons_Hg2_3T__ns2(:,2));
% plot(lons_new,latLons_Hg2_3T__ns2(:,1),'k.')
% [lons_new_impact] = convert_lon180_to_lon360(latLons_apses_Hg2_3T__ns2((apses_alts < 8/rNorm),2));
% p_landing = plot(lons_new_impact, latLons_apses_Hg2_3T__ns2((apses_alts < 8/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.drkred, 'markerfacecolor', colors.red2);
% legend([p_landing], 'Tangent Impact', 'fontsize', 20, 'fontname', 'Times New Roman', 'location', 'best')
% title('Hg2\_3T ... ns2','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)
% 
% 
% % -------------------------------------------------
% %%% Se7__s1
% % -------------------------------------------------
% [T_Se7__s1, X_Se7__s1, t_ev, x_ev, in_ev] = ode113(@Int_CR3BnSTM, [0, PO_Se7__s1(c_Tp)], [PO_Se7__s1(c_x0:c_zd0)'; stm0_colVec], options_apsis, prms);
% figure; hold all
% plot3(X_Se7__s1(:,1),X_Se7__s1(:,2),X_Se7__s1(:,3), 'linewidth', 2, 'color', colors.black)
% plotSecondary(secondary)
% PlotBoi3_CR3Bn(26)
% axis equal
% title('Se7 ... s1','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)
% 
% apses_alts = rowNorm(x_ev(:,1:3) - [1-prms.u,0,0]) - prms.R2;
% plot3(x_ev((apses_alts < 5/rNorm),1), x_ev((apses_alts < 5/rNorm),2), x_ev((apses_alts < 5/rNorm),3), '.','markersize',20, 'color', colors.red2)
% plot3(x_ev(((79/rNorm < apses_alts) & (apses_alts < 80/rNorm)),1), x_ev(((79/rNorm < apses_alts) & (apses_alts < 80/rNorm)),2), x_ev(((79/rNorm < apses_alts) & (apses_alts < 80/rNorm)),3), '.','markersize',20, 'color', colors.purp)
% 
% latLons_Se7__s1 = zeros(length(T_Se7__s1),2);
% for kk = 1:length(T_Se7__s1)
%     [lat_deg, lon_deg] = BCR2latlon(X_Se7__s1(kk,1:3)', 'secondary', prms.u);
%     latLons_Se7__s1(kk,:) = [lat_deg, lon_deg];
% end
% latLons_apses_Se7__s1 = zeros(length(t_ev),2);
% for kk = 1:length(t_ev)
%     [lat_deg, lon_deg] = BCR2latlon(x_ev(kk,1:3)', 'secondary', prms.u);
%     latLons_apses_Se7__s1(kk,:) = [lat_deg, lon_deg];
% end
% 
% figure; hold all
% xlim([0 360])
% ylim([-90 90])
% PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')
% h = image(xlim, -ylim, bodies.europa.img);
% [lons_new] = convert_lon180_to_lon360(latLons_Se7__s1(:,2));
% plot(lons_new,latLons_Se7__s1(:,1),'k.')
% [lons_new_impact] = convert_lon180_to_lon360(latLons_apses_Se7__s1((apses_alts < 5/rNorm),2));
% p_landing = plot(lons_new_impact, latLons_apses_Se7__s1((apses_alts < 5/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.drkred, 'markerfacecolor', colors.red2);
% % [lons_new_impact_close] = convert_lon180_to_lon360(latLons_apses_Se7__s1(((116/rNorm < apses_alts) & (apses_alts < 117/rNorm)),2));
% [lons_new_impact_close] = convert_lon180_to_lon360(latLons_apses_Se7__s1(((79/rNorm < apses_alts) & (apses_alts < 80/rNorm)),2));
% p_close = plot(lons_new_impact_close, latLons_apses_Se7__s1(((79/rNorm < apses_alts) & (apses_alts < 80/rNorm)),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.purp, 'markerfacecolor', colors.ltpurp);
% legend([p_landing, p_close], 'Tangent Impact', '80 km Periapsis', 'fontsize', 20, 'fontname', 'Times New Roman')
% title('Se7 ... s1','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)

% % -------------------------------------------------
% %%% Se7__ns1
% % -------------------------------------------------
% [T_Se7__ns1, X_Se7__ns1, t_ev, x_ev, in_ev] = ode113(@Int_CR3BnSTM, [0, PO_Se7__ns1(c_Tp)], [PO_Se7__ns1(c_x0:c_zd0)'; stm0_colVec], options_apsis, prms);
% figure; hold all
% plot3(X_Se7__ns1(:,1),X_Se7__ns1(:,2),X_Se7__ns1(:,3), 'linewidth', 2, 'color', colors.black)
% plotSecondary(secondary)
% PlotBoi3_CR3Bn(26)
% axis equal
% title('Se7 ... ns1','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)
% 
% 
% apses_alts = rowNorm(x_ev(:,1:3) - [1-prms.u,0,0]) - prms.R2;
% plot3(x_ev((apses_alts < 5/rNorm),1), x_ev((apses_alts < 5/rNorm),2), x_ev((apses_alts < 5/rNorm),3), '.','markersize',20, 'color', colors.red2)
% plot3(x_ev(((116/rNorm < apses_alts) & (apses_alts < 117/rNorm)),1), x_ev(((116/rNorm < apses_alts) & (apses_alts < 117/rNorm)),2), x_ev(((116/rNorm < apses_alts) & (apses_alts < 117/rNorm)),3), '.','markersize',20, 'color', colors.purp)
% 
% latLons_Se7__ns1 = zeros(length(T_Se7__ns1),2);
% for kk = 1:length(T_Se7__ns1)
%     [lat_deg, lon_deg] = BCR2latlon(X_Se7__ns1(kk,1:3)', 'secondary', prms.u);
%     latLons_Se7__ns1(kk,:) = [lat_deg, lon_deg];
% end
% latLons_apses_Se7__ns1 = zeros(length(t_ev),2);
% for kk = 1:length(t_ev)
%     [lat_deg, lon_deg] = BCR2latlon(x_ev(kk,1:3)', 'secondary', prms.u);
%     latLons_apses_Se7__ns1(kk,:) = [lat_deg, lon_deg];
% end
% 
% figure; hold all
% xlim([0 360])
% ylim([-90 90])
% PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')
% h = image(xlim, -ylim, bodies.europa.img);
% [lons_new] = convert_lon180_to_lon360(latLons_Se7__ns1(:,2));
% plot(lons_new,latLons_Se7__ns1(:,1),'k.')
% [lons_new_impact] = convert_lon180_to_lon360(latLons_apses_Se7__ns1((apses_alts < 5/rNorm),2));
% p_landing = plot(lons_new_impact, latLons_apses_Se7__ns1((apses_alts < 5/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.drkred, 'markerfacecolor', colors.red2);
% [lons_new_impact_close] = convert_lon180_to_lon360(latLons_apses_Se7__ns1(((116/rNorm < apses_alts) & (apses_alts < 117/rNorm)),2));
% p_close = plot(lons_new_impact_close, latLons_apses_Se7__ns1(((116/rNorm < apses_alts) & (apses_alts < 117/rNorm)),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.purp, 'markerfacecolor', colors.ltpurp);
% legend([p_landing, p_close], 'Tangent Impact', '116 km Periapsis', 'fontsize', 20, 'fontname', 'Times New Roman')
% title('Se7 ... ns1','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)
% 
% -------------------------------------------------
%%% Se7_6P2__ns1
% -------------------------------------------------
[T_Se7_6P2__ns1, X_Se7_6P2__ns1, t_ev, x_ev, in_ev] = ode113(@Int_CR3BnSTM, [0, PO_Se7_6P2__ns1(c_Tp)], [PO_Se7_6P2__ns1(c_x0:c_zd0)'; stm0_colVec], options_apsis, prms);
figure; hold all
plot3(X_Se7_6P2__ns1(:,1),X_Se7_6P2__ns1(:,2),X_Se7_6P2__ns1(:,3), 'linewidth', 2, 'color', colors.red2)
plotTrajShadows(X_Se7_6P2__ns1, 2, colors.grey, 'x', 1.025, 'y', 0.018, 'z', -1.4e-2, 'bodyshadow', [1-prms.u, prms.R2])
plotSecondary(secondary)
PlotBoi3_CR3Bn(26)
axis equal
title('Se7\_6P2 ... ns1','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)

apses_alts = rowNorm(x_ev(:,1:3) - [1-prms.u,0,0]) - prms.R2;
plot3(x_ev((apses_alts < 5/rNorm),1), x_ev((apses_alts < 5/rNorm),2), x_ev((apses_alts < 5/rNorm),3), '.','markersize',20, 'color', colors.black)
plot3(x_ev(((124/rNorm < apses_alts) & (apses_alts < 125/rNorm)),1), x_ev(((124/rNorm < apses_alts) & (apses_alts < 125/rNorm)),2), x_ev(((124/rNorm < apses_alts) & (apses_alts < 125/rNorm)),3), '.','markersize',20, 'color', colors.blue2)

latLons_Se7_6P2__ns1 = zeros(length(T_Se7_6P2__ns1),2);
for kk = 1:length(T_Se7_6P2__ns1)
    [lat_deg, lon_deg] = BCR2latlon(X_Se7_6P2__ns1(kk,1:3)', 'secondary', prms.u);
    latLons_Se7_6P2__ns1(kk,:) = [lat_deg, lon_deg];
end
latLons_apses_Se7_6P2__ns1 = zeros(length(t_ev),2);
for kk = 1:length(t_ev)
    [lat_deg, lon_deg] = BCR2latlon(x_ev(kk,1:3)', 'secondary', prms.u);
    latLons_apses_Se7_6P2__ns1(kk,:) = [lat_deg, lon_deg];
end

figure; hold all
xlim([0 360])
ylim([-90 90])
PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')
h = image(xlim, -ylim, bodies.europa.img);
[lons_new] = convert_lon180_to_lon360(latLons_Se7_6P2__ns1(:,2));
plot(lons_new,latLons_Se7_6P2__ns1(:,1),'.', 'color', colors.red2)
[lons_new_impact] = convert_lon180_to_lon360(latLons_apses_Se7_6P2__ns1((apses_alts < 5/rNorm),2));
p_landing = plot(lons_new_impact, latLons_apses_Se7_6P2__ns1((apses_alts < 5/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.black, 'markerfacecolor', colors.drkgrey);
[lons_new_impact_close] = convert_lon180_to_lon360(latLons_apses_Se7_6P2__ns1(((124/rNorm < apses_alts) & (apses_alts < 125/rNorm)),2));
p_close = plot(lons_new_impact_close, latLons_apses_Se7_6P2__ns1(((124/rNorm < apses_alts) & (apses_alts < 125/rNorm)),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.drkblue, 'markerfacecolor', colors.blue2);
legend([p_landing, p_close], 'Tangent Impact', '124 km Periapsis', 'fontsize', 20, 'fontname', 'Times New Roman')
title('Se7\_6P2 ... ns1','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)
% 
% 
% % -------------------------------------------------
% %%% LoPO__s1
% % -------------------------------------------------
% [T_LoPO__s1, X_LoPO__s1, t_ev, x_ev, in_ev] = ode113(@Int_CR3BnSTM, [0, PO_LoPO__s1(c_Tp)], [PO_LoPO__s1(c_x0:c_zd0)'; stm0_colVec], options_apsis, prms);
% figure; hold all
% plot3(X_LoPO__s1(:,1),X_LoPO__s1(:,2),X_LoPO__s1(:,3), 'linewidth', 2, 'color', colors.black)
% plotSecondary(secondary)
% PlotBoi3_CR3Bn(26)
% axis equal
% title('LoPO ... s1','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)
% 
% 
% apses_alts = rowNorm(x_ev(:,1:3) - [1-prms.u,0,0]) - prms.R2;
% plot3(x_ev((apses_alts < 5/rNorm),1), x_ev((apses_alts < 5/rNorm),2), x_ev((apses_alts < 5/rNorm),3), '.','markersize',20, 'color', colors.red2)
% 
% latLons_LoPO__s1 = zeros(length(T_LoPO__s1),2);
% for kk = 1:length(T_LoPO__s1)
%     [lat_deg, lon_deg] = BCR2latlon(X_LoPO__s1(kk,1:3)', 'secondary', prms.u);
%     latLons_LoPO__s1(kk,:) = [lat_deg, lon_deg];
% end
% latLons_apses_LoPO__s1 = zeros(length(t_ev),2);
% for kk = 1:length(t_ev)
%     [lat_deg, lon_deg] = BCR2latlon(x_ev(kk,1:3)', 'secondary', prms.u);
%     latLons_apses_LoPO__s1(kk,:) = [lat_deg, lon_deg];
% end
% 
% figure; hold all
% xlim([0 360])
% ylim([-90 90])
% PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')
% h = image(xlim, -ylim, bodies.europa.img);
% [lons_new] = convert_lon180_to_lon360(latLons_LoPO__s1(:,2));
% plot(lons_new,latLons_LoPO__s1(:,1),'k.')
% [lons_new_impact] = convert_lon180_to_lon360(latLons_apses_LoPO__s1((apses_alts < 5/rNorm),2));
% p_landing = plot(lons_new_impact, latLons_apses_LoPO__s1((apses_alts < 5/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.drkred, 'markerfacecolor', colors.red2);
% legend([p_landing], 'Tangent Impact', 'fontsize', 20, 'fontname', 'Times New Roman', 'location', 'best')
% title('LoPO ... s1','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)
% 
% % -------------------------------------------------
% %%% LoPO_2P2_1T__ns1
% % -------------------------------------------------
% [T_LoPO_2P2_1T__ns1, X_LoPO_2P2_1T__ns1, t_ev, x_ev, in_ev] = ode113(@Int_CR3BnSTM, [0, PO_LoPO_2P2_1T__ns1(c_Tp)], [PO_LoPO_2P2_1T__ns1(c_x0:c_zd0)'; stm0_colVec], options_apsis, prms);
% figure; hold all
% plot3(X_LoPO_2P2_1T__ns1(:,1),X_LoPO_2P2_1T__ns1(:,2),X_LoPO_2P2_1T__ns1(:,3), 'linewidth', 2, 'color', colors.black)
% plotSecondary(secondary)
% PlotBoi3_CR3Bn(26)
% axis equal
% title('LoPO\_2P2\_1T ... ns1','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)
% 
% apses_alts = rowNorm(x_ev(:,1:3) - [1-prms.u,0,0]) - prms.R2;
% plot3(x_ev((apses_alts < 5/rNorm),1), x_ev((apses_alts < 5/rNorm),2), x_ev((apses_alts < 5/rNorm),3), '.','markersize',20, 'color', colors.red2)
% 
% latLons_LoPO_2P2_1T__ns1 = zeros(length(T_LoPO_2P2_1T__ns1),2);
% for kk = 1:length(T_LoPO_2P2_1T__ns1)
%     [lat_deg, lon_deg] = BCR2latlon(X_LoPO_2P2_1T__ns1(kk,1:3)', 'secondary', prms.u);
%     latLons_LoPO_2P2_1T__ns1(kk,:) = [lat_deg, lon_deg];
% end
% latLons_apses_LoPO_2P2_1T__ns1 = zeros(length(t_ev),2);
% for kk = 1:length(t_ev)
%     [lat_deg, lon_deg] = BCR2latlon(x_ev(kk,1:3)', 'secondary', prms.u);
%     latLons_apses_LoPO_2P2_1T__ns1(kk,:) = [lat_deg, lon_deg];
% end
% 
% figure; hold all
% xlim([0 360])
% ylim([-90 90])
% PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')
% h = image(xlim, -ylim, bodies.europa.img);
% [lons_new] = convert_lon180_to_lon360(latLons_LoPO_2P2_1T__ns1(:,2));
% plot(lons_new,latLons_LoPO_2P2_1T__ns1(:,1),'k.')
% [lons_new_impact] = convert_lon180_to_lon360(latLons_apses_LoPO_2P2_1T__ns1((apses_alts < 5/rNorm),2));
% p_landing = plot(lons_new_impact, latLons_apses_LoPO_2P2_1T__ns1((apses_alts < 5/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.drkred, 'markerfacecolor', colors.red2);
% legend([p_landing], 'Tangent Impact', 'fontsize', 20, 'fontname', 'Times New Roman', 'location', 'best')
% title('LoPO\_2P2\_1T ... ns1','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)
% 
% % -------------------------------------------------
% %%% LoPO_2P4__ns1
% % -------------------------------------------------
% [T_LoPO_2P4__ns1, X_LoPO_2P4__ns1, t_ev, x_ev, in_ev] = ode113(@Int_CR3BnSTM, [0, PO_LoPO_2P4__ns1(c_Tp)], [PO_LoPO_2P4__ns1(c_x0:c_zd0)'; stm0_colVec], options_apsis, prms);
% figure; hold all
% plot3(X_LoPO_2P4__ns1(:,1),X_LoPO_2P4__ns1(:,2),X_LoPO_2P4__ns1(:,3), 'linewidth', 2, 'color', colors.black)
% plotSecondary(secondary)
% PlotBoi3_CR3Bn(26)
% axis equal
% title('LoPO\_2P4 ... ns1','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)
% 
% apses_alts = rowNorm(x_ev(:,1:3) - [1-prms.u,0,0]) - prms.R2;
% plot3(x_ev((apses_alts < 5/rNorm),1), x_ev((apses_alts < 5/rNorm),2), x_ev((apses_alts < 5/rNorm),3), '.','markersize',20, 'color', colors.red2)
% 
% latLons_LoPO_2P4__ns1 = zeros(length(T_LoPO_2P4__ns1),2);
% for kk = 1:length(T_LoPO_2P4__ns1)
%     [lat_deg, lon_deg] = BCR2latlon(X_LoPO_2P4__ns1(kk,1:3)', 'secondary', prms.u);
%     latLons_LoPO_2P4__ns1(kk,:) = [lat_deg, lon_deg];
% end
% latLons_apses_LoPO_2P4__ns1 = zeros(length(t_ev),2);
% for kk = 1:length(t_ev)
%     [lat_deg, lon_deg] = BCR2latlon(x_ev(kk,1:3)', 'secondary', prms.u);
%     latLons_apses_LoPO_2P4__ns1(kk,:) = [lat_deg, lon_deg];
% end
% 
% figure; hold all
% xlim([0 360])
% ylim([-90 90])
% PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')
% h = image(xlim, -ylim, bodies.europa.img);
% [lons_new] = convert_lon180_to_lon360(latLons_LoPO_2P4__ns1(:,2));
% plot(lons_new,latLons_LoPO_2P4__ns1(:,1),'k.')
% [lons_new_impact] = convert_lon180_to_lon360(latLons_apses_LoPO_2P4__ns1((apses_alts < 5/rNorm),2));
% p_landing = plot(lons_new_impact, latLons_apses_LoPO_2P4__ns1((apses_alts < 5/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.drkred, 'markerfacecolor', colors.red2);
% legend([p_landing], 'Tangent Impact', 'fontsize', 20, 'fontname', 'Times New Roman', 'location', 'best')
% title('LoPO\_2P4 ... ns1','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)








