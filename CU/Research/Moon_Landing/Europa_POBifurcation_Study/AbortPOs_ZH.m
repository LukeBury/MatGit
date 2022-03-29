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
%%% Run Switches
% ========================================================================
run_L2_L_1T_ns1          = true;

run_L2_L_1T_4P2_1P3__ns1 = true;

run_L2_L_1T_4P2_1P3__ns2 = true;

run_DRO__s1              = true;

run_DPO__s1              = true;

run_DPO_2P3__ns1         = true;

run_DPO_2P4__ns1         = true;

run_Hg2__ns1             = true;

run_Hg2_2P2__ns1         = true;

run_Hg2_3T__ns2          = true;

run_Se7__s1              = true;

run_Se7__ns1             = true;

run_LoPO__s1             = true;

run_LoPO_2P2_1T__ns1     = true;

run_LoPO_2P4__ns1        = true;

% ========================================================================
%%% Importing Data
% ========================================================================
%%% General data on solar system bodies
bodies = getBodyData(mbinPath);

%%% Color options/schemes
colors = get_colors();

% --------------------------
% Actually load families
% --------------------------

%%% Grabbing specific POs
PO_SHalo__ns1 = [1.0078265893081151;
 0.0000000000000000;
 -0.0248694462406022;
 0.0000000000000004;
 -0.0191250403149298;
 -0.0000000000000027;
 1.7721271619371814];

PO_SHalo__ns1_ZH = [1.0076999551872323;
 0.0000000000000000;
 -0.0246753241341064;
 0.0000000000000000;
 -0.0191251649021168;
 0.0000000000000000;
 1.7572886949755335];
% ----------------------------------
PO_L2_L_1T_4P2_1P3__ns1 = [0.9831916400446440;
 0.0000000000000000;
 -0.0131282883062513;
 0.0000000000000000;
 0.0349348518195761;
 0.0000000000000000;
 8.5802195411897149];

PO_L2_L_1T_4P2_1P3__ns1_ZH = [0.9832247139071548;
 0.0000000000000000;
 -0.0131390111396222;
 0.0000000000000000;
 0.0349474926253687;
 0.0000000000000000;
 8.5961561262183626];
% ----------------------------------
PO_L2_L_1T_4P2_1P3__ns2 = [0.9828179111434627;
 0.0000000000000000;
 -0.0141515651995302;
 0.0000000000000000;
 0.0299207531685869;
 0.0000000000000000;
 8.4292195411897151];

PO_L2_L_1T_4P2_1P3__ns2_ZH = [0.9828124891692662;
 0.0000000000000000;
 -0.0142857305887770;
 0.0000000000000000;
 0.0292061374680624;
 0.0000000000000000;
 8.4217024872517037];
% ----------------------------------
PO_DRO__s1 = [1.0023018142036717;
 0.0000000000000000;
 -0.0000000000000000;
 0.0000000000000000;
 -0.1066018066574531;
 0.0000000000000000;
 0.1373228169374367];

PO_DRO__s1_ZH = [1.0023015545697822;
 0.0000000000000000;
 0.0000000000000000;
 0.0000000000000000;
 -0.1069490020675646;
 0.0000000000000000;
 0.1384168861131381];
% ----------------------------------
PO_DPO__s1 = [1.0171603993238323;
 0.0000000000000000;
 -0.0000000000000000;
 0.0000000000000000;
 0.0122593337884390;
 0.0000000000000000;
 2.2659437264148679];

PO_DPO__s1_ZH = [1.0171155971613994;
 0.0000000000000000;
 0.0000000000000000;
 0.0000000000000000;
 0.0121028330728206;
 0.0000000000000000;
 2.2364455212280565];
% ----------------------------------
PO_DPO_2P3__ns1 = [1.0167129648757285;
 0.0000000000000000;
 0.0000000000000000;
 0.0000000000000000;
 0.0133577732382800;
 0.0000000000000000;
 6.4937848037059496];

PO_DPO_2P3__ns1_ZH = [1.0166796432704273;
 0.0000000000000000;
 0.0000000000000000;
 0.0000000000000000;
 0.0131816849419036;
 0.0000000000000000;
 6.4163889593664960];
% ----------------------------------
PO_DPO_2P4__ns1 = [1.0173707164079508;
 0.0000000000000000;
 0.0000000000000000;
 0.0000000000000000;
 0.0110944350116362;
 0.0000000000000000;
 8.0510016147908576];

PO_DPO_2P4__ns1_ZH = [1.0173542848831678;
 0.0000000000000000;
 0.0000000000000000;
 0.0000000000000000;
 0.0107493073574068;
 0.0000000000000000;
 7.8732025856505565];
% ----------------------------------
PO_Hg2__ns1 = [1.0192145642709021;
 0.0000000000000000;
 0.0000000000000000;
 0.0000000000004866;
 -0.0021050574269589;
 0.0000000000000000;
 8.7727347725916509];

PO_Hg2__ns1_ZH = [1.0195240977789048;
 0.0000000000000000;
 0.0000000000000000;
 0.0000000000000000;
 -0.0039626379971638;
 0.0000000000000000;
 8.4941025786407423];
% ----------------------------------
PO_Hg2_2P2__ns1 = [0.9833728629023240;
 0.0000000000000000;
 -0.0078930029006791;
 0.0000000000000000;
 0.0022065217289831;
 0.0000000000000000;
 14.0364627012535141];

PO_Hg2_2P2__ns1_ZH = [0.9835692394907115;
 0.0000000000000000;
 -0.0078419165999285;
 0.0000000000000000;
 0.0013821728404153;
 0.0000000000000000;
 14.0196704496260747];
% ----------------------------------
% % % % PO_Hg2_2T__ns1 
% % % % 
% % % % PO_Hg2_2T__ns1_ZH = 
% ----------------------------------
% % % % PO_Hg2_3T__ns1 
% % % % 
% % % % PO_Hg2_3T__ns1_ZH = 
% ----------------------------------
PO_Hg2_3T__ns2 = [1.0015141007363810;
 0.0000000000000000;
 0.0049107290403815;
 0.0000000000000000;
 0.0844775451517511;
 0.0000000000000000;
 8.9957123336049634];

PO_Hg2_3T__ns2_ZH = [1.0017094313883201;
 0.0000000000000000;
 0.0050808009212432;
 0.0000000000000000;
 0.0819542487576313;
 0.0000000000000000;
 9.0705207080894628];
% ----------------------------------
PO_Se7__s1 = [1.0157804758355187;
 0.0000000000000000;
 0.0000000000000000;
 0.0000000000000000;
 0.0078109160830706;
 0.0000000000000000;
 9.0527283586993299];

PO_Se7__s1_ZH = [1.0163550867348314;
 0.0000000000000000;
 0.0000000000000000;
 0.0000000000000000;
 0.0060707779930108;
 0.0000000000000000;
 9.2623025095436482];
% ----------------------------------
PO_Se7__ns1 = [1.0172567504140277;
 0.0000000000000000;
 -0.0000000000000000;
 0.0000000000000000;
 0.0057715970121136;
 0.0000000000000000;
 10.6118252000616664];

PO_Se7__ns1_ZH = [1.0189463385307829;
 0.0000000000000000;
 -0.0000000000000000;
 0.0000000000000000;
 0.0021242260174039;
 0.0000000000000000;
 11.6367562910424418];
% ----------------------------------
% % % % PO_Se7_6P2__ns1 
% % % % 
% % % % PO_Se7_6P2__ns1_ZH = 
% ----------------------------------
PO_LoPO__s1 = [0.9830705789651824;
 0.0000000000000000;
 0.0000000000000000;
 0.0000000000000000;
 -0.0125129914442458;
 0.0000000000000000;
 2.1570155042999755];

PO_LoPO__s1_ZH = [0.9830172629743472;
 0.0000000000000000;
 0.0000000000000000;
 0.0000000000000000;
 -0.0126073846535726;
 0.0000000000000000;
 2.1864224844101146];
% ----------------------------------
PO_LoPO_2P2_1T__ns1 = [0.9880134301736027;
 0.0000000000000000;
 -0.0168339215530371;
 0.0000000000000000;
 -0.0040733760074276;
 0.0000000000000000;
 4.2714646720470677];

PO_LoPO_2P2_1T__ns1_ZH = [0.9881008483546246;
 0.0000000000000000;
 -0.0166893988130436;
 0.0000000000000000;
 -0.0054083136887970;
 0.0000000000000000;
 4.3226982284263702];
% ----------------------------------
PO_LoPO_2P4__ns1 = [0.9844892944067551;
 0.0000000000000000;
 0.0000000000000000;
 0.0000000000000000;
 -0.0167432916869640;
 0.0000000000000000;
 7.8983179359852862];

PO_LoPO_2P4__ns1_ZH = [0.9842665260679727;
 0.0000000000000000;
 0.0000000000000000;
 0.0000000000000000;
 -0.0162524623480647;
 0.0000000000000000;
 8.0772539037724957];

% ========================================================================
%%% Integrate and plot the POs
% ========================================================================
% -------------------------------------------------
%%% Set up parameters
% -------------------------------------------------
%%% Set primary and secondary bodies
[primary, secondary] = assignPrimaryAndSecondary_CR3BP('Jupiter_Europa.CR3BP_J2pJ4pJ6pJ2s', bodies);

%%% Normalizing constants
[rNorm, tNorm, vNorm] = cr3bp_norms(primary, secondary, bodies.constants.G);
clear tNorm vNorm

prms.u    = secondary.MR;
prms.R2 = secondary.R_n;
prms.n  = 1;
%%% Collinear equillibrium points
rLPs_n = EquilibriumPoints(prms.u, prms.n);
%%% Set integrator handle
integratorHandle = @Int_CR3BnSTM;
    

%%% prms for integration
prms_ZH.u    = secondary.MR;
prms_ZH.R2 = secondary.R_n;
prms_ZH.J2p = primary.J2;
prms_ZH.J4p = primary.J4;
prms_ZH.J6p = primary.J6;
prms_ZH.J2s = secondary.J2;

prms_ZH.R1 = primary.R / rNorm;

%%% Determine the mean motion via the ephemeris method
tN_ephemeris = sqrt((secondary.a^3) / (bodies.constants.G*(primary.mass+secondary.mass)));
vN_ephemeris = rNorm / tN_ephemeris;

prms_ZH.n = secondary.meanMot*tN_ephemeris;

%%% Collinear equillibrium points
rLPs_n_ZH = collinearEquilibriumPoints_ZH(prms_ZH);

%%% Set integrator handle
integratorHandle_ZH = @Int_CR3BnSTM_J2pJ4pJ6pJ2s;


%%% Integration options
tol = 1e-13;
options = odeset('RelTol',tol,'AbsTol',tol);
options_apsis = odeset('Events', @event_Apsis_CR3BP, 'RelTol',tol,'AbsTol',tol);


%%% Initialize STM as column vector to be added with state
stm0_colVec = reshape(eye(6),36,1);
        

%% ------------------------------------------------
%%% Southern Halo, ns1 .... (L2_L_1T  ns1)
% -------------------------------------------------
if run_L2_L_1T_ns1
fprintf('---------------------------------')
fprintf('L2_L_1T_ns1\n')

[T_SHalo__ns1, X_SHalo__ns1, t_ev, x_ev, in_ev] = ode113(@Int_CR3BnSTM, [0, PO_SHalo__ns1(end)], [PO_SHalo__ns1(1:6); stm0_colVec], options_apsis, prms);
apses_alts = rowNorm(x_ev(:,1:3) - [1-prms.u,0,0]) - prms.R2;
% figure; hold all
% plot3(X_SHalo__ns1(:,1),X_SHalo__ns1(:,2),X_SHalo__ns1(:,3), 'linewidth', 2, 'color', colors.black)
% plotSecondary(secondary)
% PlotBoi3_CR3Bn(26)
% axis equal
% title('SHalo ... ns1','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)
% 
% plot3(x_ev((apses_alts < 5/rNorm),1), x_ev((apses_alts < 5/rNorm),2), x_ev((apses_alts < 5/rNorm),3), '.','markersize',20, 'color', colors.red2)
% set(gca,'yticklabel',num2str(get(gca,'ytick')','%1.3f'))
% set(gca,'zticklabel',num2str(get(gca,'ztick')','%1.3f'))


[T_SHalo__ns1_ZH, X_SHalo__ns1_ZH, t_ev_ZH, x_ev_ZH, in_ev_ZH] = ode113(@Int_CR3BnSTM_J2pJ4pJ6pJ2s, [0, PO_SHalo__ns1_ZH(end)], [PO_SHalo__ns1_ZH(1:6); stm0_colVec], options_apsis, prms_ZH);
apses_alts_ZH = rowNorm(x_ev_ZH(:,1:3) - [1-prms_ZH.u,0,0]) - prms_ZH.R2;
fprintf('Prop Error = %1.2e\n', norm(X_SHalo__ns1_ZH(end,1:6) - X_SHalo__ns1_ZH(1,1:6)) / norm(X_SHalo__ns1_ZH(1,1:6)))
% figure; hold all
% plot3(X_SHalo__ns1_ZH(:,1),X_SHalo__ns1_ZH(:,2),X_SHalo__ns1_ZH(:,3), 'linewidth', 2, 'color', colors.black)
% plotSecondary(secondary)
% PlotBoi3_CR3Bn(26)
% axis equal
% title('SHalo ... ns1 ... ZH','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)
% 
% plot3(x_ev_ZH((apses_alts_ZH < 5/rNorm),1), x_ev_ZH((apses_alts_ZH < 5/rNorm),2), x_ev_ZH((apses_alts_ZH < 5/rNorm),3), '.','markersize',20, 'color', colors.red2)
% set(gca,'yticklabel',num2str(get(gca,'ytick')','%1.3f'))
% set(gca,'zticklabel',num2str(get(gca,'ztick')','%1.3f'))


% figure; hold all
% p_3B = plot3(X_SHalo__ns1(:,1),X_SHalo__ns1(:,2),X_SHalo__ns1(:,3), 'linewidth', 2, 'color', colors.blue2);
% p_ZH = plot3(X_SHalo__ns1_ZH(:,1),X_SHalo__ns1_ZH(:,2),X_SHalo__ns1_ZH(:,3), 'linewidth', 2, 'color', colors.red2);
% plotSecondary(secondary)
% PlotBoi3_CR3Bn(26)
% axis equal
% legend([p_3B, p_ZH], 'CR3BP', 'CR3BP w ZH')
% set(gca,'yticklabel',num2str(get(gca,'ytick')','%1.3f'))
% set(gca,'zticklabel',num2str(get(gca,'ztick')','%1.3f'))


stm_tf_t0                           = reshape(X_SHalo__ns1(end,7:42),6,6);
monodromy                           = stm_tf_t0;
[eigenVectors_new, eigenValues_new] = eig(monodromy);
[S1, S2]                            = getStabilityIndices(diag(eigenValues_new));

stm_tf_t0_ZH                              = reshape(X_SHalo__ns1_ZH(end,7:42),6,6);
monodromy_ZH                              = stm_tf_t0_ZH;
[eigenVectors_new_ZH, eigenValues_new_ZH] = eig(monodromy_ZH);
[S1_ZH, S2_ZH]                            = getStabilityIndices(diag(eigenValues_new_ZH));

fprintf('[S1, S2]       = [%1.3f, %1.3f]\n', S1, S2)
fprintf('[S1_ZH, S2_ZH] = [%1.3f, %1.3f]\n', S1_ZH, S2_ZH)


latLons_SHalo__ns1 = zeros(length(T_SHalo__ns1),2);
for kk = 1:length(T_SHalo__ns1)
    [lat_deg, lon_deg] = BCR2latlon(X_SHalo__ns1(kk,1:3)', 'secondary', prms.u);
    latLons_SHalo__ns1(kk,:) = [lat_deg, lon_deg];
end
latLons_apses_SHalo__ns1 = zeros(length(t_ev),2);
for kk = 1:length(t_ev)
    [lat_deg, lon_deg] = BCR2latlon(x_ev(kk,1:3)', 'secondary', prms.u);
    latLons_apses_SHalo__ns1(kk,:) = [lat_deg, lon_deg];
end


latLons_SHalo__ns1_ZH = zeros(length(T_SHalo__ns1_ZH),2);
for kk = 1:length(T_SHalo__ns1_ZH)
    [lat_deg, lon_deg] = BCR2latlon(X_SHalo__ns1_ZH(kk,1:3)', 'secondary', prms_ZH.u);
    latLons_SHalo__ns1_ZH(kk,:) = [lat_deg, lon_deg];
end
latLons_apses_SHalo__ns1_ZH = zeros(length(t_ev_ZH),2);
for kk = 1:length(t_ev_ZH)
    [lat_deg, lon_deg] = BCR2latlon(x_ev_ZH(kk,1:3)', 'secondary', prms_ZH.u);
    latLons_apses_SHalo__ns1_ZH(kk,:) = [lat_deg, lon_deg];
end


% figure; hold all
% xlim([0 360])
% ylim([-90 90])
% PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')
% h = image(xlim, -ylim, bodies.europa.img);
% [lons_new] = convert_lon180_to_lon360(latLons_SHalo__ns1(:,2));
% plot(lons_new,latLons_SHalo__ns1(:,1),'k.')
% [lons_new_impact] = convert_lon180_to_lon360(latLons_apses_SHalo__ns1((apses_alts < 5/rNorm),2));
% p_landing = plot(lons_new_impact, latLons_apses_SHalo__ns1((apses_alts < 5/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.drkred, 'markerfacecolor', colors.red2);
% legend([p_landing], 'Tangent Impact', 'fontsize', 20, 'fontname', 'Times New Roman', 'location', 'northeast')
% title('SHalo ... ns1','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)
% 
% 
figure; hold all
xlim([0 360])
ylim([-90 90])
PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')
h = image(xlim, -ylim, bodies.europa.img);
[lons_new_ZH] = convert_lon180_to_lon360(latLons_SHalo__ns1_ZH(:,2));
plot(lons_new_ZH,latLons_SHalo__ns1_ZH(:,1),'.', 'color', colors.red2)
[lons_new_impact_ZH] = convert_lon180_to_lon360(latLons_apses_SHalo__ns1_ZH((apses_alts_ZH < 5/rNorm),2));
p_landing_ZH = plot(lons_new_impact_ZH, latLons_apses_SHalo__ns1_ZH((apses_alts_ZH < 5/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.black, 'markerfacecolor', colors.drkgrey);
legend([p_landing_ZH], 'Tangent Impact', 'fontsize', 20, 'fontname', 'Times New Roman', 'location', 'southeast')
title('SHalo ... ns1 ... ZH','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)
% 
% figure; hold all
% xlim([0 360])
% ylim([-90 90])
% PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')
% h = image(xlim, -ylim, bodies.europa.img);
% plot(lons_new,latLons_SHalo__ns1(:,1),'.','color',colors.blue2)
% plot(lons_new_ZH,latLons_SHalo__ns1_ZH(:,1),'.', 'color', colors.red2)
% p_landing = plot(lons_new_impact, latLons_apses_SHalo__ns1((apses_alts < 5/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.drkblue, 'markerfacecolor', colors.blue2);
% p_landing_ZH = plot(lons_new_impact_ZH, latLons_apses_SHalo__ns1_ZH((apses_alts_ZH < 5/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.drkred, 'markerfacecolor', colors.red2);
% legend([p_landing, p_landing_ZH], 'Tangent Impact - CR3BP', 'Tangent Impact - CR3BP w ZH', 'fontsize', 20, 'fontname', 'Times New Roman', 'location', 'northeast')
% title('SHalo ... ns1 ... Both','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)

figure; hold all
% plot3(X_SHalo__ns1_ZH(:,1),X_SHalo__ns1_ZH(:,2),X_SHalo__ns1_ZH(:,3), 'linewidth', 2, 'color', colors.red2)
plot3(X_SHalo__ns1_ZH(:,1),X_SHalo__ns1_ZH(:,2),X_SHalo__ns1_ZH(:,3), 'linewidth', 2, 'color', colors.mag)
plotTrajShadows(X_SHalo__ns1_ZH, 2, colors.grey, 'x', 1.015, 'y', 1.6e-2, 'z', -26e-3, 'bodyshadow', [1-prms_ZH.u, prms_ZH.R2])
plotSecondary(secondary)
PlotBoi3_CR3Bn(26)
axis equal
plot3(x_ev_ZH((apses_alts_ZH < 5/rNorm),1), x_ev_ZH((apses_alts_ZH < 5/rNorm),2), x_ev_ZH((apses_alts_ZH < 5/rNorm),3), '.','markersize',20, 'color', colors.black)
set(gca,'yticklabel',num2str(get(gca,'ytick')','%1.3f'))
set(gca,'zticklabel',num2str(get(gca,'ztick')','%1.3f'))



end % run_L2_L_1T
%% ------------------------------------------------
%%% L2_L_1T_4P2_1P3__ns1
% -------------------------------------------------
if run_L2_L_1T_4P2_1P3__ns1
fprintf('---------------------------------')
fprintf('L2_L_1T_4P2_1P3__ns1\n')

[T_L2_L_1T_4P2_1P3__ns1, X_L2_L_1T_4P2_1P3__ns1, t_ev, x_ev, in_ev] = ode113(@Int_CR3BnSTM, [0, PO_L2_L_1T_4P2_1P3__ns1(end)], [PO_L2_L_1T_4P2_1P3__ns1(1:6); stm0_colVec], options_apsis, prms);
apses_alts = rowNorm(x_ev(:,1:3) - [1-prms.u,0,0]) - prms.R2;
% figure; hold all
% plot3(X_L2_L_1T_4P2_1P3__ns1(:,1),X_L2_L_1T_4P2_1P3__ns1(:,2),X_L2_L_1T_4P2_1P3__ns1(:,3), 'linewidth', 2, 'color', colors.black)
% plotSecondary(secondary)
% PlotBoi3_CR3Bn(26)
% axis equal
% title('L2\_L\_1T\_4P2\_1P3 ... ns1','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)
% 
% plot3(x_ev((apses_alts < 5/rNorm),1), x_ev((apses_alts < 5/rNorm),2), x_ev((apses_alts < 5/rNorm),3), '.','markersize',20, 'color', colors.red2)
% set(gca,'yticklabel',num2str(get(gca,'ytick')','%1.3f'))
% set(gca,'zticklabel',num2str(get(gca,'ztick')','%1.3f'))


[T_L2_L_1T_4P2_1P3__ns1_ZH, X_L2_L_1T_4P2_1P3__ns1_ZH, t_ev_ZH, x_ev_ZH, in_ev_ZH] = ode113(@Int_CR3BnSTM_J2pJ4pJ6pJ2s, [0, PO_L2_L_1T_4P2_1P3__ns1_ZH(end)], [PO_L2_L_1T_4P2_1P3__ns1_ZH(1:6); stm0_colVec], options_apsis, prms_ZH);
apses_alts_ZH = rowNorm(x_ev_ZH(:,1:3) - [1-prms_ZH.u,0,0]) - prms_ZH.R2;
fprintf('Prop Error = %1.2e\n', norm(X_L2_L_1T_4P2_1P3__ns1_ZH(end,1:6) - X_L2_L_1T_4P2_1P3__ns1_ZH(1,1:6)) / norm(X_L2_L_1T_4P2_1P3__ns1_ZH(1,1:6)))
% figure; hold all
% plot3(X_L2_L_1T_4P2_1P3__ns1_ZH(:,1),X_L2_L_1T_4P2_1P3__ns1_ZH(:,2),X_L2_L_1T_4P2_1P3__ns1_ZH(:,3), 'linewidth', 2, 'color', colors.black)
% plotSecondary(secondary)
% PlotBoi3_CR3Bn(26)
% axis equal
% title('L2\_L\_1T\_4P2\_1P3 ... ns1 ... ZH','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)
% 
% plot3(x_ev_ZH((apses_alts_ZH < 5/rNorm),1), x_ev_ZH((apses_alts_ZH < 5/rNorm),2), x_ev_ZH((apses_alts_ZH < 5/rNorm),3), '.','markersize',20, 'color', colors.red2)
% set(gca,'yticklabel',num2str(get(gca,'ytick')','%1.3f'))
% set(gca,'zticklabel',num2str(get(gca,'ztick')','%1.3f'))
% 
% 
% figure; hold all
% p_3B = plot3(X_L2_L_1T_4P2_1P3__ns1(:,1),X_L2_L_1T_4P2_1P3__ns1(:,2),X_L2_L_1T_4P2_1P3__ns1(:,3), 'linewidth', 2, 'color', colors.blue2);
% p_ZH = plot3(X_L2_L_1T_4P2_1P3__ns1_ZH(:,1),X_L2_L_1T_4P2_1P3__ns1_ZH(:,2),X_L2_L_1T_4P2_1P3__ns1_ZH(:,3), 'linewidth', 2, 'color', colors.red2);
% plotSecondary(secondary)
% PlotBoi3_CR3Bn(26)
% axis equal
% legend([p_3B, p_ZH], 'CR3BP', 'CR3BP w ZH')
% set(gca,'yticklabel',num2str(get(gca,'ytick')','%1.3f'))
% set(gca,'zticklabel',num2str(get(gca,'ztick')','%1.3f'))


stm_tf_t0                           = reshape(X_L2_L_1T_4P2_1P3__ns1(end,7:42),6,6);
monodromy                           = stm_tf_t0;
[eigenVectors_new, eigenValues_new] = eig(monodromy);
[S1, S2]                            = getStabilityIndices(diag(eigenValues_new));

stm_tf_t0_ZH                              = reshape(X_L2_L_1T_4P2_1P3__ns1_ZH(end,7:42),6,6);
monodromy_ZH                              = stm_tf_t0_ZH;
[eigenVectors_new_ZH, eigenValues_new_ZH] = eig(monodromy_ZH);
[S1_ZH, S2_ZH]                            = getStabilityIndices(diag(eigenValues_new_ZH));

fprintf('[S1, S2]       = [%1.3f, %1.3f]\n', S1, S2)
fprintf('[S1_ZH, S2_ZH] = [%1.3f, %1.3f]\n', S1_ZH, S2_ZH)

latLons_L2_L_1T_4P2_1P3__ns1 = zeros(length(T_L2_L_1T_4P2_1P3__ns1),2);
for kk = 1:length(T_L2_L_1T_4P2_1P3__ns1)
    [lat_deg, lon_deg] = BCR2latlon(X_L2_L_1T_4P2_1P3__ns1(kk,1:3)', 'secondary', prms.u);
    latLons_L2_L_1T_4P2_1P3__ns1(kk,:) = [lat_deg, lon_deg];
end
latLons_apses_L2_L_1T_4P2_1P3__ns1 = zeros(length(t_ev),2);
for kk = 1:length(t_ev)
    [lat_deg, lon_deg] = BCR2latlon(x_ev(kk,1:3)', 'secondary', prms.u);
    latLons_apses_L2_L_1T_4P2_1P3__ns1(kk,:) = [lat_deg, lon_deg];
end


latLons_L2_L_1T_4P2_1P3__ns1_ZH = zeros(length(T_L2_L_1T_4P2_1P3__ns1_ZH),2);
for kk = 1:length(T_L2_L_1T_4P2_1P3__ns1_ZH)
    [lat_deg, lon_deg] = BCR2latlon(X_L2_L_1T_4P2_1P3__ns1_ZH(kk,1:3)', 'secondary', prms_ZH.u);
    latLons_L2_L_1T_4P2_1P3__ns1_ZH(kk,:) = [lat_deg, lon_deg];
end
latLons_apses_L2_L_1T_4P2_1P3__ns1_ZH = zeros(length(t_ev_ZH),2);
for kk = 1:length(t_ev_ZH)
    [lat_deg, lon_deg] = BCR2latlon(x_ev_ZH(kk,1:3)', 'secondary', prms_ZH.u);
    latLons_apses_L2_L_1T_4P2_1P3__ns1_ZH(kk,:) = [lat_deg, lon_deg];
end

% figure; hold all
% xlim([0 360])
% ylim([-90 90])
% PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')
% h = image(xlim, -ylim, bodies.europa.img);
% [lons_new] = convert_lon180_to_lon360(latLons_L2_L_1T_4P2_1P3__ns1(:,2));
% plot(lons_new,latLons_L2_L_1T_4P2_1P3__ns1(:,1),'k.')
% [lons_new_impact] = convert_lon180_to_lon360(latLons_apses_L2_L_1T_4P2_1P3__ns1((apses_alts < 5/rNorm),2));
% p_landing = plot(lons_new_impact, latLons_apses_L2_L_1T_4P2_1P3__ns1((apses_alts < 5/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.drkred, 'markerfacecolor', colors.red2);
% legend([p_landing], 'Tangent Impact', 'fontsize', 20, 'fontname', 'Times New Roman', 'location', 'northeast')
% title('L2\_L\_1T\_4P2\_1P3 ... ns1','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)


figure; hold all
xlim([0 360])
ylim([-90 90])
PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')
h = image(xlim, -ylim, bodies.europa.img);
[lons_new_ZH] = convert_lon180_to_lon360(latLons_L2_L_1T_4P2_1P3__ns1_ZH(:,2));
plot(lons_new_ZH,latLons_L2_L_1T_4P2_1P3__ns1_ZH(:,1),'.', 'color', colors.red2)
[lons_new_impact_ZH] = convert_lon180_to_lon360(latLons_apses_L2_L_1T_4P2_1P3__ns1_ZH((apses_alts_ZH < 5/rNorm),2));
p_landing_ZH = plot(lons_new_impact_ZH, latLons_apses_L2_L_1T_4P2_1P3__ns1_ZH((apses_alts_ZH < 5/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.black, 'markerfacecolor', colors.drkgrey);
legend([p_landing_ZH], 'Tangent Impact', 'fontsize', 20, 'fontname', 'Times New Roman', 'location', 'northeast')
title('L2\_L\_1T\_4P2\_1P3 ... ns1 ... ZH','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)

% figure; hold all
% xlim([0 360])
% ylim([-90 90])
% PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')
% h = image(xlim, -ylim, bodies.europa.img);
% plot(lons_new,latLons_L2_L_1T_4P2_1P3__ns1(:,1),'.','color',colors.blue2)
% plot(lons_new_ZH,latLons_L2_L_1T_4P2_1P3__ns1_ZH(:,1),'.', 'color', colors.red2)
% p_landing = plot(lons_new_impact, latLons_apses_L2_L_1T_4P2_1P3__ns1((apses_alts < 5/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.drkblue, 'markerfacecolor', colors.blue2);
% p_landing_ZH = plot(lons_new_impact_ZH, latLons_apses_L2_L_1T_4P2_1P3__ns1_ZH((apses_alts_ZH < 5/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.drkred, 'markerfacecolor', colors.red2);
% legend([p_landing, p_landing_ZH], 'Tangent Impact - CR3BP', 'Tangent Impact - CR3BP w ZH', 'fontsize', 20, 'fontname', 'Times New Roman', 'location', 'northeast')
% title('L2\_L\_1T\_4P2\_1P3 ... ns1 ... Both','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)

figure; hold all
% plot3(X_L2_L_1T_4P2_1P3__ns1_ZH(:,1),X_L2_L_1T_4P2_1P3__ns1_ZH(:,2),X_L2_L_1T_4P2_1P3__ns1_ZH(:,3), 'linewidth', 2, 'color', colors.red2)
plot3(X_L2_L_1T_4P2_1P3__ns1_ZH(:,1),X_L2_L_1T_4P2_1P3__ns1_ZH(:,2),X_L2_L_1T_4P2_1P3__ns1_ZH(:,3), 'linewidth', 2, 'color', colors.mag)
plotTrajShadows(X_L2_L_1T_4P2_1P3__ns1_ZH, 2, colors.grey, 'x', 1.025, 'y', 2.9e-2, 'z', -2.8e-2, 'bodyshadow', [1-prms_ZH.u, prms_ZH.R2])
plotSecondary(secondary)
PlotBoi3_CR3Bn(26)
axis equal
plot3(x_ev_ZH((apses_alts_ZH < 5/rNorm),1), x_ev_ZH((apses_alts_ZH < 5/rNorm),2), x_ev_ZH((apses_alts_ZH < 5/rNorm),3), '.','markersize',20, 'color', colors.black)
set(gca,'yticklabel',num2str(get(gca,'ytick')','%1.3f'))
set(gca,'zticklabel',num2str(get(gca,'ztick')','%1.3f'))



end % run_L2_L_1T_4P2_1P3__ns1
%% ------------------------------------------------
%%% L2_L_1T_4P2_1P3__ns2
% -------------------------------------------------
if run_L2_L_1T_4P2_1P3__ns2
fprintf('---------------------------------')
fprintf('L2_L_1T_4P2_1P3__ns2\n')

[T_L2_L_1T_4P2_1P3__ns2, X_L2_L_1T_4P2_1P3__ns2, t_ev, x_ev, in_ev] = ode113(@Int_CR3BnSTM, [0, PO_L2_L_1T_4P2_1P3__ns2(end)], [PO_L2_L_1T_4P2_1P3__ns2(1:6); stm0_colVec], options_apsis, prms);
apses_alts = rowNorm(x_ev(:,1:3) - [1-prms.u,0,0]) - prms.R2;
% figure; hold all
% plot3(X_L2_L_1T_4P2_1P3__ns2(:,1),X_L2_L_1T_4P2_1P3__ns2(:,2),X_L2_L_1T_4P2_1P3__ns2(:,3), 'linewidth', 2, 'color', colors.black)
% plotSecondary(secondary)
% PlotBoi3_CR3Bn(26)
% axis equal
% title('L2\_L\_1T\_4P2\_1P3 ... ns2','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)
% 
% plot3(x_ev((apses_alts < 5/rNorm),1), x_ev((apses_alts < 5/rNorm),2), x_ev((apses_alts < 5/rNorm),3), '.','markersize',20, 'color', colors.red2)
% set(gca,'yticklabel',num2str(get(gca,'ytick')','%1.3f'))
% set(gca,'zticklabel',num2str(get(gca,'ztick')','%1.3f'))


[T_L2_L_1T_4P2_1P3__ns2_ZH, X_L2_L_1T_4P2_1P3__ns2_ZH, t_ev_ZH, x_ev_ZH, in_ev_ZH] = ode113(@Int_CR3BnSTM_J2pJ4pJ6pJ2s, [0, PO_L2_L_1T_4P2_1P3__ns2_ZH(end)], [PO_L2_L_1T_4P2_1P3__ns2_ZH(1:6); stm0_colVec], options_apsis, prms_ZH);
apses_alts_ZH = rowNorm(x_ev_ZH(:,1:3) - [1-prms_ZH.u,0,0]) - prms_ZH.R2;
fprintf('Prop Error = %1.2e\n', norm(X_L2_L_1T_4P2_1P3__ns2_ZH(end,1:6) - X_L2_L_1T_4P2_1P3__ns2_ZH(1,1:6)) / norm(X_L2_L_1T_4P2_1P3__ns2_ZH(1,1:6)))
% figure; hold all
% plot3(X_L2_L_1T_4P2_1P3__ns2_ZH(:,1),X_L2_L_1T_4P2_1P3__ns2_ZH(:,2),X_L2_L_1T_4P2_1P3__ns2_ZH(:,3), 'linewidth', 2, 'color', colors.black)
% plotSecondary(secondary)
% PlotBoi3_CR3Bn(26)
% axis equal
% title('L2\_L\_1T\_4P2\_1P3 ... ns2 ... ZH','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)
% 
% plot3(x_ev_ZH((apses_alts_ZH < 5/rNorm),1), x_ev_ZH((apses_alts_ZH < 5/rNorm),2), x_ev_ZH((apses_alts_ZH < 5/rNorm),3), '.','markersize',20, 'color', colors.red2)
% set(gca,'yticklabel',num2str(get(gca,'ytick')','%1.3f'))
% set(gca,'zticklabel',num2str(get(gca,'ztick')','%1.3f'))


% figure; hold all
% p_3B = plot3(X_L2_L_1T_4P2_1P3__ns2(:,1),X_L2_L_1T_4P2_1P3__ns2(:,2),X_L2_L_1T_4P2_1P3__ns2(:,3), 'linewidth', 2, 'color', colors.blue2);
% p_ZH = plot3(X_L2_L_1T_4P2_1P3__ns2_ZH(:,1),X_L2_L_1T_4P2_1P3__ns2_ZH(:,2),X_L2_L_1T_4P2_1P3__ns2_ZH(:,3), 'linewidth', 2, 'color', colors.red2);
% plotSecondary(secondary)
% PlotBoi3_CR3Bn(26)
% axis equal
% legend([p_3B, p_ZH], 'CR3BP', 'CR3BP w ZH')
% set(gca,'yticklabel',num2str(get(gca,'ytick')','%1.3f'))
% set(gca,'zticklabel',num2str(get(gca,'ztick')','%1.3f'))


stm_tf_t0                           = reshape(X_L2_L_1T_4P2_1P3__ns2(end,7:42),6,6);
monodromy                           = stm_tf_t0;
[eigenVectors_new, eigenValues_new] = eig(monodromy);
[S1, S2]                            = getStabilityIndices(diag(eigenValues_new));

stm_tf_t0_ZH                              = reshape(X_L2_L_1T_4P2_1P3__ns2_ZH(end,7:42),6,6);
monodromy_ZH                              = stm_tf_t0_ZH;
[eigenVectors_new_ZH, eigenValues_new_ZH] = eig(monodromy_ZH);
[S1_ZH, S2_ZH]                            = getStabilityIndices(diag(eigenValues_new_ZH));

fprintf('[S1, S2]       = [%1.3f, %1.3f]\n', S1, S2)
fprintf('[S1_ZH, S2_ZH] = [%1.3f, %1.3f]\n', S1_ZH, S2_ZH)


latLons_L2_L_1T_4P2_1P3__ns2 = zeros(length(T_L2_L_1T_4P2_1P3__ns2),2);
for kk = 1:length(T_L2_L_1T_4P2_1P3__ns2)
    [lat_deg, lon_deg] = BCR2latlon(X_L2_L_1T_4P2_1P3__ns2(kk,1:3)', 'secondary', prms.u);
    latLons_L2_L_1T_4P2_1P3__ns2(kk,:) = [lat_deg, lon_deg];
end
latLons_apses_L2_L_1T_4P2_1P3__ns2 = zeros(length(t_ev),2);
for kk = 1:length(t_ev)
    [lat_deg, lon_deg] = BCR2latlon(x_ev(kk,1:3)', 'secondary', prms.u);
    latLons_apses_L2_L_1T_4P2_1P3__ns2(kk,:) = [lat_deg, lon_deg];
end


latLons_L2_L_1T_4P2_1P3__ns2_ZH = zeros(length(T_L2_L_1T_4P2_1P3__ns2_ZH),2);
for kk = 1:length(T_L2_L_1T_4P2_1P3__ns2_ZH)
    [lat_deg, lon_deg] = BCR2latlon(X_L2_L_1T_4P2_1P3__ns2_ZH(kk,1:3)', 'secondary', prms_ZH.u);
    latLons_L2_L_1T_4P2_1P3__ns2_ZH(kk,:) = [lat_deg, lon_deg];
end
latLons_apses_L2_L_1T_4P2_1P3__ns2_ZH = zeros(length(t_ev_ZH),2);
for kk = 1:length(t_ev_ZH)
    [lat_deg, lon_deg] = BCR2latlon(x_ev_ZH(kk,1:3)', 'secondary', prms_ZH.u);
    latLons_apses_L2_L_1T_4P2_1P3__ns2_ZH(kk,:) = [lat_deg, lon_deg];
end

% figure; hold all
% xlim([0 360])
% ylim([-90 90])
% PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')
% h = image(xlim, -ylim, bodies.europa.img);
% [lons_new] = convert_lon180_to_lon360(latLons_L2_L_1T_4P2_1P3__ns2(:,2));
% plot(lons_new,latLons_L2_L_1T_4P2_1P3__ns2(:,1),'k.')
% [lons_new_impact] = convert_lon180_to_lon360(latLons_apses_L2_L_1T_4P2_1P3__ns2((apses_alts < 5/rNorm),2));
% p_landing = plot(lons_new_impact, latLons_apses_L2_L_1T_4P2_1P3__ns2((apses_alts < 5/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.drkred, 'markerfacecolor', colors.red2);
% legend([p_landing], 'Tangent Impact', 'fontsize', 20, 'fontname', 'Times New Roman', 'location', 'northeast')
% title('L2\_L\_1T\_4P2\_1P3 ... ns2','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)


figure; hold all
xlim([0 360])
ylim([-90 90])
PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')
h = image(xlim, -ylim, bodies.europa.img);
[lons_new_ZH] = convert_lon180_to_lon360(latLons_L2_L_1T_4P2_1P3__ns2_ZH(:,2));
plot(lons_new_ZH,latLons_L2_L_1T_4P2_1P3__ns2_ZH(:,1),'.', 'color', colors.red2)
[lons_new_impact_ZH] = convert_lon180_to_lon360(latLons_apses_L2_L_1T_4P2_1P3__ns2_ZH((apses_alts_ZH < 5/rNorm),2));
p_landing_ZH = plot(lons_new_impact_ZH, latLons_apses_L2_L_1T_4P2_1P3__ns2_ZH((apses_alts_ZH < 5/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.black, 'markerfacecolor', colors.drkgrey);
legend([p_landing_ZH], 'Tangent Impact', 'fontsize', 20, 'fontname', 'Times New Roman', 'location', 'northeast')
title('L2\_L\_1T\_4P2\_1P3 ... ns2 ... ZH','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)

% figure; hold all
% xlim([0 360])
% ylim([-90 90])
% PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')
% h = image(xlim, -ylim, bodies.europa.img);
% plot(lons_new,latLons_L2_L_1T_4P2_1P3__ns2(:,1),'.','color',colors.blue2)
% plot(lons_new_ZH,latLons_L2_L_1T_4P2_1P3__ns2_ZH(:,1),'.', 'color', colors.red2)
% p_landing = plot(lons_new_impact, latLons_apses_L2_L_1T_4P2_1P3__ns2((apses_alts < 5/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.drkblue, 'markerfacecolor', colors.blue2);
% p_landing_ZH = plot(lons_new_impact_ZH, latLons_apses_L2_L_1T_4P2_1P3__ns2_ZH((apses_alts_ZH < 5/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.drkred, 'markerfacecolor', colors.red2);
% legend([p_landing, p_landing_ZH], 'Tangent Impact - CR3BP', 'Tangent Impact - CR3BP w ZH', 'fontsize', 20, 'fontname', 'Times New Roman', 'location', 'northeast')
% title('L2\_L\_1T\_4P2\_1P3 ... ns2 ... Both','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)


figure; hold all
% plot3(X_L2_L_1T_4P2_1P3__ns2_ZH(:,1),X_L2_L_1T_4P2_1P3__ns2_ZH(:,2),X_L2_L_1T_4P2_1P3__ns2_ZH(:,3), 'linewidth', 2, 'color', colors.red2)
plot3(X_L2_L_1T_4P2_1P3__ns2_ZH(:,1),X_L2_L_1T_4P2_1P3__ns2_ZH(:,2),X_L2_L_1T_4P2_1P3__ns2_ZH(:,3), 'linewidth', 2, 'color', colors.mag)
plotTrajShadows(X_L2_L_1T_4P2_1P3__ns2_ZH, 2, colors.grey, 'x', 1.025, 'y', 2.9e-2, 'z', -2.8e-2, 'bodyshadow', [1-prms_ZH.u, prms_ZH.R2])
plotSecondary(secondary)
PlotBoi3_CR3Bn(26)
axis equal
plot3(x_ev_ZH((apses_alts_ZH < 5/rNorm),1), x_ev_ZH((apses_alts_ZH < 5/rNorm),2), x_ev_ZH((apses_alts_ZH < 5/rNorm),3), '.','markersize',20, 'color', colors.black)
set(gca,'yticklabel',num2str(get(gca,'ytick')','%1.3f'))
set(gca,'zticklabel',num2str(get(gca,'ztick')','%1.3f'))


end % run_L2_L_1T_4P2_1P3__ns2
%% ------------------------------------------------
%%% DRO__s1
% -------------------------------------------------
if run_DRO__s1
fprintf('---------------------------------')
fprintf('DRO__s1 (planar)\n')

[T_DRO__s1, X_DRO__s1, t_ev, x_ev, in_ev] = ode113(@Int_CR3BnSTM, [0, PO_DRO__s1(end)], [PO_DRO__s1(1:6); stm0_colVec], options_apsis, prms);
apses_alts = rowNorm(x_ev(:,1:3) - [1-prms.u,0,0]) - prms.R2;
% figure; hold all
% plot3(X_DRO__s1(:,1),X_DRO__s1(:,2),X_DRO__s1(:,3), 'linewidth', 2, 'color', colors.black)
% plotSecondary(secondary)
% PlotBoi3_CR3Bn(26)
% axis equal
% title('DRO ... s1','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)
% 
% plot3(x_ev((apses_alts < 1/rNorm),1), x_ev((apses_alts < 1/rNorm),2), x_ev((apses_alts < 1/rNorm),3), '.','markersize',20, 'color', colors.red2)
% set(gca,'yticklabel',num2str(get(gca,'ytick')','%1.3f'))
% set(gca,'zticklabel',num2str(get(gca,'ztick')','%1.3f'))

[T_DRO__s1_ZH, X_DRO__s1_ZH, t_ev_ZH, x_ev_ZH, in_ev_ZH] = ode113(@Int_CR3BnSTM_J2pJ4pJ6pJ2s, [0, PO_DRO__s1_ZH(end)], [PO_DRO__s1_ZH(1:6); stm0_colVec], options_apsis, prms_ZH);
fprintf('Prop Error = %1.2e\n', norm(X_DRO__s1_ZH(end,1:6) - X_DRO__s1_ZH(1,1:6)) / norm(X_DRO__s1_ZH(1,1:6)))
figure; hold all
% plot3(X_DRO__s1_ZH(:,1),X_DRO__s1_ZH(:,2),X_DRO__s1_ZH(:,3), 'linewidth', 2, 'color', colors.red2)
plot3(X_DRO__s1_ZH(:,1),X_DRO__s1_ZH(:,2),X_DRO__s1_ZH(:,3), 'linewidth', 2, 'color', colors.mag)
plotSecondary(secondary)
PlotBoi3_CR3Bn(26)
axis equal
title('DRO ... s1 ... ZH','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)

apses_alts_ZH = rowNorm(x_ev_ZH(:,1:3) - [1-prms_ZH.u,0,0]) - prms_ZH.R2;
plot3(x_ev_ZH((apses_alts_ZH < 1/rNorm),1), x_ev_ZH((apses_alts_ZH < 1/rNorm),2), x_ev_ZH((apses_alts_ZH < 1/rNorm),3), '.','markersize',20, 'color', colors.black)
close_indices_ZH = ((18/rNorm < apses_alts_ZH) & (apses_alts_ZH < 19/rNorm));
plot3(x_ev_ZH(close_indices_ZH, 1), x_ev_ZH(close_indices_ZH, 2), x_ev_ZH(close_indices_ZH, 3), '.','markersize',20, 'color', colors.blue2)
set(gca,'yticklabel',num2str(get(gca,'ytick')','%1.3f'))
set(gca,'zticklabel',num2str(get(gca,'ztick')','%1.3f'))


view(0,90)

% 
% figure; hold all
% p_3B = plot3(X_DRO__s1(:,1),X_DRO__s1(:,2),X_DRO__s1(:,3), 'linewidth', 2, 'color', colors.blue2);
% p_ZH = plot3(X_DRO__s1_ZH(:,1),X_DRO__s1_ZH(:,2),X_DRO__s1_ZH(:,3), 'linewidth', 2, 'color', colors.red2);
% plotSecondary(secondary)
% PlotBoi3_CR3Bn(26)
% axis equal
% legend([p_3B, p_ZH], 'CR3BP', 'CR3BP w ZH')
% set(gca,'yticklabel',num2str(get(gca,'ytick')','%1.3f'))
% set(gca,'zticklabel',num2str(get(gca,'ztick')','%1.3f'))


stm_tf_t0                           = reshape(X_DRO__s1(end,7:42),6,6);
monodromy                           = stm_tf_t0;
[eigenVectors_new, eigenValues_new] = eig(monodromy);
[S1, S2]                            = getStabilityIndices(diag(eigenValues_new));

stm_tf_t0_ZH                              = reshape(X_DRO__s1_ZH(end,7:42),6,6);
monodromy_ZH                              = stm_tf_t0_ZH;
[eigenVectors_new_ZH, eigenValues_new_ZH] = eig(monodromy_ZH);
[S1_ZH, S2_ZH]                            = getStabilityIndices(diag(eigenValues_new_ZH));

fprintf('[S1, S2]       = [%1.3f, %1.3f]\n', S1, S2)
fprintf('[S1_ZH, S2_ZH] = [%1.3f, %1.3f]\n', S1_ZH, S2_ZH)


latLons_DRO__s1 = zeros(length(T_DRO__s1),2);
for kk = 1:length(T_DRO__s1)
    [lat_deg, lon_deg] = BCR2latlon(X_DRO__s1(kk,1:3)', 'secondary', prms.u);
    latLons_DRO__s1(kk,:) = [lat_deg, lon_deg];
end
latLons_apses_DRO__s1 = zeros(length(t_ev),2);
for kk = 1:length(t_ev)
    [lat_deg, lon_deg] = BCR2latlon(x_ev(kk,1:3)', 'secondary', prms.u);
    latLons_apses_DRO__s1(kk,:) = [lat_deg, lon_deg];
end


latLons_DRO__s1_ZH = zeros(length(T_DRO__s1_ZH),2);
for kk = 1:length(T_DRO__s1_ZH)
    [lat_deg, lon_deg] = BCR2latlon(X_DRO__s1_ZH(kk,1:3)', 'secondary', prms_ZH.u);
    latLons_DRO__s1_ZH(kk,:) = [lat_deg, lon_deg];
end
latLons_apses_DRO__s1_ZH = zeros(length(t_ev_ZH),2);
for kk = 1:length(t_ev_ZH)
    [lat_deg, lon_deg] = BCR2latlon(x_ev_ZH(kk,1:3)', 'secondary', prms_ZH.u);
    latLons_apses_DRO__s1_ZH(kk,:) = [lat_deg, lon_deg];
end

% figure; hold all
% xlim([0 360])
% ylim([-90 90])
% PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')
% h = image(xlim, -ylim, bodies.europa.img);
% [lons_new] = convert_lon180_to_lon360(latLons_DRO__s1(:,2));
% plot(lons_new,latLons_DRO__s1(:,1),'k.')
% [lons_new_impact] = convert_lon180_to_lon360(latLons_apses_DRO__s1((apses_alts < 1/rNorm),2));
% p_landing = plot(lons_new_impact, latLons_apses_DRO__s1((apses_alts < 1/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.drkred, 'markerfacecolor', colors.red2);
% legend([p_landing], 'Tangent Impact', 'fontsize', 20, 'fontname', 'Times New Roman', 'location', 'northeast')
% title('DRO ... s1','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)


figure; hold all
xlim([0 360])
ylim([-90 90])
PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')
h = image(xlim, -ylim, bodies.europa.img);
[lons_new_ZH] = convert_lon180_to_lon360(latLons_DRO__s1_ZH(:,2));
plot(lons_new_ZH,latLons_DRO__s1_ZH(:,1),'.', 'color', colors.red2)
[lons_new_impact_ZH] = convert_lon180_to_lon360(latLons_apses_DRO__s1_ZH((apses_alts_ZH < 1/rNorm),2));
p_landing_ZH = plot(lons_new_impact_ZH, latLons_apses_DRO__s1_ZH((apses_alts_ZH < 1/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.black, 'markerfacecolor', colors.drkgrey);
[lons_new_impact_close] = convert_lon180_to_lon360(latLons_apses_DRO__s1_ZH(close_indices_ZH,2));
p_close = plot(lons_new_impact_close, latLons_apses_DRO__s1_ZH(close_indices_ZH,1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.drkblue, 'markerfacecolor', colors.blue2);
legend([p_landing_ZH, p_close], 'Tangent Impact', '18 km Apoapsis', 'fontsize', 20, 'fontname', 'Times New Roman', 'location', 'northeast')
title('DRO ... s1 ... ZH','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)



% figure; hold all
% xlim([0 360])
% ylim([-90 90])
% PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')
% h = image(xlim, -ylim, bodies.europa.img);
% plot(lons_new,latLons_DRO__s1(:,1),'.','color',colors.blue2)
% plot(lons_new_ZH,latLons_DRO__s1_ZH(:,1),'.', 'color', colors.red2)
% p_landing = plot(lons_new_impact, latLons_apses_DRO__s1((apses_alts < 1/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.drkblue, 'markerfacecolor', colors.blue2);
% p_landing_ZH = plot(lons_new_impact_ZH, latLons_apses_DRO__s1_ZH((apses_alts_ZH < 1/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.drkred, 'markerfacecolor', colors.red2);
% legend([p_landing, p_landing_ZH], 'Tangent Impact - CR3BP', 'Tangent Impact - CR3BP w ZH', 'fontsize', 20, 'fontname', 'Times New Roman', 'location', 'northeast')
% title('DRO ... s1 ... Both','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)

end % run_DRO__s1
%% ------------------------------------------------
%%% DPO__s1
% -------------------------------------------------
if run_DPO__s1
fprintf('---------------------------------')
fprintf('DPO__s1 (planar)\n')

[T_DPO__s1, X_DPO__s1, t_ev, x_ev, in_ev] = ode113(@Int_CR3BnSTM, [0, PO_DPO__s1(end)], [PO_DPO__s1(1:6); stm0_colVec], options_apsis, prms);
apses_alts = rowNorm(x_ev(:,1:3) - [1-prms.u,0,0]) - prms.R2;
% figure; hold all
% plot3(X_DPO__s1(:,1),X_DPO__s1(:,2),X_DPO__s1(:,3), 'linewidth', 2, 'color', colors.black)
% plotSecondary(secondary)
% PlotBoi3_CR3Bn(26)
% axis equal
% title('DPO ... s1','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)
% 
% plot3(x_ev((apses_alts < 5/rNorm),1), x_ev((apses_alts < 5/rNorm),2), x_ev((apses_alts < 5/rNorm),3), '.','markersize',20, 'color', colors.red2)
% set(gca,'yticklabel',num2str(get(gca,'ytick')','%1.3f'))
% set(gca,'zticklabel',num2str(get(gca,'ztick')','%1.3f'))



[T_DPO__s1_ZH, X_DPO__s1_ZH, t_ev_ZH, x_ev_ZH, in_ev_ZH] = ode113(@Int_CR3BnSTM_J2pJ4pJ6pJ2s, [0, PO_DPO__s1_ZH(end)], [PO_DPO__s1_ZH(1:6); stm0_colVec], options_apsis, prms_ZH);
fprintf('Prop Error = %1.2e\n', norm(X_DPO__s1_ZH(end,1:6) - X_DPO__s1_ZH(1,1:6)) / norm(X_DPO__s1_ZH(1,1:6)))
figure; hold all
% plot3(X_DPO__s1_ZH(:,1),X_DPO__s1_ZH(:,2),X_DPO__s1_ZH(:,3), 'linewidth', 2, 'color', colors.red2)
plot3(X_DPO__s1_ZH(:,1),X_DPO__s1_ZH(:,2),X_DPO__s1_ZH(:,3), 'linewidth', 2, 'color', colors.mag)
plotSecondary(secondary)
PlotBoi3_CR3Bn(26)
axis equal
title('DPO ... s1 ... ZH','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)

apses_alts_ZH = rowNorm(x_ev_ZH(:,1:3) - [1-prms_ZH.u,0,0]) - prms_ZH.R2;
plot3(x_ev_ZH((apses_alts_ZH < 5/rNorm),1), x_ev_ZH((apses_alts_ZH < 5/rNorm),2), x_ev_ZH((apses_alts_ZH < 5/rNorm),3), '.','markersize',20, 'color', colors.black)
set(gca,'yticklabel',num2str(get(gca,'ytick')','%1.3f'))
set(gca,'zticklabel',num2str(get(gca,'ztick')','%1.3f'))

view(0,90)

% figure; hold all
% p_3B = plot3(X_DPO__s1(:,1),X_DPO__s1(:,2),X_DPO__s1(:,3), 'linewidth', 2, 'color', colors.blue2);
% p_ZH = plot3(X_DPO__s1_ZH(:,1),X_DPO__s1_ZH(:,2),X_DPO__s1_ZH(:,3), 'linewidth', 2, 'color', colors.red2);
% plotSecondary(secondary)
% PlotBoi3_CR3Bn(26)
% axis equal
% legend([p_3B, p_ZH], 'CR3BP', 'CR3BP w ZH')
% set(gca,'yticklabel',num2str(get(gca,'ytick')','%1.3f'))
% set(gca,'zticklabel',num2str(get(gca,'ztick')','%1.3f'))


stm_tf_t0                           = reshape(X_DPO__s1(end,7:42),6,6);
monodromy                           = stm_tf_t0;
[eigenVectors_new, eigenValues_new] = eig(monodromy);
[S1, S2]                            = getStabilityIndices(diag(eigenValues_new));

stm_tf_t0_ZH                              = reshape(X_DPO__s1_ZH(end,7:42),6,6);
monodromy_ZH                              = stm_tf_t0_ZH;
[eigenVectors_new_ZH, eigenValues_new_ZH] = eig(monodromy_ZH);
[S1_ZH, S2_ZH]                            = getStabilityIndices(diag(eigenValues_new_ZH));

fprintf('[S1, S2]       = [%1.3f, %1.3f]\n', S1, S2)
fprintf('[S1_ZH, S2_ZH] = [%1.3f, %1.3f]\n', S1_ZH, S2_ZH)


latLons_DPO__s1 = zeros(length(T_DPO__s1),2);
for kk = 1:length(T_DPO__s1)
    [lat_deg, lon_deg] = BCR2latlon(X_DPO__s1(kk,1:3)', 'secondary', prms.u);
    latLons_DPO__s1(kk,:) = [lat_deg, lon_deg];
end
latLons_apses_DPO__s1 = zeros(length(t_ev),2);
for kk = 1:length(t_ev)
    [lat_deg, lon_deg] = BCR2latlon(x_ev(kk,1:3)', 'secondary', prms.u);
    latLons_apses_DPO__s1(kk,:) = [lat_deg, lon_deg];
end


latLons_DPO__s1_ZH = zeros(length(T_DPO__s1_ZH),2);
for kk = 1:length(T_DPO__s1_ZH)
    [lat_deg, lon_deg] = BCR2latlon(X_DPO__s1_ZH(kk,1:3)', 'secondary', prms_ZH.u);
    latLons_DPO__s1_ZH(kk,:) = [lat_deg, lon_deg];
end
latLons_apses_DPO__s1_ZH = zeros(length(t_ev_ZH),2);
for kk = 1:length(t_ev_ZH)
    [lat_deg, lon_deg] = BCR2latlon(x_ev_ZH(kk,1:3)', 'secondary', prms_ZH.u);
    latLons_apses_DPO__s1_ZH(kk,:) = [lat_deg, lon_deg];
end

% figure; hold all
% xlim([0 360])
% ylim([-90 90])
% PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')
% h = image(xlim, -ylim, bodies.europa.img);
% [lons_new] = convert_lon180_to_lon360(latLons_DPO__s1(:,2));
% plot(lons_new,latLons_DPO__s1(:,1),'k.')
% [lons_new_impact] = convert_lon180_to_lon360(latLons_apses_DPO__s1((apses_alts < 5/rNorm),2));
% p_landing = plot(lons_new_impact, latLons_apses_DPO__s1((apses_alts < 5/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.drkred, 'markerfacecolor', colors.red2);
% legend([p_landing], 'Tangent Impact', 'fontsize', 20, 'fontname', 'Times New Roman', 'location', 'northeast')
% title('DPO ... s1','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)


figure; hold all
xlim([0 360])
ylim([-90 90])
PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')
h = image(xlim, -ylim, bodies.europa.img);
[lons_new_ZH] = convert_lon180_to_lon360(latLons_DPO__s1_ZH(:,2));
plot(lons_new_ZH,latLons_DPO__s1_ZH(:,1),'.', 'color', colors.red2)
[lons_new_impact_ZH] = convert_lon180_to_lon360(latLons_apses_DPO__s1_ZH((apses_alts_ZH < 5/rNorm),2));
p_landing_ZH = plot(lons_new_impact_ZH, latLons_apses_DPO__s1_ZH((apses_alts_ZH < 5/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.black, 'markerfacecolor', colors.drkgrey);
legend([p_landing_ZH], 'Tangent Impact', 'fontsize', 20, 'fontname', 'Times New Roman', 'location', 'northeast')
title('DPO ... s1 ... ZH','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)

% figure; hold all
% xlim([0 360])
% ylim([-90 90])
% PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')
% h = image(xlim, -ylim, bodies.europa.img);
% plot(lons_new,latLons_DPO__s1(:,1),'.','color',colors.blue2)
% plot(lons_new_ZH,latLons_DPO__s1_ZH(:,1),'.', 'color', colors.red2)
% p_landing = plot(lons_new_impact, latLons_apses_DPO__s1((apses_alts < 5/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.drkblue, 'markerfacecolor', colors.blue2);
% p_landing_ZH = plot(lons_new_impact_ZH, latLons_apses_DPO__s1_ZH((apses_alts_ZH < 5/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.drkred, 'markerfacecolor', colors.red2);
% legend([p_landing, p_landing_ZH], 'Tangent Impact - CR3BP', 'Tangent Impact - CR3BP w ZH', 'fontsize', 20, 'fontname', 'Times New Roman', 'location', 'northeast')
% title('DPO ... s1 ... Both','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)

end % run_DPO__s1
%% ------------------------------------------------
%%% DPO_2P3__ns1
% -------------------------------------------------
if run_DPO_2P3__ns1
fprintf('---------------------------------')
fprintf('DPO_2P3__ns1 (planar)\n')

[T_DPO_2P3__ns1, X_DPO_2P3__ns1, t_ev, x_ev, in_ev] = ode113(@Int_CR3BnSTM, [0, PO_DPO_2P3__ns1(end)], [PO_DPO_2P3__ns1(1:6); stm0_colVec], options_apsis, prms);
apses_alts = rowNorm(x_ev(:,1:3) - [1-prms.u,0,0]) - prms.R2;
% figure; hold all
% plot3(X_DPO_2P3__ns1(:,1),X_DPO_2P3__ns1(:,2),X_DPO_2P3__ns1(:,3), 'linewidth', 2, 'color', colors.black)
% plotSecondary(secondary)
% PlotBoi3_CR3Bn(26)
% axis equal
% title('DPO\_2P3 ... ns1','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)
% 
% plot3(x_ev((apses_alts < 5/rNorm),1), x_ev((apses_alts < 5/rNorm),2), x_ev((apses_alts < 5/rNorm),3), '.','markersize',20, 'color', colors.red2)
% set(gca,'yticklabel',num2str(get(gca,'ytick')','%1.3f'))
% set(gca,'zticklabel',num2str(get(gca,'ztick')','%1.3f'))



[T_DPO_2P3__ns1_ZH, X_DPO_2P3__ns1_ZH, t_ev_ZH, x_ev_ZH, in_ev_ZH] = ode113(@Int_CR3BnSTM_J2pJ4pJ6pJ2s, [0, PO_DPO_2P3__ns1_ZH(end)], [PO_DPO_2P3__ns1_ZH(1:6); stm0_colVec], options_apsis, prms_ZH);
fprintf('Prop Error = %1.2e\n', norm(X_DPO_2P3__ns1_ZH(end,1:6) - X_DPO_2P3__ns1_ZH(1,1:6)) / norm(X_DPO_2P3__ns1_ZH(1,1:6)))
figure; hold all
% plot3(X_DPO_2P3__ns1_ZH(:,1),X_DPO_2P3__ns1_ZH(:,2),X_DPO_2P3__ns1_ZH(:,3), 'linewidth', 2, 'color', colors.red2)
plot3(X_DPO_2P3__ns1_ZH(:,1),X_DPO_2P3__ns1_ZH(:,2),X_DPO_2P3__ns1_ZH(:,3), 'linewidth', 2, 'color', colors.mag)
plotSecondary(secondary)
PlotBoi3_CR3Bn(26)
axis equal
title('DPO\_2P3 ... ns1 ... ZH','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)

apses_alts_ZH = rowNorm(x_ev_ZH(:,1:3) - [1-prms_ZH.u,0,0]) - prms_ZH.R2;
close_indices_ZH = ((183/rNorm < apses_alts_ZH) & (apses_alts_ZH < 184/rNorm));
plot3(x_ev_ZH(close_indices_ZH, 1), x_ev_ZH(close_indices_ZH, 2), x_ev_ZH(close_indices_ZH, 3), '.','markersize',20, 'color', colors.blue2)

plot3(x_ev_ZH((apses_alts_ZH < 5/rNorm),1), x_ev_ZH((apses_alts_ZH < 5/rNorm),2), x_ev_ZH((apses_alts_ZH < 5/rNorm),3), '.','markersize',20, 'color', colors.black)
set(gca,'yticklabel',num2str(get(gca,'ytick')','%1.3f'))
set(gca,'zticklabel',num2str(get(gca,'ztick')','%1.3f'))


view(0,90)

% figure; hold all
% p_3B = plot3(X_DPO_2P3__ns1(:,1),X_DPO_2P3__ns1(:,2),X_DPO_2P3__ns1(:,3), 'linewidth', 2, 'color', colors.blue2);
% p_ZH = plot3(X_DPO_2P3__ns1_ZH(:,1),X_DPO_2P3__ns1_ZH(:,2),X_DPO_2P3__ns1_ZH(:,3), 'linewidth', 2, 'color', colors.red2);
% plotSecondary(secondary)
% PlotBoi3_CR3Bn(26)
% axis equal
% legend([p_3B, p_ZH], 'CR3BP', 'CR3BP w ZH')
% set(gca,'yticklabel',num2str(get(gca,'ytick')','%1.3f'))
% set(gca,'zticklabel',num2str(get(gca,'ztick')','%1.3f'))


stm_tf_t0                           = reshape(X_DPO_2P3__ns1(end,7:42),6,6);
monodromy                           = stm_tf_t0;
[eigenVectors_new, eigenValues_new] = eig(monodromy);
[S1, S2]                            = getStabilityIndices(diag(eigenValues_new));

stm_tf_t0_ZH                              = reshape(X_DPO_2P3__ns1_ZH(end,7:42),6,6);
monodromy_ZH                              = stm_tf_t0_ZH;
[eigenVectors_new_ZH, eigenValues_new_ZH] = eig(monodromy_ZH);
[S1_ZH, S2_ZH]                            = getStabilityIndices(diag(eigenValues_new_ZH));

fprintf('[S1, S2]       = [%1.3f, %1.3f]\n', S1, S2)
fprintf('[S1_ZH, S2_ZH] = [%1.3f, %1.3f]\n', S1_ZH, S2_ZH)


latLons_DPO_2P3__ns1 = zeros(length(T_DPO_2P3__ns1),2);
for kk = 1:length(T_DPO_2P3__ns1)
    [lat_deg, lon_deg] = BCR2latlon(X_DPO_2P3__ns1(kk,1:3)', 'secondary', prms.u);
    latLons_DPO_2P3__ns1(kk,:) = [lat_deg, lon_deg];
end
latLons_apses_DPO_2P3__ns1 = zeros(length(t_ev),2);
for kk = 1:length(t_ev)
    [lat_deg, lon_deg] = BCR2latlon(x_ev(kk,1:3)', 'secondary', prms.u);
    latLons_apses_DPO_2P3__ns1(kk,:) = [lat_deg, lon_deg];
end


latLons_DPO_2P3__ns1_ZH = zeros(length(T_DPO_2P3__ns1_ZH),2);
for kk = 1:length(T_DPO_2P3__ns1_ZH)
    [lat_deg, lon_deg] = BCR2latlon(X_DPO_2P3__ns1_ZH(kk,1:3)', 'secondary', prms_ZH.u);
    latLons_DPO_2P3__ns1_ZH(kk,:) = [lat_deg, lon_deg];
end
latLons_apses_DPO_2P3__ns1_ZH = zeros(length(t_ev_ZH),2);
for kk = 1:length(t_ev_ZH)
    [lat_deg, lon_deg] = BCR2latlon(x_ev_ZH(kk,1:3)', 'secondary', prms_ZH.u);
    latLons_apses_DPO_2P3__ns1_ZH(kk,:) = [lat_deg, lon_deg];
end

% figure; hold all
% xlim([0 360])
% ylim([-90 90])
% PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')
% h = image(xlim, -ylim, bodies.europa.img);
% [lons_new] = convert_lon180_to_lon360(latLons_DPO_2P3__ns1(:,2));
% plot(lons_new,latLons_DPO_2P3__ns1(:,1),'k.')
% [lons_new_impact] = convert_lon180_to_lon360(latLons_apses_DPO_2P3__ns1((apses_alts < 5/rNorm),2));
% p_landing = plot(lons_new_impact, latLons_apses_DPO_2P3__ns1((apses_alts < 5/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.drkred, 'markerfacecolor', colors.red2);
% legend([p_landing], 'Tangent Impact', 'fontsize', 20, 'fontname', 'Times New Roman', 'location', 'northeast')
% title('DPO\_2P3 ... ns1','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)


figure; hold all
xlim([0 360])
ylim([-90 90])
PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')
h = image(xlim, -ylim, bodies.europa.img);
[lons_new_ZH] = convert_lon180_to_lon360(latLons_DPO_2P3__ns1_ZH(:,2));
plot(lons_new_ZH,latLons_DPO_2P3__ns1_ZH(:,1),'.', 'color', colors.red2)
[lons_new_impact_close] = convert_lon180_to_lon360(latLons_apses_DPO_2P3__ns1_ZH(close_indices_ZH,2));
p_close = plot(lons_new_impact_close, latLons_apses_DPO_2P3__ns1_ZH(close_indices_ZH,1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.drkblue, 'markerfacecolor', colors.blue2);
[lons_new_impact_ZH] = convert_lon180_to_lon360(latLons_apses_DPO_2P3__ns1_ZH((apses_alts_ZH < 5/rNorm),2));
p_landing_ZH = plot(lons_new_impact_ZH, latLons_apses_DPO_2P3__ns1_ZH((apses_alts_ZH < 5/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.black, 'markerfacecolor', colors.drkgrey);
legend([p_landing_ZH, p_close], 'Tangent Impact', '183 km Periapsis', 'fontsize', 20, 'fontname', 'Times New Roman', 'location', 'northeast')
title('DPO\_2P3 ... ns1 ... ZH','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)

% figure; hold all
% xlim([0 360])
% ylim([-90 90])
% PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')
% h = image(xlim, -ylim, bodies.europa.img);
% plot(lons_new,latLons_DPO_2P3__ns1(:,1),'.','color',colors.blue2)
% plot(lons_new_ZH,latLons_DPO_2P3__ns1_ZH(:,1),'.', 'color', colors.red2)
% p_landing = plot(lons_new_impact, latLons_apses_DPO_2P3__ns1((apses_alts < 5/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.drkblue, 'markerfacecolor', colors.blue2);
% p_landing_ZH = plot(lons_new_impact_ZH, latLons_apses_DPO_2P3__ns1_ZH((apses_alts_ZH < 5/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.drkred, 'markerfacecolor', colors.red2);
% legend([p_landing, p_landing_ZH], 'Tangent Impact - CR3BP', 'Tangent Impact - CR3BP w ZH', 'fontsize', 20, 'fontname', 'Times New Roman', 'location', 'northeast')
% title('DPO\_2P3 ... ns1 ... Both','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)

end % run_DPO_2P3__ns1

%% ------------------------------------------------
%%% DPO_2P4__ns1
% -------------------------------------------------
if run_DPO_2P4__ns1
fprintf('---------------------------------')
fprintf('DPO_2P4__ns1 (planar)\n')

[T_DPO_2P4__ns1, X_DPO_2P4__ns1, t_ev, x_ev, in_ev] = ode113(@Int_CR3BnSTM, [0, PO_DPO_2P4__ns1(end)], [PO_DPO_2P4__ns1(1:6); stm0_colVec], options_apsis, prms);
apses_alts = rowNorm(x_ev(:,1:3) - [1-prms.u,0,0]) - prms.R2;
% figure; hold all
% plot3(X_DPO_2P4__ns1(:,1),X_DPO_2P4__ns1(:,2),X_DPO_2P4__ns1(:,3), 'linewidth', 2, 'color', colors.black)
% plotSecondary(secondary)
% PlotBoi3_CR3Bn(26)
% axis equal
% title('DPO\_2P4 ... ns1','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)
% 
% plot3(x_ev((apses_alts < 5/rNorm),1), x_ev((apses_alts < 5/rNorm),2), x_ev((apses_alts < 5/rNorm),3), '.','markersize',20, 'color', colors.red2)
% set(gca,'yticklabel',num2str(get(gca,'ytick')','%1.3f'))
% set(gca,'zticklabel',num2str(get(gca,'ztick')','%1.3f'))

[T_DPO_2P4__ns1_ZH, X_DPO_2P4__ns1_ZH, t_ev_ZH, x_ev_ZH, in_ev_ZH] = ode113(@Int_CR3BnSTM_J2pJ4pJ6pJ2s, [0, PO_DPO_2P4__ns1_ZH(end)], [PO_DPO_2P4__ns1_ZH(1:6); stm0_colVec], options_apsis, prms_ZH);
fprintf('Prop Error = %1.2e\n', norm(X_DPO_2P4__ns1_ZH(end,1:6) - X_DPO_2P4__ns1_ZH(1,1:6)) / norm(X_DPO_2P4__ns1_ZH(1,1:6)))
figure; hold all
% plot3(X_DPO_2P4__ns1_ZH(:,1),X_DPO_2P4__ns1_ZH(:,2),X_DPO_2P4__ns1_ZH(:,3), 'linewidth', 2, 'color', colors.red2)
plot3(X_DPO_2P4__ns1_ZH(:,1),X_DPO_2P4__ns1_ZH(:,2),X_DPO_2P4__ns1_ZH(:,3), 'linewidth', 2, 'color', colors.mag)
plotSecondary(secondary)
PlotBoi3_CR3Bn(26)
axis equal
title('DPO\_2P4 ... ns1 ... ZH','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)

apses_alts_ZH = rowNorm(x_ev_ZH(:,1:3) - [1-prms_ZH.u,0,0]) - prms_ZH.R2;
plot3(x_ev_ZH((apses_alts_ZH < 5/rNorm),1), x_ev_ZH((apses_alts_ZH < 5/rNorm),2), x_ev_ZH((apses_alts_ZH < 5/rNorm),3), '.','markersize',20, 'color', colors.black)
set(gca,'yticklabel',num2str(get(gca,'ytick')','%1.3f'))
set(gca,'zticklabel',num2str(get(gca,'ztick')','%1.3f'))

view(0,90)

% 
% figure; hold all
% p_3B = plot3(X_DPO_2P4__ns1(:,1),X_DPO_2P4__ns1(:,2),X_DPO_2P4__ns1(:,3), 'linewidth', 2, 'color', colors.blue2);
% p_ZH = plot3(X_DPO_2P4__ns1_ZH(:,1),X_DPO_2P4__ns1_ZH(:,2),X_DPO_2P4__ns1_ZH(:,3), 'linewidth', 2, 'color', colors.red2);
% plotSecondary(secondary)
% PlotBoi3_CR3Bn(26)
% axis equal
% legend([p_3B, p_ZH], 'CR3BP', 'CR3BP w ZH')
% set(gca,'yticklabel',num2str(get(gca,'ytick')','%1.3f'))
% set(gca,'zticklabel',num2str(get(gca,'ztick')','%1.3f'))


stm_tf_t0                           = reshape(X_DPO_2P4__ns1(end,7:42),6,6);
monodromy                           = stm_tf_t0;
[eigenVectors_new, eigenValues_new] = eig(monodromy);
[S1, S2]                            = getStabilityIndices(diag(eigenValues_new));

stm_tf_t0_ZH                              = reshape(X_DPO_2P4__ns1_ZH(end,7:42),6,6);
monodromy_ZH                              = stm_tf_t0_ZH;
[eigenVectors_new_ZH, eigenValues_new_ZH] = eig(monodromy_ZH);
[S1_ZH, S2_ZH]                            = getStabilityIndices(diag(eigenValues_new_ZH));

fprintf('[S1, S2]       = [%1.3f, %1.3f]\n', S1, S2)
fprintf('[S1_ZH, S2_ZH] = [%1.3f, %1.3f]\n', S1_ZH, S2_ZH)


latLons_DPO_2P4__ns1 = zeros(length(T_DPO_2P4__ns1),2);
for kk = 1:length(T_DPO_2P4__ns1)
    [lat_deg, lon_deg] = BCR2latlon(X_DPO_2P4__ns1(kk,1:3)', 'secondary', prms.u);
    latLons_DPO_2P4__ns1(kk,:) = [lat_deg, lon_deg];
end
latLons_apses_DPO_2P4__ns1 = zeros(length(t_ev),2);
for kk = 1:length(t_ev)
    [lat_deg, lon_deg] = BCR2latlon(x_ev(kk,1:3)', 'secondary', prms.u);
    latLons_apses_DPO_2P4__ns1(kk,:) = [lat_deg, lon_deg];
end


latLons_DPO_2P4__ns1_ZH = zeros(length(T_DPO_2P4__ns1_ZH),2);
for kk = 1:length(T_DPO_2P4__ns1_ZH)
    [lat_deg, lon_deg] = BCR2latlon(X_DPO_2P4__ns1_ZH(kk,1:3)', 'secondary', prms_ZH.u);
    latLons_DPO_2P4__ns1_ZH(kk,:) = [lat_deg, lon_deg];
end
latLons_apses_DPO_2P4__ns1_ZH = zeros(length(t_ev_ZH),2);
for kk = 1:length(t_ev_ZH)
    [lat_deg, lon_deg] = BCR2latlon(x_ev_ZH(kk,1:3)', 'secondary', prms_ZH.u);
    latLons_apses_DPO_2P4__ns1_ZH(kk,:) = [lat_deg, lon_deg];
end

% figure; hold all
% xlim([0 360])
% ylim([-90 90])
% PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')
% h = image(xlim, -ylim, bodies.europa.img);
% [lons_new] = convert_lon180_to_lon360(latLons_DPO_2P4__ns1(:,2));
% plot(lons_new,latLons_DPO_2P4__ns1(:,1),'k.')
% [lons_new_impact] = convert_lon180_to_lon360(latLons_apses_DPO_2P4__ns1((apses_alts < 5/rNorm),2));
% p_landing = plot(lons_new_impact, latLons_apses_DPO_2P4__ns1((apses_alts < 5/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.drkred, 'markerfacecolor', colors.red2);
% legend([p_landing], 'Tangent Impact', 'fontsize', 20, 'fontname', 'Times New Roman', 'location', 'northeast')
% title('DPO\_2P4 ... ns1','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)


figure; hold all
xlim([0 360])
ylim([-90 90])
PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')
h = image(xlim, -ylim, bodies.europa.img);
[lons_new_ZH] = convert_lon180_to_lon360(latLons_DPO_2P4__ns1_ZH(:,2));
plot(lons_new_ZH,latLons_DPO_2P4__ns1_ZH(:,1),'.', 'color', colors.red2)
[lons_new_impact_ZH] = convert_lon180_to_lon360(latLons_apses_DPO_2P4__ns1_ZH((apses_alts_ZH < 5/rNorm),2));
p_landing_ZH = plot(lons_new_impact_ZH, latLons_apses_DPO_2P4__ns1_ZH((apses_alts_ZH < 5/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.black, 'markerfacecolor', colors.drkgrey);
legend([p_landing_ZH], 'Tangent Impact', 'fontsize', 20, 'fontname', 'Times New Roman', 'location', 'northeast')
title('DPO\_2P4 ... ns1 ... ZH','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)

% figure; hold all
% xlim([0 360])
% ylim([-90 90])
% PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')
% h = image(xlim, -ylim, bodies.europa.img);
% plot(lons_new,latLons_DPO_2P4__ns1(:,1),'.','color',colors.blue2)
% plot(lons_new_ZH,latLons_DPO_2P4__ns1_ZH(:,1),'.', 'color', colors.red2)
% p_landing = plot(lons_new_impact, latLons_apses_DPO_2P4__ns1((apses_alts < 5/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.drkblue, 'markerfacecolor', colors.blue2);
% p_landing_ZH = plot(lons_new_impact_ZH, latLons_apses_DPO_2P4__ns1_ZH((apses_alts_ZH < 5/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.drkred, 'markerfacecolor', colors.red2);
% legend([p_landing, p_landing_ZH], 'Tangent Impact - CR3BP', 'Tangent Impact - CR3BP w ZH', 'fontsize', 20, 'fontname', 'Times New Roman', 'location', 'northeast')
% title('DPO\_2P4 ... ns1 ... Both','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)

end % run_DPO_2P4__ns1
%% ------------------------------------------------
%%% Hg2__ns1
% -------------------------------------------------
if run_Hg2__ns1
fprintf('---------------------------------')
fprintf('Hg2__ns1 (planar)\n')

[T_Hg2__ns1, X_Hg2__ns1, t_ev, x_ev, in_ev] = ode113(@Int_CR3BnSTM, [0, PO_Hg2__ns1(end)], [PO_Hg2__ns1(1:6); stm0_colVec], options_apsis, prms);
apses_alts = rowNorm(x_ev(:,1:3) - [1-prms.u,0,0]) - prms.R2;
% figure; hold all
% plot3(X_Hg2__ns1(:,1),X_Hg2__ns1(:,2),X_Hg2__ns1(:,3), 'linewidth', 2, 'color', colors.black)
% plotSecondary(secondary)
% PlotBoi3_CR3Bn(26)
% axis equal
% title('Hg2 ... ns1','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)
% 
% plot3(x_ev((apses_alts < 5/rNorm),1), x_ev((apses_alts < 5/rNorm),2), x_ev((apses_alts < 5/rNorm),3), '.','markersize',20, 'color', colors.red2)
% set(gca,'yticklabel',num2str(get(gca,'ytick')','%1.3f'))
% set(gca,'zticklabel',num2str(get(gca,'ztick')','%1.3f'))

[T_Hg2__ns1_ZH, X_Hg2__ns1_ZH, t_ev_ZH, x_ev_ZH, in_ev_ZH] = ode113(@Int_CR3BnSTM_J2pJ4pJ6pJ2s, [0, PO_Hg2__ns1_ZH(end)], [PO_Hg2__ns1_ZH(1:6); stm0_colVec], options_apsis, prms_ZH);
fprintf('Prop Error = %1.2e\n', norm(X_Hg2__ns1_ZH(end,1:6) - X_Hg2__ns1_ZH(1,1:6)) / norm(X_Hg2__ns1_ZH(1,1:6)))
figure; hold all
% plot3(X_Hg2__ns1_ZH(:,1),X_Hg2__ns1_ZH(:,2),X_Hg2__ns1_ZH(:,3), 'linewidth', 2, 'color', colors.red2)
plot3(X_Hg2__ns1_ZH(:,1),X_Hg2__ns1_ZH(:,2),X_Hg2__ns1_ZH(:,3), 'linewidth', 2, 'color', colors.mag)
plotSecondary(secondary)
PlotBoi3_CR3Bn(26)
axis equal
title('Hg2 ... ns1 ... ZH','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)

apses_alts_ZH = rowNorm(x_ev_ZH(:,1:3) - [1-prms_ZH.u,0,0]) - prms_ZH.R2;
plot3(x_ev_ZH((apses_alts_ZH < 5/rNorm),1), x_ev_ZH((apses_alts_ZH < 5/rNorm),2), x_ev_ZH((apses_alts_ZH < 5/rNorm),3), '.','markersize',20, 'color', colors.black)
set(gca,'yticklabel',num2str(get(gca,'ytick')','%1.3f'))
set(gca,'zticklabel',num2str(get(gca,'ztick')','%1.3f'))

close_indices_ZH = ((59/rNorm < apses_alts_ZH) & (apses_alts_ZH < 60/rNorm));
plot3(x_ev_ZH(close_indices_ZH, 1), x_ev_ZH(close_indices_ZH, 2), x_ev_ZH(close_indices_ZH, 3), '.','markersize',20, 'color', colors.blue2)


view(0,90)

% figure; hold all
% p_3B = plot3(X_Hg2__ns1(:,1),X_Hg2__ns1(:,2),X_Hg2__ns1(:,3), 'linewidth', 2, 'color', colors.blue2);
% p_ZH = plot3(X_Hg2__ns1_ZH(:,1),X_Hg2__ns1_ZH(:,2),X_Hg2__ns1_ZH(:,3), 'linewidth', 2, 'color', colors.red2);
% plotSecondary(secondary)
% PlotBoi3_CR3Bn(26)
% axis equal
% legend([p_3B, p_ZH], 'CR3BP', 'CR3BP w ZH')
% set(gca,'yticklabel',num2str(get(gca,'ytick')','%1.3f'))
% set(gca,'zticklabel',num2str(get(gca,'ztick')','%1.3f'))


stm_tf_t0                           = reshape(X_Hg2__ns1(end,7:42),6,6);
monodromy                           = stm_tf_t0;
[eigenVectors_new, eigenValues_new] = eig(monodromy);
[S1, S2]                            = getStabilityIndices(diag(eigenValues_new));

stm_tf_t0_ZH                              = reshape(X_Hg2__ns1_ZH(end,7:42),6,6);
monodromy_ZH                              = stm_tf_t0_ZH;
[eigenVectors_new_ZH, eigenValues_new_ZH] = eig(monodromy_ZH);
[S1_ZH, S2_ZH]                            = getStabilityIndices(diag(eigenValues_new_ZH));

fprintf('[S1, S2]       = [%1.3f, %1.3f]\n', S1, S2)
fprintf('[S1_ZH, S2_ZH] = [%1.3f, %1.3f]\n', S1_ZH, S2_ZH)


latLons_Hg2__ns1 = zeros(length(T_Hg2__ns1),2);
for kk = 1:length(T_Hg2__ns1)
    [lat_deg, lon_deg] = BCR2latlon(X_Hg2__ns1(kk,1:3)', 'secondary', prms.u);
    latLons_Hg2__ns1(kk,:) = [lat_deg, lon_deg];
end
latLons_apses_Hg2__ns1 = zeros(length(t_ev),2);
for kk = 1:length(t_ev)
    [lat_deg, lon_deg] = BCR2latlon(x_ev(kk,1:3)', 'secondary', prms.u);
    latLons_apses_Hg2__ns1(kk,:) = [lat_deg, lon_deg];
end


latLons_Hg2__ns1_ZH = zeros(length(T_Hg2__ns1_ZH),2);
for kk = 1:length(T_Hg2__ns1_ZH)
    [lat_deg, lon_deg] = BCR2latlon(X_Hg2__ns1_ZH(kk,1:3)', 'secondary', prms_ZH.u);
    latLons_Hg2__ns1_ZH(kk,:) = [lat_deg, lon_deg];
end
latLons_apses_Hg2__ns1_ZH = zeros(length(t_ev_ZH),2);
for kk = 1:length(t_ev_ZH)
    [lat_deg, lon_deg] = BCR2latlon(x_ev_ZH(kk,1:3)', 'secondary', prms_ZH.u);
    latLons_apses_Hg2__ns1_ZH(kk,:) = [lat_deg, lon_deg];
end

% figure; hold all
% xlim([0 360])
% ylim([-90 90])
% PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')
% h = image(xlim, -ylim, bodies.europa.img);
% [lons_new] = convert_lon180_to_lon360(latLons_Hg2__ns1(:,2));
% plot(lons_new,latLons_Hg2__ns1(:,1),'k.')
% [lons_new_impact] = convert_lon180_to_lon360(latLons_apses_Hg2__ns1((apses_alts < 5/rNorm),2));
% p_landing = plot(lons_new_impact, latLons_apses_Hg2__ns1((apses_alts < 5/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.drkred, 'markerfacecolor', colors.red2);
% legend([p_landing], 'Tangent Impact', 'fontsize', 20, 'fontname', 'Times New Roman', 'location', 'northeast')
% title('Hg2 ... ns1','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)


figure; hold all
xlim([0 360])
ylim([-90 90])
PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')
h = image(xlim, -ylim, bodies.europa.img);
[lons_new_ZH] = convert_lon180_to_lon360(latLons_Hg2__ns1_ZH(:,2));
plot(lons_new_ZH,latLons_Hg2__ns1_ZH(:,1),'.', 'color', colors.red2)
[lons_new_impact_ZH] = convert_lon180_to_lon360(latLons_apses_Hg2__ns1_ZH((apses_alts_ZH < 5/rNorm),2));
p_landing_ZH = plot(lons_new_impact_ZH, latLons_apses_Hg2__ns1_ZH((apses_alts_ZH < 5/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.black, 'markerfacecolor', colors.drkgrey);
[lons_new_impact_close] = convert_lon180_to_lon360(latLons_apses_Hg2__ns1_ZH(close_indices_ZH,2));
p_close = plot(lons_new_impact_close, latLons_apses_Hg2__ns1_ZH(close_indices_ZH,1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.drkblue, 'markerfacecolor', colors.blue2);
legend([p_landing_ZH, p_close], 'Tangent Impact', '59 km Periapsis', 'fontsize', 20, 'fontname', 'Times New Roman', 'location', 'northeast')
title('Hg2 ... ns1 ... ZH','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)


%legend([p_landing_ZH, p_close], 'Tangent Impact', '18 km Periapsis', 'fontsize', 20, 'fontname', 'Times New Roman', 'location', 'northeast')

% figure; hold all
% xlim([0 360])
% ylim([-90 90])
% PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')
% h = image(xlim, -ylim, bodies.europa.img);
% plot(lons_new,latLons_Hg2__ns1(:,1),'.','color',colors.blue2)
% plot(lons_new_ZH,latLons_Hg2__ns1_ZH(:,1),'.', 'color', colors.red2)
% p_landing = plot(lons_new_impact, latLons_apses_Hg2__ns1((apses_alts < 5/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.drkblue, 'markerfacecolor', colors.blue2);
% p_landing_ZH = plot(lons_new_impact_ZH, latLons_apses_Hg2__ns1_ZH((apses_alts_ZH < 5/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.drkred, 'markerfacecolor', colors.red2);
% legend([p_landing, p_landing_ZH], 'Tangent Impact - CR3BP', 'Tangent Impact - CR3BP w ZH', 'fontsize', 20, 'fontname', 'Times New Roman', 'location', 'northeast')
% title('Hg2 ... ns1 ... Both','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)

end % run_Hg2__ns1
%% ------------------------------------------------
%%% Hg2_2P2__ns1
% -------------------------------------------------
if run_Hg2_2P2__ns1
fprintf('---------------------------------')
fprintf('Hg2_2P2__ns1\n')

[T_Hg2_2P2__ns1, X_Hg2_2P2__ns1, t_ev, x_ev, in_ev] = ode113(@Int_CR3BnSTM, [0, PO_Hg2_2P2__ns1(end)], [PO_Hg2_2P2__ns1(1:6); stm0_colVec], options_apsis, prms);
apses_alts = rowNorm(x_ev(:,1:3) - [1-prms.u,0,0]) - prms.R2;
% figure; hold all
% plot3(X_Hg2_2P2__ns1(:,1),X_Hg2_2P2__ns1(:,2),X_Hg2_2P2__ns1(:,3), 'linewidth', 2, 'color', colors.black)
% plotSecondary(secondary)
% PlotBoi3_CR3Bn(26)
% axis equal
% title('Hg2\_2P2 ... ns1','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)
% 
% 
% plot3(x_ev((apses_alts < 5/rNorm),1), x_ev((apses_alts < 5/rNorm),2), x_ev((apses_alts < 5/rNorm),3), '.','markersize',20, 'color', colors.red2)
% set(gca,'yticklabel',num2str(get(gca,'ytick')','%1.3f'))
% set(gca,'zticklabel',num2str(get(gca,'ztick')','%1.3f'))



[T_Hg2_2P2__ns1_ZH, X_Hg2_2P2__ns1_ZH, t_ev_ZH, x_ev_ZH, in_ev_ZH] = ode113(@Int_CR3BnSTM_J2pJ4pJ6pJ2s, [0, PO_Hg2_2P2__ns1_ZH(end)], [PO_Hg2_2P2__ns1_ZH(1:6); stm0_colVec], options_apsis, prms_ZH);
fprintf('Prop Error = %1.2e\n', norm(X_Hg2_2P2__ns1_ZH(end,1:6) - X_Hg2_2P2__ns1_ZH(1,1:6)) / norm(X_Hg2_2P2__ns1_ZH(1,1:6)))
apses_alts_ZH = rowNorm(x_ev_ZH(:,1:3) - [1-prms_ZH.u,0,0]) - prms_ZH.R2;
close_indices_ZH = ((74/rNorm < apses_alts_ZH) & (apses_alts_ZH < 75/rNorm));
% figure; hold all
% plot3(X_Hg2_2P2__ns1_ZH(:,1),X_Hg2_2P2__ns1_ZH(:,2),X_Hg2_2P2__ns1_ZH(:,3), 'linewidth', 2, 'color', colors.black)
% plotSecondary(secondary)
% PlotBoi3_CR3Bn(26)
% axis equal
% title('Hg2\_2P2 ... ns1 ... ZH','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)
% 
% 
% plot3(x_ev_ZH((apses_alts_ZH < 5/rNorm),1), x_ev_ZH((apses_alts_ZH < 5/rNorm),2), x_ev_ZH((apses_alts_ZH < 5/rNorm),3), '.','markersize',20, 'color', colors.red2)
% set(gca,'yticklabel',num2str(get(gca,'ytick')','%1.3f'))
% set(gca,'zticklabel',num2str(get(gca,'ztick')','%1.3f'))


% figure; hold all
% p_3B = plot3(X_Hg2_2P2__ns1(:,1),X_Hg2_2P2__ns1(:,2),X_Hg2_2P2__ns1(:,3), 'linewidth', 2, 'color', colors.blue2);
% p_ZH = plot3(X_Hg2_2P2__ns1_ZH(:,1),X_Hg2_2P2__ns1_ZH(:,2),X_Hg2_2P2__ns1_ZH(:,3), 'linewidth', 2, 'color', colors.red2);
% plotSecondary(secondary)
% PlotBoi3_CR3Bn(26)
% axis equal
% legend([p_3B, p_ZH], 'CR3BP', 'CR3BP w ZH')
% set(gca,'yticklabel',num2str(get(gca,'ytick')','%1.3f'))
% set(gca,'zticklabel',num2str(get(gca,'ztick')','%1.3f'))


stm_tf_t0                           = reshape(X_Hg2_2P2__ns1(end,7:42),6,6);
monodromy                           = stm_tf_t0;
[eigenVectors_new, eigenValues_new] = eig(monodromy);
[S1, S2]                            = getStabilityIndices(diag(eigenValues_new));

stm_tf_t0_ZH                              = reshape(X_Hg2_2P2__ns1_ZH(end,7:42),6,6);
monodromy_ZH                              = stm_tf_t0_ZH;
[eigenVectors_new_ZH, eigenValues_new_ZH] = eig(monodromy_ZH);
[S1_ZH, S2_ZH]                            = getStabilityIndices(diag(eigenValues_new_ZH));

fprintf('[S1, S2]       = [%1.3f, %1.3f]\n', S1, S2)
fprintf('[S1_ZH, S2_ZH] = [%1.3f, %1.3f]\n', S1_ZH, S2_ZH)

latLons_Hg2_2P2__ns1 = zeros(length(T_Hg2_2P2__ns1),2);
for kk = 1:length(T_Hg2_2P2__ns1)
    [lat_deg, lon_deg] = BCR2latlon(X_Hg2_2P2__ns1(kk,1:3)', 'secondary', prms.u);
    latLons_Hg2_2P2__ns1(kk,:) = [lat_deg, lon_deg];
end
latLons_apses_Hg2_2P2__ns1 = zeros(length(t_ev),2);
for kk = 1:length(t_ev)
    [lat_deg, lon_deg] = BCR2latlon(x_ev(kk,1:3)', 'secondary', prms.u);
    latLons_apses_Hg2_2P2__ns1(kk,:) = [lat_deg, lon_deg];
end


latLons_Hg2_2P2__ns1_ZH = zeros(length(T_Hg2_2P2__ns1_ZH),2);
for kk = 1:length(T_Hg2_2P2__ns1_ZH)
    [lat_deg, lon_deg] = BCR2latlon(X_Hg2_2P2__ns1_ZH(kk,1:3)', 'secondary', prms_ZH.u);
    latLons_Hg2_2P2__ns1_ZH(kk,:) = [lat_deg, lon_deg];
end
latLons_apses_Hg2_2P2__ns1_ZH = zeros(length(t_ev_ZH),2);
for kk = 1:length(t_ev_ZH)
    [lat_deg, lon_deg] = BCR2latlon(x_ev_ZH(kk,1:3)', 'secondary', prms_ZH.u);
    latLons_apses_Hg2_2P2__ns1_ZH(kk,:) = [lat_deg, lon_deg];
end

% figure; hold all
% xlim([0 360])
% ylim([-90 90])
% PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')
% h = image(xlim, -ylim, bodies.europa.img);
% [lons_new] = convert_lon180_to_lon360(latLons_Hg2_2P2__ns1(:,2));
% plot(lons_new,latLons_Hg2_2P2__ns1(:,1),'k.')
% [lons_new_impact] = convert_lon180_to_lon360(latLons_apses_Hg2_2P2__ns1((apses_alts < 5/rNorm),2));
% p_landing = plot(lons_new_impact, latLons_apses_Hg2_2P2__ns1((apses_alts < 5/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.drkred, 'markerfacecolor', colors.red2);
% legend([p_landing], 'Tangent Impact', 'fontsize', 20, 'fontname', 'Times New Roman', 'location', 'northeast')
% title('Hg2\_2P2 ... ns1','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)


figure; hold all
xlim([0 360])
ylim([-90 90])
PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')
h = image(xlim, -ylim, bodies.europa.img);
[lons_new_ZH] = convert_lon180_to_lon360(latLons_Hg2_2P2__ns1_ZH(:,2));
plot(lons_new_ZH,latLons_Hg2_2P2__ns1_ZH(:,1),'.', 'color', colors.red2)
[lons_new_impact_ZH] = convert_lon180_to_lon360(latLons_apses_Hg2_2P2__ns1_ZH((apses_alts_ZH < 5/rNorm),2));
p_landing_ZH = plot(lons_new_impact_ZH, latLons_apses_Hg2_2P2__ns1_ZH((apses_alts_ZH < 5/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.black, 'markerfacecolor', colors.drkgrey);
[lons_new_impact_close] = convert_lon180_to_lon360(latLons_apses_Hg2_2P2__ns1_ZH(close_indices_ZH,2));
p_close = plot(lons_new_impact_close, latLons_apses_Hg2_2P2__ns1_ZH(close_indices_ZH,1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.drkblue, 'markerfacecolor', colors.blue2);
legend([p_landing_ZH, p_close], 'Tangent Impact', '74 km Periapsis', 'fontsize', 20, 'fontname', 'Times New Roman', 'location', 'northeast')
title('Hg2\_2P2 ... ns1 ... ZH','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)


% figure; hold all
% xlim([0 360])
% ylim([-90 90])
% PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')
% h = image(xlim, -ylim, bodies.europa.img);
% plot(lons_new,latLons_Hg2_2P2__ns1(:,1),'.','color',colors.blue2)
% plot(lons_new_ZH,latLons_Hg2_2P2__ns1_ZH(:,1),'.', 'color', colors.red2)
% p_landing = plot(lons_new_impact, latLons_apses_Hg2_2P2__ns1((apses_alts < 5/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.drkblue, 'markerfacecolor', colors.blue2);
% p_landing_ZH = plot(lons_new_impact_ZH, latLons_apses_Hg2_2P2__ns1_ZH((apses_alts_ZH < 5/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.drkred, 'markerfacecolor', colors.red2);
% legend([p_landing, p_landing_ZH], 'Tangent Impact - CR3BP', 'Tangent Impact - CR3BP w ZH', 'fontsize', 20, 'fontname', 'Times New Roman', 'location', 'northeast')
% title('Hg2\_2P2 ... ns1 ... Both','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)


figure; hold all
% plot3(X_Hg2_2P2__ns1_ZH(:,1),X_Hg2_2P2__ns1_ZH(:,2),X_Hg2_2P2__ns1_ZH(:,3), 'linewidth', 2, 'color', colors.red2)
plot3(X_Hg2_2P2__ns1_ZH(:,1),X_Hg2_2P2__ns1_ZH(:,2),X_Hg2_2P2__ns1_ZH(:,3), 'linewidth', 2, 'color', colors.mag)
plotTrajShadows(X_Hg2_2P2__ns1_ZH, 2, colors.grey, 'x', 0.975, 'y', 2.7e-2, 'z', -2.2e-2, 'bodyshadow', [1-prms_ZH.u, prms_ZH.R2])
plotSecondary(secondary)
PlotBoi3_CR3Bn(26)
axis equal
plot3(x_ev_ZH((apses_alts_ZH < 5/rNorm),1), x_ev_ZH((apses_alts_ZH < 5/rNorm),2), x_ev_ZH((apses_alts_ZH < 5/rNorm),3), '.','markersize',20, 'color', colors.black)
plot3(x_ev_ZH(close_indices_ZH, 1), x_ev_ZH(close_indices_ZH, 2), x_ev_ZH(close_indices_ZH, 3), '.','markersize',20, 'color', colors.blue2)
set(gca,'yticklabel',num2str(get(gca,'ytick')','%1.3f'))
set(gca,'zticklabel',num2str(get(gca,'ztick')','%1.3f'))
view(22, 32)


end % run_Hg2_2P2__ns1
%% ------------------------------------------------
%%% Hg2_3T__ns2
% -------------------------------------------------
if run_Hg2_3T__ns2
fprintf('---------------------------------')
fprintf('Hg2_3T__ns2\n')

[T_Hg2_3T__ns2, X_Hg2_3T__ns2, t_ev, x_ev, in_ev] = ode113(@Int_CR3BnSTM, [0, PO_Hg2_3T__ns2(end)], [PO_Hg2_3T__ns2(1:6); stm0_colVec], options_apsis, prms);
apses_alts = rowNorm(x_ev(:,1:3) - [1-prms.u,0,0]) - prms.R2;
% figure; hold all
% plot3(X_Hg2_3T__ns2(:,1),X_Hg2_3T__ns2(:,2),X_Hg2_3T__ns2(:,3), 'linewidth', 2, 'color', colors.black)
% plotSecondary(secondary)
% PlotBoi3_CR3Bn(26)
% axis equal
% title('Hg2\_3T ... ns2','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)
% 
% 
% plot3(x_ev((apses_alts < 8/rNorm),1), x_ev((apses_alts < 8/rNorm),2), x_ev((apses_alts < 8/rNorm),3), '.','markersize',20, 'color', colors.red2)
% set(gca,'yticklabel',num2str(get(gca,'ytick')','%1.3f'))
% set(gca,'zticklabel',num2str(get(gca,'ztick')','%1.3f'))

[T_Hg2_3T__ns2_ZH, X_Hg2_3T__ns2_ZH, t_ev_ZH, x_ev_ZH, in_ev_ZH] = ode113(@Int_CR3BnSTM_J2pJ4pJ6pJ2s, [0, PO_Hg2_3T__ns2_ZH(end)], [PO_Hg2_3T__ns2_ZH(1:6); stm0_colVec], options_apsis, prms_ZH);
apses_alts_ZH = rowNorm(x_ev_ZH(:,1:3) - [1-prms_ZH.u,0,0]) - prms_ZH.R2;
fprintf('Prop Error = %1.2e\n', norm(X_Hg2_3T__ns2_ZH(end,1:6) - X_Hg2_3T__ns2_ZH(1,1:6)) / norm(X_Hg2_3T__ns2_ZH(1,1:6)))
% figure; hold all
% plot3(X_Hg2_3T__ns2_ZH(:,1),X_Hg2_3T__ns2_ZH(:,2),X_Hg2_3T__ns2_ZH(:,3), 'linewidth', 2, 'color', colors.black)
% plotSecondary(secondary)
% PlotBoi3_CR3Bn(26)
% axis equal
% title('Hg2\_3T ... ns2 ... ZH','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)
% 
% plot3(x_ev_ZH((apses_alts_ZH < 8/rNorm),1), x_ev_ZH((apses_alts_ZH < 8/rNorm),2), x_ev_ZH((apses_alts_ZH < 8/rNorm),3), '.','markersize',20, 'color', colors.red2)
% set(gca,'yticklabel',num2str(get(gca,'ytick')','%1.3f'))
% set(gca,'zticklabel',num2str(get(gca,'ztick')','%1.3f'))
% 
% 
% figure; hold all
% p_3B = plot3(X_Hg2_3T__ns2(:,1),X_Hg2_3T__ns2(:,2),X_Hg2_3T__ns2(:,3), 'linewidth', 2, 'color', colors.blue2);
% p_ZH = plot3(X_Hg2_3T__ns2_ZH(:,1),X_Hg2_3T__ns2_ZH(:,2),X_Hg2_3T__ns2_ZH(:,3), 'linewidth', 2, 'color', colors.red2);
% plotSecondary(secondary)
% PlotBoi3_CR3Bn(26)
% axis equal
% legend([p_3B, p_ZH], 'CR3BP', 'CR3BP w ZH')
% set(gca,'yticklabel',num2str(get(gca,'ytick')','%1.3f'))
% set(gca,'zticklabel',num2str(get(gca,'ztick')','%1.3f'))



stm_tf_t0                           = reshape(X_Hg2_3T__ns2(end,7:42),6,6);
monodromy                           = stm_tf_t0;
[eigenVectors_new, eigenValues_new] = eig(monodromy);
[S1, S2]                            = getStabilityIndices(diag(eigenValues_new));

stm_tf_t0_ZH                              = reshape(X_Hg2_3T__ns2_ZH(end,7:42),6,6);
monodromy_ZH                              = stm_tf_t0_ZH;
[eigenVectors_new_ZH, eigenValues_new_ZH] = eig(monodromy_ZH);
[S1_ZH, S2_ZH]                            = getStabilityIndices(diag(eigenValues_new_ZH));

fprintf('[S1, S2]       = [%1.3f, %1.3f]\n', S1, S2)
fprintf('[S1_ZH, S2_ZH] = [%1.3f, %1.3f]\n', S1_ZH, S2_ZH)


latLons_Hg2_3T__ns2 = zeros(length(T_Hg2_3T__ns2),2);
for kk = 1:length(T_Hg2_3T__ns2)
    [lat_deg, lon_deg] = BCR2latlon(X_Hg2_3T__ns2(kk,1:3)', 'secondary', prms.u);
    latLons_Hg2_3T__ns2(kk,:) = [lat_deg, lon_deg];
end
latLons_apses_Hg2_3T__ns2 = zeros(length(t_ev),2);
for kk = 1:length(t_ev)
    [lat_deg, lon_deg] = BCR2latlon(x_ev(kk,1:3)', 'secondary', prms.u);
    latLons_apses_Hg2_3T__ns2(kk,:) = [lat_deg, lon_deg];
end


latLons_Hg2_3T__ns2_ZH = zeros(length(T_Hg2_3T__ns2_ZH),2);
for kk = 1:length(T_Hg2_3T__ns2_ZH)
    [lat_deg, lon_deg] = BCR2latlon(X_Hg2_3T__ns2_ZH(kk,1:3)', 'secondary', prms_ZH.u);
    latLons_Hg2_3T__ns2_ZH(kk,:) = [lat_deg, lon_deg];
end
latLons_apses_Hg2_3T__ns2_ZH = zeros(length(t_ev_ZH),2);
for kk = 1:length(t_ev_ZH)
    [lat_deg, lon_deg] = BCR2latlon(x_ev_ZH(kk,1:3)', 'secondary', prms_ZH.u);
    latLons_apses_Hg2_3T__ns2_ZH(kk,:) = [lat_deg, lon_deg];
end

% figure; hold all
% xlim([0 360])
% ylim([-90 90])
% PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')
% h = image(xlim, -ylim, bodies.europa.img);
% [lons_new] = convert_lon180_to_lon360(latLons_Hg2_3T__ns2(:,2));
% plot(lons_new,latLons_Hg2_3T__ns2(:,1),'k.')
% [lons_new_impact] = convert_lon180_to_lon360(latLons_apses_Hg2_3T__ns2((apses_alts < 8/rNorm),2));
% p_landing = plot(lons_new_impact, latLons_apses_Hg2_3T__ns2((apses_alts < 8/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.drkred, 'markerfacecolor', colors.red2);
% legend([p_landing], 'Tangent Impact', 'fontsize', 20, 'fontname', 'Times New Roman', 'location', 'northeast')
% title('Hg2\_3T ... ns2','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)


figure; hold all
xlim([0 360])
ylim([-90 90])
PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')
h = image(xlim, -ylim, bodies.europa.img);
[lons_new_ZH] = convert_lon180_to_lon360(latLons_Hg2_3T__ns2_ZH(:,2));
plot(lons_new_ZH,latLons_Hg2_3T__ns2_ZH(:,1),'.', 'color', colors.red2)
[lons_new_impact_ZH] = convert_lon180_to_lon360(latLons_apses_Hg2_3T__ns2_ZH((apses_alts_ZH < 8/rNorm),2));
p_landing_ZH = plot(lons_new_impact_ZH, latLons_apses_Hg2_3T__ns2_ZH((apses_alts_ZH < 8/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.black, 'markerfacecolor', colors.drkgrey);
legend([p_landing_ZH], 'Tangent Impact', 'fontsize', 20, 'fontname', 'Times New Roman', 'location', 'northeast')
title('Hg2\_3T ... ns2 ... ZH','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)

% figure; hold all
% xlim([0 360])
% ylim([-90 90])
% PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')
% h = image(xlim, -ylim, bodies.europa.img);
% plot(lons_new,latLons_Hg2_3T__ns2(:,1),'.','color',colors.blue2)
% plot(lons_new_ZH,latLons_Hg2_3T__ns2_ZH(:,1),'.', 'color', colors.red2)
% p_landing = plot(lons_new_impact, latLons_apses_Hg2_3T__ns2((apses_alts < 8/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.drkblue, 'markerfacecolor', colors.blue2);
% p_landing_ZH = plot(lons_new_impact_ZH, latLons_apses_Hg2_3T__ns2_ZH((apses_alts_ZH < 8/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.drkred, 'markerfacecolor', colors.red2);
% legend([p_landing, p_landing_ZH], 'Tangent Impact - CR3BP', 'Tangent Impact - CR3BP w ZH', 'fontsize', 20, 'fontname', 'Times New Roman', 'location', 'northeast')
% title('Hg2\_3T ... ns2 ... Both','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)


figure; hold all
% plot3(X_Hg2_3T__ns2_ZH(:,1),X_Hg2_3T__ns2_ZH(:,2),X_Hg2_3T__ns2_ZH(:,3), 'linewidth', 2, 'color', colors.red2)
plot3(X_Hg2_3T__ns2_ZH(:,1),X_Hg2_3T__ns2_ZH(:,2),X_Hg2_3T__ns2_ZH(:,3), 'linewidth', 2, 'color', colors.mag)
plotTrajShadows(X_Hg2_3T__ns2_ZH, 2, colors.grey, 'x', 1.02, 'y', 2.6e-2, 'z', -2.2e-2, 'bodyshadow', [1-prms_ZH.u, prms_ZH.R2])
plotSecondary(secondary)
PlotBoi3_CR3Bn(26)
axis equal
plot3(x_ev_ZH((apses_alts_ZH < 8/rNorm),1), x_ev_ZH((apses_alts_ZH < 8/rNorm),2), x_ev_ZH((apses_alts_ZH < 8/rNorm),3), '.','markersize',20, 'color', colors.black)
set(gca,'yticklabel',num2str(get(gca,'ytick')','%1.3f'))
set(gca,'zticklabel',num2str(get(gca,'ztick')','%1.3f'))




end % run_Hg2_3T__ns2
%% ------------------------------------------------
%%% Se7__s1
% -------------------------------------------------
if run_Se7__s1
fprintf('---------------------------------')
fprintf('Se7__s1 (planar)\n')

[T_Se7__s1, X_Se7__s1, t_ev, x_ev, in_ev] = ode113(@Int_CR3BnSTM, [0, PO_Se7__s1(end)], [PO_Se7__s1(1:6); stm0_colVec], options_apsis, prms);
apses_alts = rowNorm(x_ev(:,1:3) - [1-prms.u,0,0]) - prms.R2;
% figure; hold all
% plot3(X_Se7__s1(:,1),X_Se7__s1(:,2),X_Se7__s1(:,3), 'linewidth', 2, 'color', colors.black)
% plotSecondary(secondary)
% PlotBoi3_CR3Bn(26)
% axis equal
% title('Se7 ... s1','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)
% 
% plot3(x_ev((apses_alts < 5/rNorm),1), x_ev((apses_alts < 5/rNorm),2), x_ev((apses_alts < 5/rNorm),3), '.','markersize',20, 'color', colors.red2)
% set(gca,'yticklabel',num2str(get(gca,'ytick')','%1.3f'))
% set(gca,'zticklabel',num2str(get(gca,'ztick')','%1.3f'))


[T_Se7__s1_ZH, X_Se7__s1_ZH, t_ev_ZH, x_ev_ZH, in_ev_ZH] = ode113(@Int_CR3BnSTM_J2pJ4pJ6pJ2s, [0, PO_Se7__s1_ZH(end)], [PO_Se7__s1_ZH(1:6); stm0_colVec], options_apsis, prms_ZH);
fprintf('Prop Error = %1.2e\n', norm(X_Se7__s1_ZH(end,1:6) - X_Se7__s1_ZH(1,1:6)) / norm(X_Se7__s1_ZH(1,1:6)))
figure; hold all
% plot3(X_Se7__s1_ZH(:,1),X_Se7__s1_ZH(:,2),X_Se7__s1_ZH(:,3), 'linewidth', 2, 'color', colors.red2)
plot3(X_Se7__s1_ZH(:,1),X_Se7__s1_ZH(:,2),X_Se7__s1_ZH(:,3), 'linewidth', 2, 'color', colors.mag)
plotSecondary(secondary)
PlotBoi3_CR3Bn(26)
axis equal
title('Se7 ... s1 ... ZH','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)
set(gca,'yticklabel',num2str(get(gca,'ytick')','%1.3f'))
set(gca,'zticklabel',num2str(get(gca,'ztick')','%1.3f'))

apses_alts_ZH = rowNorm(x_ev_ZH(:,1:3) - [1-prms_ZH.u,0,0]) - prms_ZH.R2;
plot3(x_ev_ZH((apses_alts_ZH < 5/rNorm),1), x_ev_ZH((apses_alts_ZH < 5/rNorm),2), x_ev_ZH((apses_alts_ZH < 5/rNorm),3), '.','markersize',20, 'color', colors.black)

close_indices_ZH = ((110/rNorm < apses_alts_ZH) & (apses_alts_ZH < 111/rNorm));
plot3(x_ev_ZH(close_indices_ZH, 1), x_ev_ZH(close_indices_ZH, 2), x_ev_ZH(close_indices_ZH, 3), '.','markersize',20, 'color', colors.blue2)

view(0,90)

% figure; hold all
% p_3B = plot3(X_Se7__s1(:,1),X_Se7__s1(:,2),X_Se7__s1(:,3), 'linewidth', 2, 'color', colors.blue2);
% p_ZH = plot3(X_Se7__s1_ZH(:,1),X_Se7__s1_ZH(:,2),X_Se7__s1_ZH(:,3), 'linewidth', 2, 'color', colors.red2);
% plotSecondary(secondary)
% PlotBoi3_CR3Bn(26)
% axis equal
% legend([p_3B, p_ZH], 'CR3BP', 'CR3BP w ZH')
% set(gca,'yticklabel',num2str(get(gca,'ytick')','%1.3f'))
% set(gca,'zticklabel',num2str(get(gca,'ztick')','%1.3f'))



stm_tf_t0                           = reshape(X_Se7__s1(end,7:42),6,6);
monodromy                           = stm_tf_t0;
[eigenVectors_new, eigenValues_new] = eig(monodromy);
[S1, S2]                            = getStabilityIndices(diag(eigenValues_new));

stm_tf_t0_ZH                              = reshape(X_Se7__s1_ZH(end,7:42),6,6);
monodromy_ZH                              = stm_tf_t0_ZH;
[eigenVectors_new_ZH, eigenValues_new_ZH] = eig(monodromy_ZH);
[S1_ZH, S2_ZH]                            = getStabilityIndices(diag(eigenValues_new_ZH));

fprintf('[S1, S2]       = [%1.3f, %1.3f]\n', S1, S2)
fprintf('[S1_ZH, S2_ZH] = [%1.3f, %1.3f]\n', S1_ZH, S2_ZH)


latLons_Se7__s1 = zeros(length(T_Se7__s1),2);
for kk = 1:length(T_Se7__s1)
    [lat_deg, lon_deg] = BCR2latlon(X_Se7__s1(kk,1:3)', 'secondary', prms.u);
    latLons_Se7__s1(kk,:) = [lat_deg, lon_deg];
end
latLons_apses_Se7__s1 = zeros(length(t_ev),2);
for kk = 1:length(t_ev)
    [lat_deg, lon_deg] = BCR2latlon(x_ev(kk,1:3)', 'secondary', prms.u);
    latLons_apses_Se7__s1(kk,:) = [lat_deg, lon_deg];
end


latLons_Se7__s1_ZH = zeros(length(T_Se7__s1_ZH),2);
for kk = 1:length(T_Se7__s1_ZH)
    [lat_deg, lon_deg] = BCR2latlon(X_Se7__s1_ZH(kk,1:3)', 'secondary', prms_ZH.u);
    latLons_Se7__s1_ZH(kk,:) = [lat_deg, lon_deg];
end
latLons_apses_Se7__s1_ZH = zeros(length(t_ev_ZH),2);
for kk = 1:length(t_ev_ZH)
    [lat_deg, lon_deg] = BCR2latlon(x_ev_ZH(kk,1:3)', 'secondary', prms_ZH.u);
    latLons_apses_Se7__s1_ZH(kk,:) = [lat_deg, lon_deg];
end

% figure; hold all
% xlim([0 360])
% ylim([-90 90])
% PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')
% h = image(xlim, -ylim, bodies.europa.img);
% [lons_new] = convert_lon180_to_lon360(latLons_Se7__s1(:,2));
% plot(lons_new,latLons_Se7__s1(:,1),'k.')
% [lons_new_impact] = convert_lon180_to_lon360(latLons_apses_Se7__s1((apses_alts < 5/rNorm),2));
% p_landing = plot(lons_new_impact, latLons_apses_Se7__s1((apses_alts < 5/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.drkred, 'markerfacecolor', colors.red2);
% legend([p_landing], 'Tangent Impact', 'fontsize', 20, 'fontname', 'Times New Roman', 'location', 'northeast')
% title('Se7 ... s1','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)


figure; hold all
xlim([0 360])
ylim([-90 90])
PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')
h = image(xlim, -ylim, bodies.europa.img);
[lons_new_ZH] = convert_lon180_to_lon360(latLons_Se7__s1_ZH(:,2));
plot(lons_new_ZH,latLons_Se7__s1_ZH(:,1),'.', 'color', colors.red2)
[lons_new_impact_ZH] = convert_lon180_to_lon360(latLons_apses_Se7__s1_ZH((apses_alts_ZH < 5/rNorm),2));
p_landing_ZH = plot(lons_new_impact_ZH, latLons_apses_Se7__s1_ZH((apses_alts_ZH < 5/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.black, 'markerfacecolor', colors.drkgrey);
[lons_new_impact_close] = convert_lon180_to_lon360(latLons_apses_Se7__s1_ZH(close_indices_ZH,2));
p_close = plot(lons_new_impact_close, latLons_apses_Se7__s1_ZH(close_indices_ZH,1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.drkblue, 'markerfacecolor', colors.blue2);
legend([p_landing_ZH, p_close], 'Tangent Impact', '110 km Periapsis', 'fontsize', 20, 'fontname', 'Times New Roman', 'location', 'northeast')
title('Se7 ... s1 ... ZH','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)



% figure; hold all
% xlim([0 360])
% ylim([-90 90])
% PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')
% h = image(xlim, -ylim, bodies.europa.img);
% plot(lons_new,latLons_Se7__s1(:,1),'.','color',colors.blue2)
% plot(lons_new_ZH,latLons_Se7__s1_ZH(:,1),'.', 'color', colors.red2)
% p_landing = plot(lons_new_impact, latLons_apses_Se7__s1((apses_alts < 5/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.drkblue, 'markerfacecolor', colors.blue2);
% p_landing_ZH = plot(lons_new_impact_ZH, latLons_apses_Se7__s1_ZH((apses_alts_ZH < 5/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.drkred, 'markerfacecolor', colors.red2);
% legend([p_landing, p_landing_ZH], 'Tangent Impact - CR3BP', 'Tangent Impact - CR3BP w ZH', 'fontsize', 20, 'fontname', 'Times New Roman', 'location', 'northeast')
% title('Se7 ... s1 ... Both','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)

end % run_Se7__s1
%% ------------------------------------------------
%%% Se7__ns1
% -------------------------------------------------
if run_Se7__ns1
fprintf('---------------------------------')
fprintf('Se7__ns1 (planar)\n')

[T_Se7__ns1, X_Se7__ns1, t_ev, x_ev, in_ev] = ode113(@Int_CR3BnSTM, [0, PO_Se7__ns1(end)], [PO_Se7__ns1(1:6); stm0_colVec], options_apsis, prms);
apses_alts = rowNorm(x_ev(:,1:3) - [1-prms.u,0,0]) - prms.R2;
% figure; hold all
% plot3(X_Se7__ns1(:,1),X_Se7__ns1(:,2),X_Se7__ns1(:,3), 'linewidth', 2, 'color', colors.black)
% plotSecondary(secondary)
% PlotBoi3_CR3Bn(26)
% axis equal
% title('Se7 ... ns1','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)
% 
% plot3(x_ev((apses_alts < 5/rNorm),1), x_ev((apses_alts < 5/rNorm),2), x_ev((apses_alts < 5/rNorm),3), '.','markersize',20, 'color', colors.red2)
% set(gca,'yticklabel',num2str(get(gca,'ytick')','%1.3f'))
% set(gca,'zticklabel',num2str(get(gca,'ztick')','%1.3f'))


[T_Se7__ns1_ZH, X_Se7__ns1_ZH, t_ev_ZH, x_ev_ZH, in_ev_ZH] = ode113(@Int_CR3BnSTM_J2pJ4pJ6pJ2s, [0, PO_Se7__ns1_ZH(end)], [PO_Se7__ns1_ZH(1:6); stm0_colVec], options_apsis, prms_ZH);
fprintf('Prop Error = %1.2e\n', norm(X_Se7__ns1_ZH(end,1:6) - X_Se7__ns1_ZH(1,1:6)) / norm(X_Se7__ns1_ZH(1,1:6)))
figure; hold all
% plot3(X_Se7__ns1_ZH(:,1),X_Se7__ns1_ZH(:,2),X_Se7__ns1_ZH(:,3), 'linewidth', 2, 'color', colors.red2)
plot3(X_Se7__ns1_ZH(:,1),X_Se7__ns1_ZH(:,2),X_Se7__ns1_ZH(:,3), 'linewidth', 2, 'color', colors.mag)
plotSecondary(secondary)
PlotBoi3_CR3Bn(26)
axis equal
title('Se7 ... ns1 ... ZH','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)

apses_alts_ZH = rowNorm(x_ev_ZH(:,1:3) - [1-prms_ZH.u,0,0]) - prms_ZH.R2;
plot3(x_ev_ZH((apses_alts_ZH < 5/rNorm),1), x_ev_ZH((apses_alts_ZH < 5/rNorm),2), x_ev_ZH((apses_alts_ZH < 5/rNorm),3), '.','markersize',20, 'color', colors.black)
set(gca,'yticklabel',num2str(get(gca,'ytick')','%1.3f'))
set(gca,'zticklabel',num2str(get(gca,'ztick')','%1.3f'))

close_indices_ZH = ((43/rNorm < apses_alts_ZH) & (apses_alts_ZH < 44/rNorm));
plot3(x_ev_ZH(close_indices_ZH, 1), x_ev_ZH(close_indices_ZH, 2), x_ev_ZH(close_indices_ZH, 3), '.','markersize',20, 'color', colors.blue2)

view(0,90)

% figure; hold all
% p_3B = plot3(X_Se7__ns1(:,1),X_Se7__ns1(:,2),X_Se7__ns1(:,3), 'linewidth', 2, 'color', colors.blue2);
% p_ZH = plot3(X_Se7__ns1_ZH(:,1),X_Se7__ns1_ZH(:,2),X_Se7__ns1_ZH(:,3), 'linewidth', 2, 'color', colors.red2);
% plotSecondary(secondary)
% PlotBoi3_CR3Bn(26)
% axis equal
% legend([p_3B, p_ZH], 'CR3BP', 'CR3BP w ZH')
% set(gca,'yticklabel',num2str(get(gca,'ytick')','%1.3f'))
% set(gca,'zticklabel',num2str(get(gca,'ztick')','%1.3f'))



stm_tf_t0                           = reshape(X_Se7__ns1(end,7:42),6,6);
monodromy                           = stm_tf_t0;
[eigenVectors_new, eigenValues_new] = eig(monodromy);
[S1, S2]                            = getStabilityIndices(diag(eigenValues_new));

stm_tf_t0_ZH                              = reshape(X_Se7__ns1_ZH(end,7:42),6,6);
monodromy_ZH                              = stm_tf_t0_ZH;
[eigenVectors_new_ZH, eigenValues_new_ZH] = eig(monodromy_ZH);
[S1_ZH, S2_ZH]                            = getStabilityIndices(diag(eigenValues_new_ZH));

fprintf('[S1, S2]       = [%1.3f, %1.3f]\n', S1, S2)
fprintf('[S1_ZH, S2_ZH] = [%1.3f, %1.3f]\n', S1_ZH, S2_ZH)


latLons_Se7__ns1 = zeros(length(T_Se7__ns1),2);
for kk = 1:length(T_Se7__ns1)
    [lat_deg, lon_deg] = BCR2latlon(X_Se7__ns1(kk,1:3)', 'secondary', prms.u);
    latLons_Se7__ns1(kk,:) = [lat_deg, lon_deg];
end
latLons_apses_Se7__ns1 = zeros(length(t_ev),2);
for kk = 1:length(t_ev)
    [lat_deg, lon_deg] = BCR2latlon(x_ev(kk,1:3)', 'secondary', prms.u);
    latLons_apses_Se7__ns1(kk,:) = [lat_deg, lon_deg];
end


latLons_Se7__ns1_ZH = zeros(length(T_Se7__ns1_ZH),2);
for kk = 1:length(T_Se7__ns1_ZH)
    [lat_deg, lon_deg] = BCR2latlon(X_Se7__ns1_ZH(kk,1:3)', 'secondary', prms_ZH.u);
    latLons_Se7__ns1_ZH(kk,:) = [lat_deg, lon_deg];
end
latLons_apses_Se7__ns1_ZH = zeros(length(t_ev_ZH),2);
for kk = 1:length(t_ev_ZH)
    [lat_deg, lon_deg] = BCR2latlon(x_ev_ZH(kk,1:3)', 'secondary', prms_ZH.u);
    latLons_apses_Se7__ns1_ZH(kk,:) = [lat_deg, lon_deg];
end

% figure; hold all
% xlim([0 360])
% ylim([-90 90])
% PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')
% h = image(xlim, -ylim, bodies.europa.img);
% [lons_new] = convert_lon180_to_lon360(latLons_Se7__ns1(:,2));
% plot(lons_new,latLons_Se7__ns1(:,1),'k.')
% [lons_new_impact] = convert_lon180_to_lon360(latLons_apses_Se7__ns1((apses_alts < 5/rNorm),2));
% p_landing = plot(lons_new_impact, latLons_apses_Se7__ns1((apses_alts < 5/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.drkred, 'markerfacecolor', colors.red2);
% legend([p_landing], 'Tangent Impact', 'fontsize', 20, 'fontname', 'Times New Roman', 'location', 'northeast')
% title('Se7 ... ns1','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)


figure; hold all
xlim([0 360])
ylim([-90 90])
PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')
h = image(xlim, -ylim, bodies.europa.img);
[lons_new_ZH] = convert_lon180_to_lon360(latLons_Se7__ns1_ZH(:,2));
plot(lons_new_ZH,latLons_Se7__ns1_ZH(:,1),'.', 'color', colors.red2)
[lons_new_impact_ZH] = convert_lon180_to_lon360(latLons_apses_Se7__ns1_ZH((apses_alts_ZH < 5/rNorm),2));
p_landing_ZH = plot(lons_new_impact_ZH, latLons_apses_Se7__ns1_ZH((apses_alts_ZH < 5/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.black, 'markerfacecolor', colors.drkgrey);
[lons_new_impact_close] = convert_lon180_to_lon360(latLons_apses_Se7__ns1_ZH(close_indices_ZH,2));
p_close = plot(lons_new_impact_close, latLons_apses_Se7__ns1_ZH(close_indices_ZH,1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.drkblue, 'markerfacecolor', colors.blue2);
legend([p_landing_ZH, p_close], 'Tangent Impact', '43 km Periapsis', 'fontsize', 20, 'fontname', 'Times New Roman', 'location', 'northeast')
title('Se7 ... ns1 ... ZH','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)



% figure; hold all
% xlim([0 360])
% ylim([-90 90])
% PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')
% h = image(xlim, -ylim, bodies.europa.img);
% plot(lons_new,latLons_Se7__ns1(:,1),'.','color',colors.blue2)
% plot(lons_new_ZH,latLons_Se7__ns1_ZH(:,1),'.', 'color', colors.red2)
% p_landing = plot(lons_new_impact, latLons_apses_Se7__ns1((apses_alts < 5/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.drkblue, 'markerfacecolor', colors.blue2);
% p_landing_ZH = plot(lons_new_impact_ZH, latLons_apses_Se7__ns1_ZH((apses_alts_ZH < 5/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.drkred, 'markerfacecolor', colors.red2);
% legend([p_landing, p_landing_ZH], 'Tangent Impact - CR3BP', 'Tangent Impact - CR3BP w ZH', 'fontsize', 20, 'fontname', 'Times New Roman', 'location', 'northeast')
% title('Se7 ... ns1 ... Both','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)

end % run_Se7__ns1
%% ------------------------------------------------
%%% LoPO__s1
% -------------------------------------------------
if run_LoPO__s1
fprintf('---------------------------------')
fprintf('LoPO__s1 (planar)\n')

[T_LoPO__s1, X_LoPO__s1, t_ev, x_ev, in_ev] = ode113(@Int_CR3BnSTM, [0, PO_LoPO__s1(end)], [PO_LoPO__s1(1:6); stm0_colVec], options_apsis, prms);
apses_alts = rowNorm(x_ev(:,1:3) - [1-prms.u,0,0]) - prms.R2;
% figure; hold all
% plot3(X_LoPO__s1(:,1),X_LoPO__s1(:,2),X_LoPO__s1(:,3), 'linewidth', 2, 'color', colors.black)
% plotSecondary(secondary)
% PlotBoi3_CR3Bn(26)
% axis equal
% title('LoPO ... s1','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)
% 
% plot3(x_ev((apses_alts < 5/rNorm),1), x_ev((apses_alts < 5/rNorm),2), x_ev((apses_alts < 5/rNorm),3), '.','markersize',20, 'color', colors.red2)
% set(gca,'yticklabel',num2str(get(gca,'ytick')','%1.3f'))
% set(gca,'zticklabel',num2str(get(gca,'ztick')','%1.3f'))


[T_LoPO__s1_ZH, X_LoPO__s1_ZH, t_ev_ZH, x_ev_ZH, in_ev_ZH] = ode113(@Int_CR3BnSTM_J2pJ4pJ6pJ2s, [0, PO_LoPO__s1_ZH(end)], [PO_LoPO__s1_ZH(1:6); stm0_colVec], options_apsis, prms_ZH);
fprintf('Prop Error = %1.2e\n', norm(X_LoPO__s1_ZH(end,1:6) - X_LoPO__s1_ZH(1,1:6)) / norm(X_LoPO__s1_ZH(1,1:6)))
figure; hold all
% plot3(X_LoPO__s1_ZH(:,1),X_LoPO__s1_ZH(:,2),X_LoPO__s1_ZH(:,3), 'linewidth', 2, 'color', colors.red2)
plot3(X_LoPO__s1_ZH(:,1),X_LoPO__s1_ZH(:,2),X_LoPO__s1_ZH(:,3), 'linewidth', 2, 'color', colors.mag)
plotSecondary(secondary)
PlotBoi3_CR3Bn(26)
axis equal
title('LoPO ... s1 ... ZH','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)

apses_alts_ZH = rowNorm(x_ev_ZH(:,1:3) - [1-prms_ZH.u,0,0]) - prms_ZH.R2;
plot3(x_ev_ZH((apses_alts_ZH < 5/rNorm),1), x_ev_ZH((apses_alts_ZH < 5/rNorm),2), x_ev_ZH((apses_alts_ZH < 5/rNorm),3), '.','markersize',20, 'color', colors.black)
set(gca,'yticklabel',num2str(get(gca,'ytick')','%1.3f'))
set(gca,'zticklabel',num2str(get(gca,'ztick')','%1.3f'))

view(0,90)

% figure; hold all
% p_3B = plot3(X_LoPO__s1(:,1),X_LoPO__s1(:,2),X_LoPO__s1(:,3), 'linewidth', 2, 'color', colors.blue2);
% p_ZH = plot3(X_LoPO__s1_ZH(:,1),X_LoPO__s1_ZH(:,2),X_LoPO__s1_ZH(:,3), 'linewidth', 2, 'color', colors.red2);
% plotSecondary(secondary)
% PlotBoi3_CR3Bn(26)
% axis equal
% legend([p_3B, p_ZH], 'CR3BP', 'CR3BP w ZH')
% set(gca,'yticklabel',num2str(get(gca,'ytick')','%1.3f'))
% set(gca,'zticklabel',num2str(get(gca,'ztick')','%1.3f'))



stm_tf_t0                           = reshape(X_LoPO__s1(end,7:42),6,6);
monodromy                           = stm_tf_t0;
[eigenVectors_new, eigenValues_new] = eig(monodromy);
[S1, S2]                            = getStabilityIndices(diag(eigenValues_new));

stm_tf_t0_ZH                              = reshape(X_LoPO__s1_ZH(end,7:42),6,6);
monodromy_ZH                              = stm_tf_t0_ZH;
[eigenVectors_new_ZH, eigenValues_new_ZH] = eig(monodromy_ZH);
[S1_ZH, S2_ZH]                            = getStabilityIndices(diag(eigenValues_new_ZH));

fprintf('[S1, S2]       = [%1.3f, %1.3f]\n', S1, S2)
fprintf('[S1_ZH, S2_ZH] = [%1.3f, %1.3f]\n', S1_ZH, S2_ZH)


latLons_LoPO__s1 = zeros(length(T_LoPO__s1),2);
for kk = 1:length(T_LoPO__s1)
    [lat_deg, lon_deg] = BCR2latlon(X_LoPO__s1(kk,1:3)', 'secondary', prms.u);
    latLons_LoPO__s1(kk,:) = [lat_deg, lon_deg];
end
latLons_apses_LoPO__s1 = zeros(length(t_ev),2);
for kk = 1:length(t_ev)
    [lat_deg, lon_deg] = BCR2latlon(x_ev(kk,1:3)', 'secondary', prms.u);
    latLons_apses_LoPO__s1(kk,:) = [lat_deg, lon_deg];
end


latLons_LoPO__s1_ZH = zeros(length(T_LoPO__s1_ZH),2);
for kk = 1:length(T_LoPO__s1_ZH)
    [lat_deg, lon_deg] = BCR2latlon(X_LoPO__s1_ZH(kk,1:3)', 'secondary', prms_ZH.u);
    latLons_LoPO__s1_ZH(kk,:) = [lat_deg, lon_deg];
end
latLons_apses_LoPO__s1_ZH = zeros(length(t_ev_ZH),2);
for kk = 1:length(t_ev_ZH)
    [lat_deg, lon_deg] = BCR2latlon(x_ev_ZH(kk,1:3)', 'secondary', prms_ZH.u);
    latLons_apses_LoPO__s1_ZH(kk,:) = [lat_deg, lon_deg];
end

% figure; hold all
% xlim([0 360])
% ylim([-90 90])
% PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')
% h = image(xlim, -ylim, bodies.europa.img);
% [lons_new] = convert_lon180_to_lon360(latLons_LoPO__s1(:,2));
% plot(lons_new,latLons_LoPO__s1(:,1),'k.')
% [lons_new_impact] = convert_lon180_to_lon360(latLons_apses_LoPO__s1((apses_alts < 5/rNorm),2));
% p_landing = plot(lons_new_impact, latLons_apses_LoPO__s1((apses_alts < 5/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.drkred, 'markerfacecolor', colors.red2);
% legend([p_landing], 'Tangent Impact', 'fontsize', 20, 'fontname', 'Times New Roman', 'location', 'northeast')
% title('LoPO ... s1','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)


figure; hold all
xlim([0 360])
ylim([-90 90])
PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')
h = image(xlim, -ylim, bodies.europa.img);
[lons_new_ZH] = convert_lon180_to_lon360(latLons_LoPO__s1_ZH(:,2));
plot(lons_new_ZH,latLons_LoPO__s1_ZH(:,1),'.', 'color', colors.red2)
[lons_new_impact_ZH] = convert_lon180_to_lon360(latLons_apses_LoPO__s1_ZH((apses_alts_ZH < 5/rNorm),2));
p_landing_ZH = plot(lons_new_impact_ZH, latLons_apses_LoPO__s1_ZH((apses_alts_ZH < 5/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.black, 'markerfacecolor', colors.drkgrey);
legend([p_landing_ZH], 'Tangent Impact', 'fontsize', 20, 'fontname', 'Times New Roman', 'location', 'northeast')
title('LoPO ... s1 ... ZH','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)

% figure; hold all
% xlim([0 360])
% ylim([-90 90])
% PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')
% h = image(xlim, -ylim, bodies.europa.img);
% plot(lons_new,latLons_LoPO__s1(:,1),'.','color',colors.blue2)
% plot(lons_new_ZH,latLons_LoPO__s1_ZH(:,1),'.', 'color', colors.red2)
% p_landing = plot(lons_new_impact, latLons_apses_LoPO__s1((apses_alts < 5/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.drkblue, 'markerfacecolor', colors.blue2);
% p_landing_ZH = plot(lons_new_impact_ZH, latLons_apses_LoPO__s1_ZH((apses_alts_ZH < 5/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.drkred, 'markerfacecolor', colors.red2);
% legend([p_landing, p_landing_ZH], 'Tangent Impact - CR3BP', 'Tangent Impact - CR3BP w ZH', 'fontsize', 20, 'fontname', 'Times New Roman', 'location', 'northeast')
% title('LoPO ... s1 ... Both','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)

end % run_LoPO__s1
%% ------------------------------------------------
%%% LoPO_2P2_1T__ns1
% -------------------------------------------------
if run_LoPO_2P2_1T__ns1
fprintf('---------------------------------')
fprintf('LoPO_2P2_1T__ns1\n')

[T_LoPO_2P2_1T__ns1, X_LoPO_2P2_1T__ns1, t_ev, x_ev, in_ev] = ode113(@Int_CR3BnSTM, [0, PO_LoPO_2P2_1T__ns1(end)], [PO_LoPO_2P2_1T__ns1(1:6); stm0_colVec], options_apsis, prms);
apses_alts = rowNorm(x_ev(:,1:3) - [1-prms.u,0,0]) - prms.R2;
% figure; hold all
% plot3(X_LoPO_2P2_1T__ns1(:,1),X_LoPO_2P2_1T__ns1(:,2),X_LoPO_2P2_1T__ns1(:,3), 'linewidth', 2, 'color', colors.black)
% plotSecondary(secondary)
% PlotBoi3_CR3Bn(26)
% axis equal
% title('LoPO\_2P2\_1T ... ns1','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)
% 
% plot3(x_ev((apses_alts < 5/rNorm),1), x_ev((apses_alts < 5/rNorm),2), x_ev((apses_alts < 5/rNorm),3), '.','markersize',20, 'color', colors.red2)
% set(gca,'yticklabel',num2str(get(gca,'ytick')','%1.3f'))
% set(gca,'zticklabel',num2str(get(gca,'ztick')','%1.3f'))


[T_LoPO_2P2_1T__ns1_ZH, X_LoPO_2P2_1T__ns1_ZH, t_ev_ZH, x_ev_ZH, in_ev_ZH] = ode113(@Int_CR3BnSTM_J2pJ4pJ6pJ2s, [0, PO_LoPO_2P2_1T__ns1_ZH(end)], [PO_LoPO_2P2_1T__ns1_ZH(1:6); stm0_colVec], options_apsis, prms_ZH);
apses_alts_ZH = rowNorm(x_ev_ZH(:,1:3) - [1-prms_ZH.u,0,0]) - prms_ZH.R2;
fprintf('Prop Error = %1.2e\n', norm(X_LoPO_2P2_1T__ns1_ZH(end,1:6) - X_LoPO_2P2_1T__ns1_ZH(1,1:6)) / norm(X_LoPO_2P2_1T__ns1_ZH(1,1:6)))
% figure; hold all
% plot3(X_LoPO_2P2_1T__ns1_ZH(:,1),X_LoPO_2P2_1T__ns1_ZH(:,2),X_LoPO_2P2_1T__ns1_ZH(:,3), 'linewidth', 2, 'color', colors.black)
% plotSecondary(secondary)
% PlotBoi3_CR3Bn(26)
% axis equal
% title('LoPO\_2P2\_1T ... ns1 ... ZH','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)
% 
% plot3(x_ev_ZH((apses_alts_ZH < 5/rNorm),1), x_ev_ZH((apses_alts_ZH < 5/rNorm),2), x_ev_ZH((apses_alts_ZH < 5/rNorm),3), '.','markersize',20, 'color', colors.red2)
% set(gca,'yticklabel',num2str(get(gca,'ytick')','%1.3f'))
% set(gca,'zticklabel',num2str(get(gca,'ztick')','%1.3f'))



% figure; hold all
% p_3B = plot3(X_LoPO_2P2_1T__ns1(:,1),X_LoPO_2P2_1T__ns1(:,2),X_LoPO_2P2_1T__ns1(:,3), 'linewidth', 2, 'color', colors.blue2);
% p_ZH = plot3(X_LoPO_2P2_1T__ns1_ZH(:,1),X_LoPO_2P2_1T__ns1_ZH(:,2),X_LoPO_2P2_1T__ns1_ZH(:,3), 'linewidth', 2, 'color', colors.red2);
% plotSecondary(secondary)
% PlotBoi3_CR3Bn(26)
% axis equal
% legend([p_3B, p_ZH], 'CR3BP', 'CR3BP w ZH')
% set(gca,'yticklabel',num2str(get(gca,'ytick')','%1.3f'))
% set(gca,'zticklabel',num2str(get(gca,'ztick')','%1.3f'))



stm_tf_t0                           = reshape(X_LoPO_2P2_1T__ns1(end,7:42),6,6);
monodromy                           = stm_tf_t0;
[eigenVectors_new, eigenValues_new] = eig(monodromy);
[S1, S2]                            = getStabilityIndices(diag(eigenValues_new));

stm_tf_t0_ZH                              = reshape(X_LoPO_2P2_1T__ns1_ZH(end,7:42),6,6);
monodromy_ZH                              = stm_tf_t0_ZH;
[eigenVectors_new_ZH, eigenValues_new_ZH] = eig(monodromy_ZH);
[S1_ZH, S2_ZH]                            = getStabilityIndices(diag(eigenValues_new_ZH));

fprintf('[S1, S2]       = [%1.3f, %1.3f]\n', S1, S2)
fprintf('[S1_ZH, S2_ZH] = [%1.3f, %1.3f]\n', S1_ZH, S2_ZH)


latLons_LoPO_2P2_1T__ns1 = zeros(length(T_LoPO_2P2_1T__ns1),2);
for kk = 1:length(T_LoPO_2P2_1T__ns1)
    [lat_deg, lon_deg] = BCR2latlon(X_LoPO_2P2_1T__ns1(kk,1:3)', 'secondary', prms.u);
    latLons_LoPO_2P2_1T__ns1(kk,:) = [lat_deg, lon_deg];
end
latLons_apses_LoPO_2P2_1T__ns1 = zeros(length(t_ev),2);
for kk = 1:length(t_ev)
    [lat_deg, lon_deg] = BCR2latlon(x_ev(kk,1:3)', 'secondary', prms.u);
    latLons_apses_LoPO_2P2_1T__ns1(kk,:) = [lat_deg, lon_deg];
end


latLons_LoPO_2P2_1T__ns1_ZH = zeros(length(T_LoPO_2P2_1T__ns1_ZH),2);
for kk = 1:length(T_LoPO_2P2_1T__ns1_ZH)
    [lat_deg, lon_deg] = BCR2latlon(X_LoPO_2P2_1T__ns1_ZH(kk,1:3)', 'secondary', prms_ZH.u);
    latLons_LoPO_2P2_1T__ns1_ZH(kk,:) = [lat_deg, lon_deg];
end
latLons_apses_LoPO_2P2_1T__ns1_ZH = zeros(length(t_ev_ZH),2);
for kk = 1:length(t_ev_ZH)
    [lat_deg, lon_deg] = BCR2latlon(x_ev_ZH(kk,1:3)', 'secondary', prms_ZH.u);
    latLons_apses_LoPO_2P2_1T__ns1_ZH(kk,:) = [lat_deg, lon_deg];
end

% figure; hold all
% xlim([0 360])
% ylim([-90 90])
% PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')
% h = image(xlim, -ylim, bodies.europa.img);
% [lons_new] = convert_lon180_to_lon360(latLons_LoPO_2P2_1T__ns1(:,2));
% plot(lons_new,latLons_LoPO_2P2_1T__ns1(:,1),'k.')
% [lons_new_impact] = convert_lon180_to_lon360(latLons_apses_LoPO_2P2_1T__ns1((apses_alts < 5/rNorm),2));
% p_landing = plot(lons_new_impact, latLons_apses_LoPO_2P2_1T__ns1((apses_alts < 5/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.drkred, 'markerfacecolor', colors.red2);
% legend([p_landing], 'Tangent Impact', 'fontsize', 20, 'fontname', 'Times New Roman', 'location', 'northeast')
% title('LoPO\_2P2\_1T ... ns1','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)


figure; hold all
xlim([0 360])
ylim([-90 90])
PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')
h = image(xlim, -ylim, bodies.europa.img);
[lons_new_ZH] = convert_lon180_to_lon360(latLons_LoPO_2P2_1T__ns1_ZH(:,2));
plot(lons_new_ZH,latLons_LoPO_2P2_1T__ns1_ZH(:,1),'.', 'color', colors.red2)
[lons_new_impact_ZH] = convert_lon180_to_lon360(latLons_apses_LoPO_2P2_1T__ns1_ZH((apses_alts_ZH < 5/rNorm),2));
p_landing_ZH = plot(lons_new_impact_ZH, latLons_apses_LoPO_2P2_1T__ns1_ZH((apses_alts_ZH < 5/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.black, 'markerfacecolor', colors.drkgrey);
legend([p_landing_ZH], 'Tangent Impact', 'fontsize', 20, 'fontname', 'Times New Roman', 'location', 'northeast')
title('LoPO\_2P2\_1T ... ns1 ... ZH','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)

% figure; hold all
% xlim([0 360])
% ylim([-90 90])
% PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')
% h = image(xlim, -ylim, bodies.europa.img);
% plot(lons_new,latLons_LoPO_2P2_1T__ns1(:,1),'.','color',colors.blue2)
% plot(lons_new_ZH,latLons_LoPO_2P2_1T__ns1_ZH(:,1),'.', 'color', colors.red2)
% p_landing = plot(lons_new_impact, latLons_apses_LoPO_2P2_1T__ns1((apses_alts < 5/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.drkblue, 'markerfacecolor', colors.blue2);
% p_landing_ZH = plot(lons_new_impact_ZH, latLons_apses_LoPO_2P2_1T__ns1_ZH((apses_alts_ZH < 5/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.drkred, 'markerfacecolor', colors.red2);
% legend([p_landing, p_landing_ZH], 'Tangent Impact - CR3BP', 'Tangent Impact - CR3BP w ZH', 'fontsize', 20, 'fontname', 'Times New Roman', 'location', 'northeast')
% title('LoPO\_2P2\_1T ... ns1 ... Both','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)


figure; hold all
% plot3(X_LoPO_2P2_1T__ns1_ZH(:,1),X_LoPO_2P2_1T__ns1_ZH(:,2),X_LoPO_2P2_1T__ns1_ZH(:,3), 'linewidth', 2, 'color', colors.red2)
plot3(X_LoPO_2P2_1T__ns1_ZH(:,1),X_LoPO_2P2_1T__ns1_ZH(:,2),X_LoPO_2P2_1T__ns1_ZH(:,3), 'linewidth', 2, 'color', colors.mag)
plotTrajShadows(X_LoPO_2P2_1T__ns1_ZH, 2, colors.grey, 'x', 1.017, 'y', 2.1e-2, 'z', -22e-3, 'bodyshadow', [1-prms_ZH.u, prms_ZH.R2])
plotSecondary(secondary)
PlotBoi3_CR3Bn(26)
axis equal
plot3(x_ev_ZH((apses_alts_ZH < 5/rNorm),1), x_ev_ZH((apses_alts_ZH < 5/rNorm),2), x_ev_ZH((apses_alts_ZH < 5/rNorm),3), '.','markersize',20, 'color', colors.black)
set(gca,'yticklabel',num2str(get(gca,'ytick')','%1.3f'))
set(gca,'zticklabel',num2str(get(gca,'ztick')','%1.3f'))



end % run_LoPO_2P2_1T__ns1
%% ------------------------------------------------
%%% LoPO_2P4__ns1
% -------------------------------------------------
if run_LoPO_2P4__ns1
fprintf('---------------------------------')
fprintf('LoPO_2P4__ns1 (planar)\n')

[T_LoPO_2P4__ns1, X_LoPO_2P4__ns1, t_ev, x_ev, in_ev] = ode113(@Int_CR3BnSTM, [0, PO_LoPO_2P4__ns1(end)], [PO_LoPO_2P4__ns1(1:6); stm0_colVec], options_apsis, prms);
apses_alts = rowNorm(x_ev(:,1:3) - [1-prms.u,0,0]) - prms.R2;
% figure; hold all
% plot3(X_LoPO_2P4__ns1(:,1),X_LoPO_2P4__ns1(:,2),X_LoPO_2P4__ns1(:,3), 'linewidth', 2, 'color', colors.black)
% plotSecondary(secondary)
% PlotBoi3_CR3Bn(26)
% axis equal
% title('LoPO\_2P4 ... ns1','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)
% 
% plot3(x_ev((apses_alts < 5/rNorm),1), x_ev((apses_alts < 5/rNorm),2), x_ev((apses_alts < 5/rNorm),3), '.','markersize',20, 'color', colors.red2)
% set(gca,'yticklabel',num2str(get(gca,'ytick')','%1.3f'))
% set(gca,'zticklabel',num2str(get(gca,'ztick')','%1.3f'))


[T_LoPO_2P4__ns1_ZH, X_LoPO_2P4__ns1_ZH, t_ev_ZH, x_ev_ZH, in_ev_ZH] = ode113(@Int_CR3BnSTM_J2pJ4pJ6pJ2s, [0, PO_LoPO_2P4__ns1_ZH(end)], [PO_LoPO_2P4__ns1_ZH(1:6); stm0_colVec], options_apsis, prms_ZH);
fprintf('Prop Error = %1.2e\n', norm(X_LoPO_2P4__ns1_ZH(end,1:6) - X_LoPO_2P4__ns1_ZH(1,1:6)) / norm(X_LoPO_2P4__ns1_ZH(1,1:6)))
figure; hold all
% plot3(X_LoPO_2P4__ns1_ZH(:,1),X_LoPO_2P4__ns1_ZH(:,2),X_LoPO_2P4__ns1_ZH(:,3), 'linewidth', 2, 'color', colors.red2)
plot3(X_LoPO_2P4__ns1_ZH(:,1),X_LoPO_2P4__ns1_ZH(:,2),X_LoPO_2P4__ns1_ZH(:,3), 'linewidth', 2, 'color', colors.mag)
plotSecondary(secondary)
PlotBoi3_CR3Bn(26)
axis equal
title('LoPO\_2P4 ... ns1 ... ZH','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)

apses_alts_ZH = rowNorm(x_ev_ZH(:,1:3) - [1-prms_ZH.u,0,0]) - prms_ZH.R2;
plot3(x_ev_ZH((apses_alts_ZH < 5/rNorm),1), x_ev_ZH((apses_alts_ZH < 5/rNorm),2), x_ev_ZH((apses_alts_ZH < 5/rNorm),3), '.','markersize',20, 'color', colors.black)
set(gca,'yticklabel',num2str(get(gca,'ytick')','%1.3f'))
set(gca,'zticklabel',num2str(get(gca,'ztick')','%1.3f'))

view(0,90)

% figure; hold all
% p_3B = plot3(X_LoPO_2P4__ns1(:,1),X_LoPO_2P4__ns1(:,2),X_LoPO_2P4__ns1(:,3), 'linewidth', 2, 'color', colors.blue2);
% p_ZH = plot3(X_LoPO_2P4__ns1_ZH(:,1),X_LoPO_2P4__ns1_ZH(:,2),X_LoPO_2P4__ns1_ZH(:,3), 'linewidth', 2, 'color', colors.red2);
% plotSecondary(secondary)
% PlotBoi3_CR3Bn(26)
% axis equal
% legend([p_3B, p_ZH], 'CR3BP', 'CR3BP w ZH')
% set(gca,'yticklabel',num2str(get(gca,'ytick')','%1.3f'))
% set(gca,'zticklabel',num2str(get(gca,'ztick')','%1.3f'))


stm_tf_t0                           = reshape(X_LoPO_2P4__ns1(end,7:42),6,6);
monodromy                           = stm_tf_t0;
[eigenVectors_new, eigenValues_new] = eig(monodromy);
[S1, S2]                            = getStabilityIndices(diag(eigenValues_new));

stm_tf_t0_ZH                              = reshape(X_LoPO_2P4__ns1_ZH(end,7:42),6,6);
monodromy_ZH                              = stm_tf_t0_ZH;
[eigenVectors_new_ZH, eigenValues_new_ZH] = eig(monodromy_ZH);
[S1_ZH, S2_ZH]                            = getStabilityIndices(diag(eigenValues_new_ZH));

fprintf('[S1, S2]       = [%1.3f, %1.3f]\n', S1, S2)
fprintf('[S1_ZH, S2_ZH] = [%1.3f, %1.3f]\n', S1_ZH, S2_ZH)


latLons_LoPO_2P4__ns1 = zeros(length(T_LoPO_2P4__ns1),2);
for kk = 1:length(T_LoPO_2P4__ns1)
    [lat_deg, lon_deg] = BCR2latlon(X_LoPO_2P4__ns1(kk,1:3)', 'secondary', prms.u);
    latLons_LoPO_2P4__ns1(kk,:) = [lat_deg, lon_deg];
end
latLons_apses_LoPO_2P4__ns1 = zeros(length(t_ev),2);
for kk = 1:length(t_ev)
    [lat_deg, lon_deg] = BCR2latlon(x_ev(kk,1:3)', 'secondary', prms.u);
    latLons_apses_LoPO_2P4__ns1(kk,:) = [lat_deg, lon_deg];
end


latLons_LoPO_2P4__ns1_ZH = zeros(length(T_LoPO_2P4__ns1_ZH),2);
for kk = 1:length(T_LoPO_2P4__ns1_ZH)
    [lat_deg, lon_deg] = BCR2latlon(X_LoPO_2P4__ns1_ZH(kk,1:3)', 'secondary', prms_ZH.u);
    latLons_LoPO_2P4__ns1_ZH(kk,:) = [lat_deg, lon_deg];
end
latLons_apses_LoPO_2P4__ns1_ZH = zeros(length(t_ev_ZH),2);
for kk = 1:length(t_ev_ZH)
    [lat_deg, lon_deg] = BCR2latlon(x_ev_ZH(kk,1:3)', 'secondary', prms_ZH.u);
    latLons_apses_LoPO_2P4__ns1_ZH(kk,:) = [lat_deg, lon_deg];
end

% figure; hold all
% xlim([0 360])
% ylim([-90 90])
% PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')
% h = image(xlim, -ylim, bodies.europa.img);
% [lons_new] = convert_lon180_to_lon360(latLons_LoPO_2P4__ns1(:,2));
% plot(lons_new,latLons_LoPO_2P4__ns1(:,1),'k.')
% [lons_new_impact] = convert_lon180_to_lon360(latLons_apses_LoPO_2P4__ns1((apses_alts < 5/rNorm),2));
% p_landing = plot(lons_new_impact, latLons_apses_LoPO_2P4__ns1((apses_alts < 5/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.drkred, 'markerfacecolor', colors.red2);
% legend([p_landing], 'Tangent Impact', 'fontsize', 20, 'fontname', 'Times New Roman', 'location', 'northeast')
% title('LoPO\_2P4 ... ns1','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)


figure; hold all
xlim([0 360])
ylim([-90 90])
PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')
h = image(xlim, -ylim, bodies.europa.img);
[lons_new_ZH] = convert_lon180_to_lon360(latLons_LoPO_2P4__ns1_ZH(:,2));
plot(lons_new_ZH,latLons_LoPO_2P4__ns1_ZH(:,1),'.', 'color', colors.red2)
[lons_new_impact_ZH] = convert_lon180_to_lon360(latLons_apses_LoPO_2P4__ns1_ZH((apses_alts_ZH < 5/rNorm),2));
p_landing_ZH = plot(lons_new_impact_ZH, latLons_apses_LoPO_2P4__ns1_ZH((apses_alts_ZH < 5/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.black, 'markerfacecolor', colors.drkgrey);
legend([p_landing_ZH], 'Tangent Impact', 'fontsize', 20, 'fontname', 'Times New Roman', 'location', 'northeast')
title('LoPO\_2P4 ... ns1 ... ZH','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)

% figure; hold all
% xlim([0 360])
% ylim([-90 90])
% PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')
% h = image(xlim, -ylim, bodies.europa.img);
% plot(lons_new,latLons_LoPO_2P4__ns1(:,1),'.','color',colors.blue2)
% plot(lons_new_ZH,latLons_LoPO_2P4__ns1_ZH(:,1),'.', 'color', colors.red2)
% p_landing = plot(lons_new_impact, latLons_apses_LoPO_2P4__ns1((apses_alts < 5/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.drkblue, 'markerfacecolor', colors.blue2);
% p_landing_ZH = plot(lons_new_impact_ZH, latLons_apses_LoPO_2P4__ns1_ZH((apses_alts_ZH < 5/rNorm),1), 'o','markersize',16, 'linewidth', 2, 'markeredgecolor', colors.drkred, 'markerfacecolor', colors.red2);
% legend([p_landing, p_landing_ZH], 'Tangent Impact - CR3BP', 'Tangent Impact - CR3BP w ZH', 'fontsize', 20, 'fontname', 'Times New Roman', 'location', 'northeast')
% title('LoPO\_2P4 ... ns1 ... Both','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)

end % run_LoPO_2P4__ns1





