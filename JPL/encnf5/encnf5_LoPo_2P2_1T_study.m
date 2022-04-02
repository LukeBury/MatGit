% ========================================================================
%%% Description
% ========================================================================
% 

% Created: 3/29/22
% Author : Luke Bury, luke.bury@jpl.nasa.gov
% ========================================================================
%%% Initialization
% ========================================================================
clear
clc
close all
% mbinPath = '~/CU_Google_Drive/Documents/MatGit/mbin';
mbinPath = '~/Documents/MatGit/mbin';
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
%%% Choose data
% -------------------------------------------------
PO_file = 'Saturn_Enceladus.CR3BP.LoPO_2P2_1T_newMR.txt'; 

% -------------------------------------------------
%%% Set up parameters
% -------------------------------------------------
%%% Set primary and secondary bodies
[primary, secondary] = assignPrimaryAndSecondary_CR3BP(PO_file, bodies);

%%% Normalizing constants
[rNorm, tNorm, vNorm] = cr3bp_norms(primary, secondary, bodies.constants.G);

%%% prms for integration
secondary.MR = 1.898884589251784e-07;
prms.u  = secondary.MR;
prms.R2 = secondary.R_n;
prms.n    = 1;

%Equilibrium Points
rLPs_n = EquilibriumPoints(prms.u, prms.n);

% -------------------------------------------------
%%% Load data
% -------------------------------------------------
%%% Load the data file
PO_data = dlmread(PO_file,',',1,0);

warning('Ignoring off data past 200 km ca-altitude')
PO_data = PO_data(99:end,:);

%%% Grab header line
fid = fopen(PO_file, 'rt');  %the 't' is important!
header = fgetl(fid);
fclose(fid);

%%% Number of ICs
n_POs = size(PO_data,1);

%%% Indices for plotting purposes
PO_indices = linspace(1, n_POs, n_POs);
% --------------------------
% Create column specifiers
% --------------------------
PO_header_2020 = 'x0,y0,z0,xd0,yd0,zd0,Tp,JC,stabilityIndex1,stabilityIndex2,alpha,beta,impactFlag,error';

if contains(header, PO_header_2020)
    c_x0 = 1;   c_y0 = 2;   c_z0 = 3;
    c_xd0 = 4;  c_yd0 = 5;  c_zd0 = 6;
    c_Tp = 7;   c_JC = 8;   c_S1 = 9;   c_S2 = 10;
    c_alpha = 11;   c_beta = 12;    c_impactFlag = 13;
    c_error = 14;
end


% ========================================================================
%%% Integrate orbits
% ========================================================================
%%% Integration options
tol = 1e-13;
options = odeset('RelTol',tol,'AbsTol',tol);
options_apsis = odeset('Events', @event_Apsis_CR3BP, 'RelTol',tol,'AbsTol',tol);

%%% Preallocate cell array for PO outputs
POs = cell(n_POs,1);

%%% Initialize STM as column vector to be added with state
stm0_colVec = reshape(eye(6),36,1);

% n_POs = 5
parfor PO_i = 1:n_POs
    %%% Integrate PO and track all apses
    [t_PO, X_PO, t_aps, X_aps, in_aps] = ode113(@Int_CR3BnSTM, [0, PO_data(PO_i, c_Tp)], [PO_data(PO_i, c_x0:c_zd0)'; stm0_colVec], options_apsis, prms);
    
    %%% Calculate stability info
    stm_tf_t0                           = reshape(X_PO(end,7:42),6,6);
    monodromy                           = stm_tf_t0;
    [eigenVectors_new, eigenValues_new] = eig(monodromy);
    [S1, S2]                            = getStabilityIndices(diag(eigenValues_new));
    
    %%% Determine closest approach and generate info
    aps_altitudes = rowNorm(X_aps(:,1:3) - [1-prms.u, 0, 0]) - prms.R2;
    closeApproachIndex = find(aps_altitudes == min(aps_altitudes));
    closeApproach_alt = aps_altitudes(closeApproachIndex);
    [lat, lon] = BCR2latlon(X_aps(closeApproachIndex, 1:3), 'secondary', prms.u);
    closeApproach_latlon = [lat, lon];

    %%% Store PO data in cell arrays
    POs{PO_i}.T            = t_PO;
    POs{PO_i}.traj         = X_PO(:,1:6);
    POs{PO_i}.t_aps        = t_aps;
    POs{PO_i}.X_aps        = X_aps(:,1:6);
    POs{PO_i}.aps_alts_km  = aps_altitudes.*rNorm;
    POs{PO_i}.t_ca         = t_aps(closeApproachIndex);
    POs{PO_i}.x_ca         = X_aps(closeApproachIndex, 1:6);
    POs{PO_i}.eigs         = diag(eigenValues_new);
    POs{PO_i}.S1S2         = [S1, S2];
    POs{PO_i}.JC           = getJacobiConstant_ZH(X_PO(1,1:6), prms);
    POs{PO_i}.alt_ca_km    = closeApproach_alt.*rNorm;
    POs{PO_i}.latlon_ca    = closeApproach_latlon;

end

%%% Converting data to struct array and parsing
POs_structArray = [POs{:}];

data_alt_ca_km  = [POs_structArray(:).alt_ca_km]';
data_JC         = [POs_structArray(:).JC]';
data_S1S2       = reshape([POs_structArray(:).S1S2], 2, n_POs)';
data_eigs       = reshape([POs_structArray(:).eigs], 6, n_POs)';
data_latlon     = reshape([POs_structArray(:).latlon_ca], 2, n_POs)';



% ========================================================================
%%% Plotting Energy vs Tp
% ========================================================================
%%% Calculating time period in days
data_Tp_days = PO_data(:,c_Tp).*tNorm/86400;

%%% Plot Jacobi constant and Tp
figure; hold all
plot3(data_Tp_days, data_JC, 1:n_POs,'o', 'markeredgecolor',colors.black,'markerfacecolor',colors.ltgrey)
PlotBoi2('$T_P$ (days)','Jacobi Constant',26,'LaTex')

% ========================================================================
%%% Closest approach altitude     vs      Jacobi Constant
% ========================================================================
figure; hold all

plot3(data_alt_ca_km, data_JC, 1:n_POs, 'o', 'markeredgecolor',colors.black,'markerfacecolor',colors.ltgrey)

H = gca;
H.FontSize = 14;
PlotBoi2('Closest Approach Altitude (km)', 'Jacobi Constant', 26, 'LaTex')

% ========================================================================
%%% Closest approach altitude     vs      Time Period
% ========================================================================
figure; hold all

plot3(data_alt_ca_km, data_Tp_days, 1:n_POs, 'o', 'markeredgecolor',colors.black,'markerfacecolor',colors.ltgrey)

H = gca;
H.FontSize = 14;
PlotBoi2('Closest Approach Altitude (km)', '$T_P$ (days)', 26, 'LaTex')



% ========================================================================
%%% Closest approach altitude     vs      Stability index
% ========================================================================
figure('position', [476 506 947 360]);

subplot(1,2,1); hold all
p_s1 = plot3(data_alt_ca_km, data_S1S2(:,1), 1:n_POs, 'o','markeredgecolor',colors.blue,'markerfacecolor',colors.ltblue);
p_s2 = plot3(data_alt_ca_km, data_S1S2(:,2), 1:n_POs, 'o','markeredgecolor',colors.red,'markerfacecolor',colors.ltred);
plot3([min(data_alt_ca_km) max(data_alt_ca_km)],[2 2], [1 1].*1000,'k','linewidth',1)
plot3([min(data_alt_ca_km) max(data_alt_ca_km)],[-2 -2], [1 1].*1000,'k','linewidth',1)

H = gca;
H.FontSize = 14;
PlotBoi2('Closest Approach Altitude (km)', 'Stability Indices', 26, 'LaTex')

legend([p_s1, p_s2], 'S1', 'S2', 'FontSize', 16, 'location', 'best')
xlim([min(data_alt_ca_km) max(data_alt_ca_km)])

subplot(1,2,2); hold all
p_s1 = plot3(data_alt_ca_km, data_S1S2(:,1), 1:n_POs, 'o','markeredgecolor',colors.blue,'markerfacecolor',colors.ltblue);
p_s2 = plot3(data_alt_ca_km, data_S1S2(:,2), 1:n_POs, 'o','markeredgecolor',colors.red,'markerfacecolor',colors.ltred);
plot3([min(data_alt_ca_km) max(data_alt_ca_km)],[2 2], [1 1].*1000,'k','linewidth',1)
plot3([min(data_alt_ca_km) max(data_alt_ca_km)],[-2 -2], [1 1].*1000,'k','linewidth',1)

H = gca;
H.FontSize = 14;
PlotBoi2('Closest Approach Altitude (km)', 'Stability Indices', 26, 'LaTex')
ylim([-1 1].*3)
xlim([min(data_alt_ca_km) max(data_alt_ca_km)])


% ========================================================================
%%% Closest approach altitude     vs      Latitude
% ========================================================================
figure; hold all

plot3(data_alt_ca_km, data_latlon(:,1), 1:n_POs, 'o', 'markeredgecolor',colors.black,'markerfacecolor',colors.ltgrey)

H = gca;
H.FontSize = 14;
PlotBoi2('Closest Approach Altitude (km)', 'Closest Approach Latitude ($^\circ$)', 26, 'LaTex')


% ========================================================================
%%% Looking at apses
% ========================================================================

figure; hold all
for PO_i = 1:n_POs
%     plot3(data_JC(PO_i), POs{PO_i}.aps_alts_km(1), PO_i, 'o', 'markeredgecolor',colors.black,'markerfacecolor',colors.ltgrey)
%     plot3(data_JC(PO_i), POs{PO_i}.aps_alts_km(2), PO_i, 'o', 'markeredgecolor',colors.black,'markerfacecolor',colors.ltgrey)
%     plot3(data_JC(PO_i), POs{PO_i}.aps_alts_km(3), PO_i, 'o', 'markeredgecolor',colors.black,'markerfacecolor',colors.ltgrey)
%     plot3(data_JC(PO_i), POs{PO_i}.aps_alts_km(4), PO_i, 'o', 'markeredgecolor',colors.black,'markerfacecolor',colors.ltgrey)
%     plot3(data_JC(PO_i), POs{PO_i}.aps_alts_km(5), PO_i, 'o', 'markeredgecolor',colors.black,'markerfacecolor',colors.ltgrey)
%     plot3(data_JC(PO_i), POs{PO_i}.aps_alts_km(6), PO_i, 'o', 'markeredgecolor',colors.black,'markerfacecolor',colors.ltgrey)
% 
%     if length(POs{PO_i}.aps_alts_km) >= 7
%         plot3(data_JC(PO_i), POs{PO_i}.aps_alts_km(7), PO_i, 'o', 'markeredgecolor',colors.black,'markerfacecolor',colors.ltgrey)
%         if length(POs{PO_i}.aps_alts_km) == 8
%             plot3(data_JC(PO_i), POs{PO_i}.aps_alts_km(8), PO_i, 'o', 'markeredgecolor',colors.black,'markerfacecolor',colors.ltgrey)
%         end
%     end

    plot3(POs{PO_i}.aps_alts_km(1), data_Tp_days(PO_i), PO_i, 'o', 'markeredgecolor',colors.black,'markerfacecolor',colors.ltgrey)
    plot3(POs{PO_i}.aps_alts_km(2), data_Tp_days(PO_i), PO_i, 'o', 'markeredgecolor',colors.black,'markerfacecolor',colors.ltgrey)
    plot3(POs{PO_i}.aps_alts_km(3), data_Tp_days(PO_i), PO_i, 'o', 'markeredgecolor',colors.black,'markerfacecolor',colors.ltgrey)
    plot3(POs{PO_i}.aps_alts_km(4), data_Tp_days(PO_i), PO_i, 'o', 'markeredgecolor',colors.black,'markerfacecolor',colors.ltgrey)
    plot3(POs{PO_i}.aps_alts_km(5), data_Tp_days(PO_i), PO_i, 'o', 'markeredgecolor',colors.black,'markerfacecolor',colors.ltgrey)
    plot3(POs{PO_i}.aps_alts_km(6), data_Tp_days(PO_i), PO_i, 'o', 'markeredgecolor',colors.black,'markerfacecolor',colors.ltgrey)

    if length(POs{PO_i}.aps_alts_km) >= 7
        plot3(POs{PO_i}.aps_alts_km(7), data_Tp_days(PO_i), PO_i, 'o', 'markeredgecolor',colors.black,'markerfacecolor',colors.ltgrey)
        if length(POs{PO_i}.aps_alts_km) == 8
            plot3(POs{PO_i}.aps_alts_km(8), data_Tp_days(PO_i), PO_i, 'o', 'markeredgecolor',colors.black,'markerfacecolor',colors.ltgrey)
        end
    end
end

PlotBoi2('Apsis Altitudes (km)', 'Time Period (days)', 26, 'LaTex')


% ========================================================================
%%% Orbits
% ========================================================================


n_plot_POs = 13;
plot_PO_indices = getIndices_spacedByTpJcArclength(PO_data(:,c_Tp), PO_data(:,c_JC), n_plot_POs);
% plot_PO_indices = 78;

%%% optional color spectrum
% color_spectrum = colors.blue2;
color_spectrum = colorScale([colors.blue2; colors.mag], length(plot_PO_indices));


%%% Loop and plot
figure; hold all
color_index = 0;
for PO_i = plot_PO_indices'
    color_index = color_index + 1;
    plot3(POs{PO_i}.traj(:,1), POs{PO_i}.traj(:,2), POs{PO_i}.traj(:,3), 'linewidth', 2, 'color', color_spectrum(color_index,:)) % 0.1961    0.3922    0.7843
end
ax = gca;
ax.FontSize = 14;
PlotBoi3_CR3Bn(28)
axis equal

plotSecondary(secondary)


% plot3(POs{n_POs}.X_aps(:,1), POs{n_POs}.X_aps(:,2), POs{n_POs}.X_aps(:,3), 'k.', 'markersize', 20)
% plot3(POs{1}.X_aps(:,1), POs{1}.X_aps(:,2), POs{1}.X_aps(:,3), 'k.', 'markersize', 13)
% 
% 
% plotTrajShadows(POs{PO_i}.traj, 2, colors.grey, 'x', 1.004, 'y', 10e-3, 'z', -4e-3, 'bodyshadow', [1-prms.u, prms.R2])




%%% For plotting a single orbit
if 1+1==1

%     PO_i = 206; % 20 km, 206
% % %     PO_i = 166; % 50 km, 166
%     PO_i = 102; % 100 km, 102
% % %     PO_i = 49; % 150 km, 49
    PO_i = 2; % 200 km, 2

    %%% Loop and plot
    figure; hold all
    plot3(POs{PO_i}.traj(:,1), POs{PO_i}.traj(:,2), POs{PO_i}.traj(:,3), 'linewidth', 2, 'color', colors.blue2) % 0.1961    0.3922    0.7843
    ax = gca;
    ax.FontSize = 14;
    PlotBoi3_CR3Bn(28)
    axis equal
    plotSecondary(secondary)
%     plot3(POs{PO_i}.X_aps([1,3,5],1), POs{PO_i}.X_aps([1,3,5],2), POs{PO_i}.X_aps([1,3,5],3), 'o', 'markeredgecolor',colors.black,'markerfacecolor',colors.ltgrey)
    plotTrajShadows(POs{PO_i}.traj, 2, colors.grey, 'x', 1.004, 'y', 10e-3, 'z', -4e-3, 'bodyshadow', [1-prms.u, prms.R2])

    axis normal
    axis equal
    xlim([0.997 1.004])
    view(-28,16)


    latLons = zeros(length(POs{PO_i}.T),2);
    for kk = 1:length(POs{PO_i}.T)
        [lat_deg, lon_deg] = BCR2latlon(POs{PO_i}.traj(kk,1:3)', 'secondary', prms.u);
        latLons(kk,:) = [lat_deg, lon_deg];
    end
    figure; hold all
    xlim([0 360])
    ylim([-90 90])
    PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')
    h = image(xlim, -ylim, bodies.enceladus.img);
    [lons_new] = convert_lon180_to_lon360(latLons(:,2));
    plot(lons_new,latLons(:,1),'.', 'color', colors.blue2)
%     legend([p_landing_ZH], 'Tangent Impact', 'fontsize', 20, 'fontname', 'Times New Roman', 'location', 'southeast')
%     title('SHalo ... ns1 ... ZH','Interpreter', 'LaTex', 'Fontname', 'Times New Roman', 'FontSize', 16)




end




%%% Inertial
if 1+1==1
%     PO_i = 206; % 20 km, 206
% % %     PO_i = 166; % 50 km, 166
    PO_i = 102; % 100 km, 102
% % %     PO_i = 49; % 150 km, 49
%     PO_i = 2; % 200 km, 2

    figure; hold all
    plot3(POs{PO_i}.traj(:,1), POs{PO_i}.traj(:,2), POs{PO_i}.traj(:,3), 'linewidth', 2, 'color', colors.black)
    axis equal
    PlotBoi3_CR3Bn(28)
    
    PO_x2 = [POs{PO_i}.traj; POs{PO_i}.traj];
    T_x2 = [POs{PO_i}.T; POs{PO_i}.T + POs{PO_i}.T(end)];

    PO_x3 = [POs{PO_i}.traj; POs{PO_i}.traj; POs{PO_i}.traj];
    T_x3 = [POs{PO_i}.T; POs{PO_i}.T + POs{PO_i}.T(end); POs{PO_i}.T + POs{PO_i}.T(end)*2];
    
    PO_x4 = [POs{PO_i}.traj; POs{PO_i}.traj; POs{PO_i}.traj; POs{PO_i}.traj];
    T_x4 = [POs{PO_i}.T; POs{PO_i}.T + POs{PO_i}.T(end); POs{PO_i}.T + POs{PO_i}.T(end)*2; POs{PO_i}.T + POs{PO_i}.T(end)*3];

    [X_SCI] = X_BaCR2SCI(POs{PO_i}.traj, POs{PO_i}.T, prms);
    [X_SCI_x2] = X_BaCR2SCI(PO_x2, T_x2, prms);
    [X_SCI_x3] = X_BaCR2SCI(PO_x3, T_x3, prms);
    [X_SCI_x4] = X_BaCR2SCI(PO_x4, T_x4, prms);


    figure; hold all
%     plot3(X_SCI(:,1), X_SCI(:,2), X_SCI(:,3), 'linewidth', 2, 'color', colors.black)
%     plot3(X_SCI_x2(:,1), X_SCI_x2(:,2), X_SCI_x2(:,3), 'linewidth', 2, 'color', colors.black)
%     plot3(X_SCI_x3(:,1), X_SCI_x3(:,2), X_SCI_x3(:,3), 'linewidth', 2, 'color', colors.black)
    plot3(X_SCI_x4(:,1), X_SCI_x4(:,2), X_SCI_x4(:,3), 'linewidth', 2, 'color', colors.black)
    axis equal
    PlotBoi3_CR3Bn(28)
    
    plotBody3(prms.R2, [0, 0, 0], colors.ltblue, 0.5)

end




