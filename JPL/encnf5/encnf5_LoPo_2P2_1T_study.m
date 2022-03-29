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

%%% Integration options
tol = 1e-13;
options = odeset('RelTol',tol,'AbsTol',tol);

% -------------------------------------------------
%%% Load data
% -------------------------------------------------
%%% Load the data file
PO_data = dlmread(PO_file,',',1,0);

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
    POs{PO_i}.T         = t_PO;
    POs{PO_i}.traj      = X_PO(:,1:6);
    POs{PO_i}.t_ca      = t_aps(closeApproachIndex);
    POs{PO_i}.x_ca      = X_aps(closeApproachIndex, 1:6);
    POs{PO_i}.eigs      = diag(eigenValues_new);
    POs{PO_i}.S1S2      = [S1, S2];
    POs{PO_i}.JC        = getJacobiConstant_ZH(X_PO(1,1:6), prms);
    POs{PO_i}.alt_ca_km = closeApproach_alt*rNorm;
    POs{PO_i}.latlon_ca = closeApproach_latlon;

end

%%% Converting data to struct array and parsing
POs_structArray = [POs{:}];

data_alt_ca_km  = [POs_structArray(:).alt_ca_km]';
data_JC         = [POs_structArray(:).JC]';
data_S1S2       = reshape([POs_structArray(:).S1S2], 2, n_POs)';
data_eigs       = reshape([POs_structArray(:).eigs], 6, n_POs)';
data_latlon     = reshape([POs_structArray(:).latlon_ca], 2, n_POs)';

% ========================================================================
%%% Close approach altitude     vs      Jacobi Constant
% ========================================================================
figure; hold all

plot3(data_alt_ca_km, data_JC, 1:n_POs, 'k.', 'markersize', 20)

H = gca;
H.FontSize = 14;
PlotBoi2('Close Approach Altitude (km)', 'Jacobi Constant', 26, 'LaTex')


% ========================================================================
%%% Close approach altitude     vs      Stability index
% ========================================================================
figure; hold all

p_s1 = plot3(data_alt_ca_km, data_S1S2(:,1), 1:n_POs, 'r.', 'markersize', 20);
p_s2 = plot3(data_alt_ca_km, data_S1S2(:,2), 1:n_POs, 'b.', 'markersize', 20);

H = gca;
H.FontSize = 14;
PlotBoi2('Close Approach Altitude (km)', 'Stability Indices', 26, 'LaTex')

legend([p_s1, p_s2], 'S1', 'S2', 'FontSize', 16)


% ========================================================================
%%% Close approach altitude     vs      Latitude
% ========================================================================
figure; hold all

plot3(data_alt_ca_km, data_latlon(:,1), 1:n_POs, 'k.', 'markersize', 20)

H = gca;
H.FontSize = 14;
PlotBoi2('Close Approach Altitude (km)', 'Close Approach Latitude ($^\circ$)', 26, 'LaTex')


% ========================================================================
%%% Orbits
% ========================================================================

figure; hold all







plot_PO_indices = getIndices_spacedByTpJcArclength(PO_data(:,c_Tp), PO_data(:,c_JC), 25);
% plot_PO_indices = 50;

%%% optional color spectrum
% color_spectrum = colors.blue2;
color_spectrum = colorScale([colors.blue2; colors.mag], length(plot_PO_indices));


%%% Loop and plot
figure; hold all
color_index = 0;
for PO_i = plot_PO_indices'
    color_index = color_index + 1;
    plot3(POs{PO_i}.traj(:,1), POs{PO_i}.traj(:,2), POs{PO_i}.traj(:,3), 'linewidth', 1.5, 'color', color_spectrum(color_index,:)) % 0.1961    0.3922    0.7843
end
ax = gca;
ax.FontSize = 14;
PlotBoi3_CR3Bn(28)
axis equal

plotSecondary(secondary)















