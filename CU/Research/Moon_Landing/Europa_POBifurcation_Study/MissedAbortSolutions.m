% ========================================================================
%%% Description
% ========================================================================
% For plotting families of POs and describing the bifurcations that occur
% within them. For plotting, designed to integrate the trajectories in
% parallel, but the 'parfor' can be switched to 'for' if small, fast
% batches are desired

% Created: 12/16/20
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


% ========================================================================
%%% Setup
% ========================================================================
% -------------------------------------------------
%%% Personal options
% -------------------------------------------------


% -------------------------------------------------
%%% Choose data
% -------------------------------------------------
family = 'Jupiter_Europa.CR3BP_J2pJ4pJ6pJ2s.Hg2_2T.txt'; 
% family = 'Jupiter_Europa.CR3BP_J2pJ4pJ6pJ2s.Hg2_3T_ns1.txt'; 
% family = 'Jupiter_Europa.CR3BP_J2pJ4pJ6pJ2s.Se7_6P2_ns1.txt'; 



% --------------------------
% Actually load family
% --------------------------
%%% Path from mbin to data
dataPathFromMBin = '/Data/InitialConditions/PO_Families/';

%%% PO data file
PO_datafile = [mbinPath, dataPathFromMBin, family];
% -------------------------------------------------
%%% Set up parameters
% -------------------------------------------------
%%% Set primary and secondary bodies
[primary, secondary] = assignPrimaryAndSecondary_CR3BP(family, bodies);

%%% Normalizing constants
[rNorm, tNorm, vNorm] = cr3bp_norms(primary, secondary, bodies.constants.G);
    

%%% prms for integration
prms_ZH.u   = secondary.MR;
prms_ZH.R2  = secondary.R_n;
prms_ZH.J2p = primary.J2;
prms_ZH.J4p = primary.J4;
prms_ZH.J6p = primary.J6;
prms_ZH.J2s = secondary.J2;

prms_ZH.R1 = primary.R / rNorm;

%%% Determine the mean motion via the ephemeris method
tN_ephemeris = sqrt((secondary.a^3) / (bodies.constants.G*(primary.mass+secondary.mass)));
prms_ZH.n = secondary.meanMot*tN_ephemeris;

%%% Collinear equillibrium points
rLPs_n_ZH = collinearEquilibriumPoints_ZH(prms_ZH);

%%% Integration options
tol = 1e-13;
options = odeset('RelTol',tol,'AbsTol',tol);
options_apsis = odeset('Events', @event_Apsis_CR3BP, 'RelTol',tol,'AbsTol',tol);

%%% Initialize STM as column vector to be added with state
stm0_colVec = reshape(eye(6),36,1);
% -------------------------------------------------
%%% Load data
% -------------------------------------------------
%%% Load the data file
PO_data = dlmread(PO_datafile,',',1,0);

%%% Grab header line
fid = fopen(PO_datafile, 'rt');  %the 't' is important!
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
%%% Look at apses
% ========================================================================
plot_PO_indices = 1:n_POs;



if isequal(family, 'Jupiter_Europa.CR3BP_J2pJ4pJ6pJ2s.Hg2_2T.txt')
    %%% Preallocate
    apse_1 = cell(length(plot_PO_indices),1);
    apse_2 = cell(length(plot_PO_indices),1);
    
    lat_1 = cell(length(plot_PO_indices),1);

    parfor PO_i = 1:length(plot_PO_indices)
        index = plot_PO_indices(PO_i);

        [T_PO, X_PO, t_ev, X_ev, in_ev] = ode113(@Int_CR3BnSTM_J2pJ4pJ6pJ2s, [0, PO_data(index,c_Tp)], [PO_data(index,c_x0:c_zd0)'; stm0_colVec], options_apsis, prms_ZH);

        X_ev = X_ev((0.997 < X_ev(:,1)) & ( X_ev(:,1) < 0.999), 1:6);
        
        apse_1{PO_i} = [X_ev(1,1:6), PO_data(index,c_Tp)];
        apse_2{PO_i} = [X_ev(2,1:6), PO_data(index,c_Tp)];
        
        [lat_deg, lon_deg] = BCR2latlon(apse_1{PO_i}(1:3)', 'secondary', prms_ZH.u);
        lat_1{PO_i} = lat_deg;
    end
    
    
    apse_1 = cell2mat(apse_1);
    apse_2 = cell2mat(apse_2);
    
    lat_1 = cell2mat(lat_1);

    alt_apse1 = (rowNorm(apse_1(:,1:3) - [1-prms_ZH.u, 0, 0]) - secondary.R_n).*rNorm;
    alt_apse2 = (rowNorm(apse_2(:,1:3) - [1-prms_ZH.u, 0, 0]) - secondary.R_n).*rNorm;
    
    figure; plot(apse_1(:,7), alt_apse1, 'k', 'linewidth', 1.5)
    PlotBoi2('$T_P$', 'Altitude (km)', 32, 'LaTex')
	
    figure; plot(alt_apse1, apse_1(:,3), 'k', 'linewidth', 1.5)
    PlotBoi2('Altitude (km)', '$z_n$', 32, 'LaTex')
    xlim([0, 300])
%     figure; plot(alt_apse2)

    figure; plot(alt_apse1, lat_1, 'k', 'linewidth', 1.5)
    PlotBoi2('Altitude (km)', 'Latitude ($^\circ$)', 32, 'LaTex')
    xlim([0, 300])

end


if isequal(family, 'Jupiter_Europa.CR3BP_J2pJ4pJ6pJ2s.Hg2_3T_ns1.txt')
    %%% Preallocate
    apse_1 = cell(length(plot_PO_indices),1);
    apse_2 = cell(length(plot_PO_indices),1);
    
    lat_1 = cell(length(plot_PO_indices),1);

    parfor PO_i = 1:length(plot_PO_indices)
        index = plot_PO_indices(PO_i);

        [T_PO, X_PO, t_ev, X_ev, in_ev] = ode113(@Int_CR3BnSTM_J2pJ4pJ6pJ2s, [0, PO_data(index,c_Tp)], [PO_data(index,c_x0:c_zd0)'; stm0_colVec], options_apsis, prms_ZH);

        X_ev = X_ev((0.997 < X_ev(:,1)) & ( X_ev(:,1) < 0.999), 1:6);

        apse_1{PO_i} = [X_ev(1,1:6), PO_data(index,c_Tp)];
        apse_2{PO_i} = [X_ev(2,1:6), PO_data(index,c_Tp)];
        
        [lat_deg, lon_deg] = BCR2latlon(apse_1{PO_i}(1:3)', 'secondary', prms_ZH.u);
        lat_1{PO_i} = lat_deg;
    end
    
    apse_1 = cell2mat(apse_1);
    apse_2 = cell2mat(apse_2);
    
    lat_1 = cell2mat(lat_1);
    
    alt_apse1 = (rowNorm(apse_1(:,1:3) - [1-prms_ZH.u, 0, 0]) - secondary.R_n).*rNorm;
    alt_apse2 = (rowNorm(apse_2(:,1:3) - [1-prms_ZH.u, 0, 0]) - secondary.R_n).*rNorm;
    
    figure; plot(apse_1(:,7), alt_apse1, 'k', 'linewidth', 1.5)
    PlotBoi2('$T_P$', 'Altitude (km)', 26, 'LaTex')
	
    figure; plot(alt_apse1, apse_1(:,3), 'k', 'linewidth', 1.5)
    PlotBoi2('Altitude (km)', '$z_n$', 32, 'LaTex')
    xlim([0, 300])
%     figure; plot(alt_apse2)

    figure; plot(alt_apse1, lat_1, 'k', 'linewidth', 1.5)
    PlotBoi2('Altitude (km)', 'Latitude ($^\circ$)', 32, 'LaTex')
    xlim([0, 300])
end


if isequal(family, 'Jupiter_Europa.CR3BP_J2pJ4pJ6pJ2s.Se7_6P2_ns1.txt')
    %%% Preallocate
    apse_1 = cell(length(plot_PO_indices),1);
    apse_2 = cell(length(plot_PO_indices),1);
    
    lat_1 = cell(length(plot_PO_indices),1);

    parfor PO_i = 1:length(plot_PO_indices)
        index = plot_PO_indices(PO_i);

        [T_PO, X_PO, t_ev, X_ev, in_ev] = ode113(@Int_CR3BnSTM_J2pJ4pJ6pJ2s, [0, PO_data(index,c_Tp)], [PO_data(index,c_x0:c_zd0)'; stm0_colVec], options_apsis, prms_ZH);

        X_ev = X_ev((1.002 < X_ev(:,1)) & ( X_ev(:,1) < 1.003), 1:6);
        X_ev = X_ev((0 < X_ev(:,3)), 1:6);

        apse_1{PO_i} = [X_ev(1,1:6), PO_data(index,c_Tp)];
        apse_2{PO_i} = [X_ev(2,1:6), PO_data(index,c_Tp)];
        
        [lat_deg, lon_deg] = BCR2latlon(apse_1{PO_i}(1:3)', 'secondary', prms_ZH.u);
        lat_1{PO_i} = lat_deg;
    end
    
    apse_1 = cell2mat(apse_1);
    apse_2 = cell2mat(apse_2);
    
    lat_1 = cell2mat(lat_1);
    
    alt_apse1 = (rowNorm(apse_1(:,1:3) - [1-prms_ZH.u, 0, 0]) - secondary.R_n).*rNorm;
    alt_apse2 = (rowNorm(apse_2(:,1:3) - [1-prms_ZH.u, 0, 0]) - secondary.R_n).*rNorm;
    
    figure; plot(apse_1(:,7), alt_apse1, 'k', 'linewidth', 1.5)
    PlotBoi2('$T_P$', 'Altitude (km)', 26, 'LaTex')
	
    figure; plot(alt_apse1, apse_1(:,3), 'k', 'linewidth', 1.5)
    PlotBoi2('Altitude (km)', '$z_n$', 32, 'LaTex')
    xlim([0, 300])

    figure; plot(alt_apse1, lat_1, 'k', 'linewidth', 1.5)
    PlotBoi2('Altitude (km)', 'Latitude ($^\circ$)', 32, 'LaTex')
    xlim([0, 300])
end







if 1+1==1
    
    
    newIndex = 98;
    
    [T_PO_new, X_PO_new] = ode113(@Int_CR3BnSTM, [0, PO_data(newIndex,c_Tp)], [PO_data(newIndex,c_x0:c_zd0)'; stm0_colVec], options, prms);
    
    stm_tf_t0                           = reshape(X_PO_new(end,7:42),6,6);
    monodromy                           = stm_tf_t0;
    [eigenVectors_new, eigenValues_new] = eig(monodromy);
    [S1, S2]                            = getStabilityIndices(diag(eigenValues_new));
    S1
    S2
    diag(eigenValues_new)
    

    figure;
    hold all
    PlotBoi3_CR3Bn(26)
    axis equal
    plot3(rLPs_n(1:2,1),[0,0],[0,0],'^','markeredgecolor',colors.black, 'markerfacecolor', colors.blue)
    plotSecondary(secondary)
    
    plot3(X_PO_new(:,1), X_PO_new(:,2), X_PO_new(:,3), 'linewidth', lw, 'color', colors.blue2) % 0.1961    0.3922    0.7843


    latLons_PO = zeros(length(T_PO_new),2);
    for kk = 1:length(latLons_PO)
        [lat_deg, lon_deg] = BCR2latlon(X_PO_new(kk,1:3)', 'secondary', prms.u);
        latLons_PO(kk,:) = [lat_deg, lon_deg];
    end

    figure; hold all
    xlim([0 360])
    ylim([-90 90])
    PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')
    h = image(xlim, -ylim, secondary.img);
    [lons_new] = convert_lon180_to_lon360(latLons_PO(:,2));
    plot(lons_new,latLons_PO(:,1),'.', 'color', colors.blue2)


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
fprintf('\nElapsed time: %1.4f seconds\n',tocWhole)






