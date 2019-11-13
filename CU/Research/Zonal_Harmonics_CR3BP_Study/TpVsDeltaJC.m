% ========================================================================
%%% Description
% ========================================================================
% For comparing the relationship between time period and Jacobi constant
% for families of periodic orbits computed in 1) the vanilla CR3BP, and 2)
% the CR3BP with J2p, J4p, J6p, and J2s

% Created: 10/07/19
% Author : Luke Bury, luke.bury@colorado.edu
% ========================================================================
%%% Initialization
% ========================================================================
clear
clc
close all
mbinPath = '~/CU_Google_Drive/Documents/MatGit/mbin';
POfamilyPath = '~/CU_Google_Drive/Documents/MatGit/mbin/Data/InitialConditions/PO_Families/';
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
run_Tp_vs_deltaJC_for_verticals    = 1;
run_ZAmplitude_vs_Tp_for_verticals = 0;

% ========================================================================
%%% Setup
% ========================================================================
% -------------------------------------------------
%%% Column identifiers for family data file
% -------------------------------------------------
c_x0_n                = 1;
c_y0_n                = 2;
c_z0_n                = 3;
c_xd0_n               = 4;
c_yd0_n               = 5;
c_zd0_n               = 6;
c_Tp_n                = 7;
c_JC                  = 8;
c_L2FlythroughVel_mps = 9;
c_landingVelocity_mps = 10;

% -------------------------------------------------
%%% Family names
% -------------------------------------------------
%%% Jupiter-Europa
famName_JupEur_Lyap     = 'Jupiter_Europa.CR3BP.L2_Lyapunov';
famName_JupEur_Vert     = 'Jupiter_Europa.CR3BP.L2_Vertical';
famName_JupEur_SHalo    = 'Jupiter_Europa.CR3BP.L2_SHalo';

famName_JupEur_ZH_Lyap  = 'Jupiter_Europa.CR3BP_J2pJ4pJ6pJ2s.L2_Lyapunov';
famName_JupEur_ZH_Vert  = 'Jupiter_Europa.CR3BP_J2pJ4pJ6pJ2s.L2_Vertical';
famName_JupEur_ZH_SHalo = 'Jupiter_Europa.CR3BP_J2pJ4pJ6pJ2s.L2_SHalo';

%%% Saturn-Enceladus
famName_SatEnc_Lyap     = 'Saturn_Enceladus.CR3BP.L2_Lyapunov';
famName_SatEnc_Vert     = 'Saturn_Enceladus.CR3BP.L2_Vertical';
famName_SatEnc_SHalo    = 'Saturn_Enceladus.CR3BP.L2_SHalo';

famName_SatEnc_ZH_Lyap  = 'Saturn_Enceladus.CR3BP_J2pJ4pJ6pJ2s.L2_Lyapunov';
famName_SatEnc_ZH_Vert  = 'Saturn_Enceladus.CR3BP_J2pJ4pJ6pJ2s.L2_Vertical';
famName_SatEnc_ZH_SHalo = 'Saturn_Enceladus.CR3BP_J2pJ4pJ6pJ2s.L2_SHalo';



% -------------------------------------------------
%%% Set bodies
% -------------------------------------------------
jupiter   = bodies.jupiter;
europa    = bodies.europa;
saturn    = bodies.saturn;
enceladus = bodies.enceladus;

%%% Normalizing constants
rNorm_JupEur = europa.a;         % n <-> km
tNorm_JupEur = 1/europa.meanMot; % n <-> sec
vNorm_JupEur= rNorm_JupEur / tNorm_JupEur;       % n <-> km/sec

rNorm_SatEnc = enceladus.a;         % n <-> km
tNorm_SatEnc = 1/enceladus.meanMot; % n <-> sec
vNorm_SatEnc = rNorm_SatEnc / tNorm_SatEnc;       % n <-> km/sec
% -------------------------------------------------
%%% Load data 
% -------------------------------------------------
POFamilyData_JupEur_Lyap     = dlmread([POfamilyPath, famName_JupEur_Lyap, '.txt'],',',1,0);
POFamilyData_JupEur_Vert     = dlmread([POfamilyPath, famName_JupEur_Vert, '.txt'],',',1,0);
POFamilyData_JupEur_SHalo    = dlmread([POfamilyPath, famName_JupEur_SHalo, '.txt'],',',1,0);
POFamilyData_JupEur_Lyap_ZH  = dlmread([POfamilyPath, famName_JupEur_ZH_Lyap, '.txt'],',',1,0);
POFamilyData_JupEur_Vert_ZH  = dlmread([POfamilyPath, famName_JupEur_ZH_Vert, '.txt'],',',1,0);
POFamilyData_JupEur_SHalo_ZH = dlmread([POfamilyPath, famName_JupEur_ZH_SHalo, '.txt'],',',1,0);

POFamilyData_SatEnc_Lyap     = dlmread([POfamilyPath, famName_SatEnc_Lyap, '.txt'],',',1,0);
POFamilyData_SatEnc_Vert     = dlmread([POfamilyPath, famName_SatEnc_Vert, '.txt'],',',1,0);
POFamilyData_SatEnc_SHalo    = dlmread([POfamilyPath, famName_SatEnc_SHalo, '.txt'],',',1,0);
POFamilyData_SatEnc_Lyap_ZH  = dlmread([POfamilyPath, famName_SatEnc_ZH_Lyap, '.txt'],',',1,0);
POFamilyData_SatEnc_Vert_ZH  = dlmread([POfamilyPath, famName_SatEnc_ZH_Vert, '.txt'],',',1,0);
POFamilyData_SatEnc_SHalo_ZH = dlmread([POfamilyPath, famName_SatEnc_ZH_SHalo, '.txt'],',',1,0);

% -------------------------------------------------
%%% Determine JC of L2 in each system
% -------------------------------------------------
prms_JupEur.u = europa.MR;
prms_SatEnc.u = enceladus.MR;

rLp_JupEur = collinearEquilibriumPoints_ZH(prms_JupEur);
rLp_SatEnc = collinearEquilibriumPoints_ZH(prms_SatEnc);


prms_JupEur_ZH = prms_JupEur;
prms_SatEnc_ZH = prms_SatEnc;
prms_JupEur_ZH.J2p = jupiter.J2; prms_JupEur_ZH.J4p = jupiter.J4; prms_JupEur_ZH.J6p = jupiter.J6; prms_JupEur_ZH.J2s = europa.J2;    prms_JupEur_ZH.R1 = jupiter.R/rNorm_JupEur; prms_JupEur_ZH.R2 = europa.R_n;
prms_SatEnc_ZH.J2p = saturn.J2;  prms_SatEnc_ZH.J4p = saturn.J4;  prms_SatEnc_ZH.J6p = saturn.J6;  prms_SatEnc_ZH.J2s = enceladus.J2; prms_SatEnc_ZH.R1 = saturn.R/rNorm_SatEnc;  prms_SatEnc_ZH.R2 = enceladus.R_n;

rLp_JupEur_ZH = collinearEquilibriumPoints_ZH(prms_JupEur_ZH);
rLp_SatEnc_ZH = collinearEquilibriumPoints_ZH(prms_SatEnc_ZH);

rL2_JupEur    = rLp_JupEur(2,:);
rL2_JupEur_ZH = rLp_JupEur_ZH(2,:);

rL2_SatEnc    = rLp_SatEnc(2,:);
rL2_SatEnc_ZH = rLp_SatEnc_ZH(2,:);

JC_L2_JupEur = getJacobiConstant_ZH([rL2_JupEur, zeros(1,3)], prms_JupEur);
JC_L2_JupEur_ZH = getJacobiConstant_ZH([rL2_JupEur_ZH, zeros(1,3)], prms_JupEur_ZH);

JC_L2_SatEnc = getJacobiConstant_ZH([rL2_SatEnc, zeros(1,3)], prms_SatEnc);
JC_L2_SatEnc_ZH = getJacobiConstant_ZH([rL2_SatEnc_ZH, zeros(1,3)], prms_SatEnc_ZH);

% ========================================================================
%%% Plotting
% ========================================================================
if run_Tp_vs_deltaJC_for_verticals
    %%% dJC vectors
    dJC_JupEur_Lyap     = POFamilyData_JupEur_Lyap(:,c_JC)     - JC_L2_JupEur;
    dJC_JupEur_Lyap_ZH  = POFamilyData_JupEur_Lyap_ZH(:,c_JC)  - JC_L2_JupEur_ZH;
    dJC_JupEur_Vert     = POFamilyData_JupEur_Vert(:,c_JC)     - JC_L2_JupEur;
    dJC_JupEur_Vert_ZH  = POFamilyData_JupEur_Vert_ZH(:,c_JC)  - JC_L2_JupEur_ZH;
    dJC_JupEur_SHalo    = POFamilyData_JupEur_SHalo(:,c_JC)    - JC_L2_JupEur;
    dJC_JupEur_SHalo_ZH = POFamilyData_JupEur_SHalo_ZH(:,c_JC) - JC_L2_JupEur_ZH;
    
    dJC_SatEnc_Lyap     = POFamilyData_SatEnc_Lyap(:,c_JC)     - JC_L2_SatEnc;
    dJC_SatEnc_Lyap_ZH  = POFamilyData_SatEnc_Lyap_ZH(:,c_JC)  - JC_L2_SatEnc_ZH;
    dJC_SatEnc_Vert     = POFamilyData_SatEnc_Vert(:,c_JC)     - JC_L2_SatEnc;
    dJC_SatEnc_Vert_ZH  = POFamilyData_SatEnc_Vert_ZH(:,c_JC)  - JC_L2_SatEnc_ZH;
    dJC_SatEnc_SHalo    = POFamilyData_SatEnc_SHalo(:,c_JC)    - JC_L2_SatEnc;
    dJC_SatEnc_SHalo_ZH = POFamilyData_SatEnc_SHalo_ZH(:,c_JC) - JC_L2_SatEnc_ZH;
    
    
    %%% Time period bounds
    minTp_JupEur_Lyap  = max([min(POFamilyData_JupEur_Lyap(:, c_Tp_n)),  min(POFamilyData_JupEur_Lyap_ZH(:, c_Tp_n))]);
    maxTp_JupEur_Lyap  = min([max(POFamilyData_JupEur_Lyap(:, c_Tp_n)),  max(POFamilyData_JupEur_Lyap_ZH(:, c_Tp_n))]);
    minTp_JupEur_Vert  = max([min(POFamilyData_JupEur_Vert(:, c_Tp_n)),  min(POFamilyData_JupEur_Vert_ZH(:, c_Tp_n))]);
    maxTp_JupEur_Vert  = min([max(POFamilyData_JupEur_Vert(:, c_Tp_n)),  max(POFamilyData_JupEur_Vert_ZH(:, c_Tp_n))]);
    minTp_JupEur_SHalo = max([min(POFamilyData_JupEur_SHalo(:, c_Tp_n)), min(POFamilyData_JupEur_SHalo_ZH(:, c_Tp_n))]);
    maxTp_JupEur_SHalo = min([max(POFamilyData_JupEur_SHalo(:, c_Tp_n)), max(POFamilyData_JupEur_SHalo_ZH(:, c_Tp_n))]);
    
    minTp_SatEnc_Lyap  = max([min(POFamilyData_SatEnc_Lyap(:, c_Tp_n)),  min(POFamilyData_SatEnc_Lyap_ZH(:, c_Tp_n))]);
    maxTp_SatEnc_Lyap  = min([max(POFamilyData_SatEnc_Lyap(:, c_Tp_n)),  max(POFamilyData_SatEnc_Lyap_ZH(:, c_Tp_n))]);
    minTp_SatEnc_Vert  = max([min(POFamilyData_SatEnc_Vert(:, c_Tp_n)),  min(POFamilyData_SatEnc_Vert_ZH(:, c_Tp_n))]);
    maxTp_SatEnc_Vert  = min([max(POFamilyData_SatEnc_Vert(:, c_Tp_n)),  max(POFamilyData_SatEnc_Vert_ZH(:, c_Tp_n))]);
    minTp_SatEnc_SHalo = max([min(POFamilyData_SatEnc_SHalo(:, c_Tp_n)), min(POFamilyData_SatEnc_SHalo_ZH(:, c_Tp_n))]);
    maxTp_SatEnc_SHalo = min([max(POFamilyData_SatEnc_SHalo(:, c_Tp_n)), max(POFamilyData_SatEnc_SHalo_ZH(:, c_Tp_n))]);
    
    %%% Time-period ranges
    nInterpPoints = 100;
    TpRange_JupEur_Lyap  = linspace(minTp_JupEur_Lyap,  maxTp_JupEur_Lyap,  nInterpPoints);
    TpRange_JupEur_Vert  = linspace(minTp_JupEur_Vert,  maxTp_JupEur_Vert,  nInterpPoints);
    TpRange_JupEur_SHalo = linspace(minTp_JupEur_SHalo, maxTp_JupEur_SHalo, nInterpPoints);
    TpRange_SatEnc_Lyap  = linspace(minTp_SatEnc_Lyap,  maxTp_SatEnc_Lyap,  nInterpPoints);
    TpRange_SatEnc_Vert  = linspace(minTp_SatEnc_Vert,  maxTp_SatEnc_Vert,  nInterpPoints);
    TpRange_SatEnc_SHalo = linspace(minTp_SatEnc_SHalo, maxTp_SatEnc_SHalo, nInterpPoints);
    
    %%% Calculate interpolated data
    interpolated_dJCs_JupEur_Lyap     = interp1(POFamilyData_JupEur_Lyap(:, c_Tp_n),     dJC_JupEur_Lyap,     TpRange_JupEur_Lyap);
    interpolated_dJCs_JupEur_Lyap_ZH  = interp1(POFamilyData_JupEur_Lyap_ZH(:, c_Tp_n),  dJC_JupEur_Lyap_ZH,  TpRange_JupEur_Lyap);
    interpolated_dJCs_JupEur_Vert     = interp1(POFamilyData_JupEur_Vert(:, c_Tp_n),     dJC_JupEur_Vert,     TpRange_JupEur_Vert);
    interpolated_dJCs_JupEur_Vert_ZH  = interp1(POFamilyData_JupEur_Vert_ZH(:, c_Tp_n),  dJC_JupEur_Vert_ZH,  TpRange_JupEur_Vert);
    interpolated_dJCs_JupEur_SHalo    = interp1(POFamilyData_JupEur_SHalo(:, c_Tp_n),    dJC_JupEur_SHalo,    TpRange_JupEur_SHalo);
    interpolated_dJCs_JupEur_SHalo_ZH = interp1(POFamilyData_JupEur_SHalo_ZH(:, c_Tp_n), dJC_JupEur_SHalo_ZH, TpRange_JupEur_SHalo);
    
    interpolated_dJCs_SatEnc_Lyap     = interp1(POFamilyData_SatEnc_Lyap(:, c_Tp_n),     dJC_SatEnc_Lyap,     TpRange_SatEnc_Lyap);
    interpolated_dJCs_SatEnc_Lyap_ZH  = interp1(POFamilyData_SatEnc_Lyap_ZH(:, c_Tp_n),  dJC_SatEnc_Lyap_ZH,  TpRange_SatEnc_Lyap);
    interpolated_dJCs_SatEnc_Vert     = interp1(POFamilyData_SatEnc_Vert(:, c_Tp_n),     dJC_SatEnc_Vert,     TpRange_SatEnc_Vert);
    interpolated_dJCs_SatEnc_Vert_ZH  = interp1(POFamilyData_SatEnc_Vert_ZH(:, c_Tp_n),  dJC_SatEnc_Vert_ZH,  TpRange_SatEnc_Vert);
    interpolated_dJCs_SatEnc_SHalo    = interp1(POFamilyData_SatEnc_SHalo(:, c_Tp_n),    dJC_SatEnc_SHalo,    TpRange_SatEnc_SHalo);
    interpolated_dJCs_SatEnc_SHalo_ZH = interp1(POFamilyData_SatEnc_SHalo_ZH(:, c_Tp_n), dJC_SatEnc_SHalo_ZH, TpRange_SatEnc_SHalo);
    
    %%% Difference interpolated data
    delta_dJCs_JupEur_Lyap  = interpolated_dJCs_JupEur_Lyap_ZH  - interpolated_dJCs_JupEur_Lyap;
    delta_dJCs_JupEur_Vert  = interpolated_dJCs_JupEur_Vert_ZH  - interpolated_dJCs_JupEur_Vert;
    delta_dJCs_JupEur_SHalo = interpolated_dJCs_JupEur_SHalo_ZH - interpolated_dJCs_JupEur_SHalo;
    
    delta_dJCs_SatEnc_Lyap  = interpolated_dJCs_SatEnc_Lyap_ZH  - interpolated_dJCs_SatEnc_Lyap;
    delta_dJCs_SatEnc_Vert  = interpolated_dJCs_SatEnc_Vert_ZH  - interpolated_dJCs_SatEnc_Vert;
    delta_dJCs_SatEnc_SHalo = interpolated_dJCs_SatEnc_SHalo_ZH - interpolated_dJCs_SatEnc_SHalo;
    
    
%     figure; hold all
%     p1 = plot(dJC_JupEur_Lyap, POFamilyData_JupEur_Lyap(:, c_Tp_n),'r', 'linewidth', 2);
%     p2 = plot(dJC_JupEur_Lyap_ZH, POFamilyData_JupEur_Lyap_ZH(:, c_Tp_n),'b', 'linewidth', 2);
%     PlotBoi2('$JC$ - $JC_{L_2}$', 'Normalized Time Period',18,'LaTex')
%     legend([p1, p2],'CR3BP','CR3BP w ZH')
%     % title('Jupiter-Europa, L_2 Lyapunov')
% 
% 
%     figure; hold all
%     p1 = plot(dJC_JupEur_Vert, POFamilyData_JupEur_Vert(:, c_Tp_n),'r', 'linewidth', 2);
%     p2 = plot(dJC_JupEur_Vert_ZH, POFamilyData_JupEur_Vert_ZH(:, c_Tp_n),'b', 'linewidth', 2);
%     PlotBoi2('$JC$ - $JC_{L_2}$', 'Normalized Time Period',18,'LaTex')
%     legend([p1, p2],'CR3BP','CR3BP w ZH')
%     % title('Jupiter-Europa, L_2 Vertical')
% 
%     figure; hold all
%     p1 = plot(dJC_JupEur_SHalo, POFamilyData_JupEur_SHalo(:, c_Tp_n),'r', 'linewidth', 2);
%     p2 = plot(dJC_JupEur_SHalo_ZH, POFamilyData_JupEur_SHalo_ZH(:, c_Tp_n),'b', 'linewidth', 2);
%     PlotBoi2('$JC$ - $JC_{L_2}$', 'Normalized Time Period',18,'LaTex')
%     legend([p1, p2],'CR3BP','CR3BP w ZH')
%     % title('Jupiter-Europa, L_2 Northern Halo')
% 
% 
%     figure; hold all
%     p1 = plot(dJC_SatEnc_Lyap, POFamilyData_SatEnc_Lyap(:, c_Tp_n),'r', 'linewidth', 2);
%     p2 = plot(dJC_SatEnc_Lyap_ZH, POFamilyData_SatEnc_Lyap_ZH(:, c_Tp_n),'b', 'linewidth', 2);
%     PlotBoi2('$JC$ - $JC_{L_2}$', 'Normalized Time Period',18,'LaTex')
%     legend([p1, p2],'CR3BP','CR3BP w ZH')
%     % title('Saturn-Enceladus, L_2 Lyapunov')
% 
% 
%     figure; hold all
%     p1 = plot(dJC_SatEnc_Vert, POFamilyData_SatEnc_Vert(:, c_Tp_n),'r', 'linewidth', 2);
%     p2 = plot(dJC_SatEnc_Vert_ZH, POFamilyData_SatEnc_Vert_ZH(:, c_Tp_n),'b', 'linewidth', 2);
%     PlotBoi2('$JC$ - $JC_{L_2}$', 'Normalized Time Period',18,'LaTex')
%     legend([p1, p2],'CR3BP','CR3BP w ZH')
%     % title('Saturn-Enceladus, L_2 Vertical')

%     figure; hold all
%     p1 = plot(dJC_SatEnc_SHalo, POFamilyData_SatEnc_SHalo(:, c_Tp_n),'r', 'linewidth', 2);
%     p2 = plot(dJC_SatEnc_SHalo_ZH, POFamilyData_SatEnc_SHalo_ZH(:, c_Tp_n),'b', 'linewidth', 2);
%     PlotBoi2('$JC$ - $JC_{L_2}$', 'Normalized Time Period',18,'LaTex')
%     legend([p1, p2],'CR3BP','CR3BP w ZH')
%     % title('Saturn-Enceladus, L_2 Northern Halo')
    
    figure; hold all
    plot(TpRange_JupEur_Lyap, delta_dJCs_JupEur_Lyap, 'm', 'linewidth', 2)
    PlotBoi2('Normalized Time Period', '$\Delta$($JC_{Lyap}$ - $JC_{L_2}$)',18,'LaTex')
    figure; hold all
    plot(TpRange_JupEur_Vert, delta_dJCs_JupEur_Vert, 'm', 'linewidth', 2)
    PlotBoi2('Normalized Time Period', '$\Delta$($JC_{Vert}$ - $JC_{L_2}$)',18,'LaTex')
    figure; hold all
    plot(TpRange_JupEur_SHalo, delta_dJCs_JupEur_SHalo, 'm', 'linewidth', 2)
    PlotBoi2('Normalized Time Period', '$\Delta$($JC_{SHalo}$ - $JC_{L_2}$)',18,'LaTex')
    
    figure; hold all
    plot(TpRange_SatEnc_Lyap, delta_dJCs_SatEnc_Lyap, 'm', 'linewidth', 2)
    PlotBoi2('Normalized Time Period', '$\Delta$($JC_{Lyap}$ - $JC_{L_2}$)',18,'LaTex')
    figure; hold all
    plot(TpRange_SatEnc_Vert, delta_dJCs_SatEnc_Vert, 'm', 'linewidth', 2)
    PlotBoi2('Normalized Time Period', '$\Delta$($JC_{Vert}$ - $JC_{L_2}$)',18,'LaTex')
    figure; hold all
    plot(TpRange_SatEnc_SHalo, delta_dJCs_SatEnc_SHalo, 'm', 'linewidth', 2)
    PlotBoi2('Normalized Time Period', '$\Delta$($JC_{SHalo}$ - $JC_{L_2}$)',18,'LaTex')
    


end % run_Tp_vs_deltaJC_for_verticals


% ========================================================================
%%% Comparing amplitude of vertical orbits to their period
% ========================================================================
if run_ZAmplitude_vs_Tp_for_verticals
    %%% Integration options
    tol     = 1e-13;
    options = odeset('RelTol',tol,'AbsTol',tol);
    
    
    amp_JupEur = NaN(size(POFamilyData_JupEur_Vert,1),1);
    for kk = 1:size(POFamilyData_JupEur_Vert,1)
        X0_n = POFamilyData_JupEur_Vert(kk, c_x0_n:c_zd0_n);
        Tp_n = POFamilyData_JupEur_Vert(kk, c_Tp_n);
        
        [~, X_JupEur] = ode113(@Int_CR3Bn, [0, Tp_n], X0_n', options, prms_JupEur);
        
        amp_JupEur(kk) = max(X_JupEur(:,3));
    end
    
    amp_JupEur_ZH = NaN(size(POFamilyData_JupEur_Vert_ZH,1),1);
    for kk = 1:size(POFamilyData_JupEur_Vert_ZH,1)
        X0_n = POFamilyData_JupEur_Vert_ZH(kk, c_x0_n:c_zd0_n);
        Tp_n = POFamilyData_JupEur_Vert_ZH(kk, c_Tp_n);
        
        [~, X_JupEur_ZH] = ode113(@Int_CR3Bn_J2pJ4pJ6pJ2s, [0, Tp_n], X0_n', options, prms_JupEur_ZH);
        
        amp_JupEur_ZH(kk) = max(X_JupEur_ZH(:,3));
    end
    
    amp_SatEnc = NaN(size(POFamilyData_SatEnc_Vert,1),1);
    for kk = 1:size(POFamilyData_SatEnc_Vert,1)
        X0_n = POFamilyData_SatEnc_Vert(kk, c_x0_n:c_zd0_n);
        Tp_n = POFamilyData_SatEnc_Vert(kk, c_Tp_n);
        
        [~, X_SatEnc] = ode113(@Int_CR3Bn, [0, Tp_n], X0_n', options, prms_SatEnc);
        
        amp_SatEnc(kk) = max(X_SatEnc(:,3));
    end
    
    amp_SatEnc_ZH = NaN(size(POFamilyData_SatEnc_Vert_ZH,1),1);
    for kk = 1:size(POFamilyData_SatEnc_Vert_ZH,1)
        X0_n = POFamilyData_SatEnc_Vert_ZH(kk, c_x0_n:c_zd0_n);
        Tp_n = POFamilyData_SatEnc_Vert_ZH(kk, c_Tp_n);
        
        [~, X_SatEnc_ZH] = ode113(@Int_CR3Bn_J2pJ4pJ6pJ2s, [0, Tp_n], X0_n', options, prms_SatEnc_ZH);
        
        amp_SatEnc_ZH(kk) = max(X_SatEnc_ZH(:,3));
    end
    
    figure; hold all
    p2 = plot(amp_JupEur_ZH, POFamilyData_JupEur_Vert_ZH(:, c_Tp_n), 'b', 'linewidth', 2);
    p1 = plot(amp_JupEur, POFamilyData_JupEur_Vert(:, c_Tp_n), 'r', 'linewidth', 2);
    PlotBoi2('Amplitude of Vertical Orbit','Normalized Time Period',18,'LaTex')
    legend([p1, p2],'CR3BP','CR3BP w ZH')
    
    figure; hold all
    p2 = plot(amp_SatEnc_ZH, POFamilyData_SatEnc_Vert_ZH(:, c_Tp_n), 'b', 'linewidth', 2);
    p1 = plot(amp_SatEnc, POFamilyData_SatEnc_Vert(:, c_Tp_n), 'r', 'linewidth', 2);
    PlotBoi2('Amplitude of Vertical Orbit','Normalized Time Period',18,'LaTex')
    legend([p1, p2],'CR3BP','CR3BP w ZH')
    
%     figure; hold all
%     p2 = plot(amp_JupEur_ZH, POFamilyData_JupEur_ZH_Vert(:, c_JC) - JC_L2_JupEur_ZH, 'b', 'linewidth', 2);
%     p1 = plot(amp_JupEur, POFamilyData_JupEur_Vert(:, c_JC) - JC_L2_JupEur, 'r', 'linewidth', 2);
%     PlotBoi2('Amplitude of Vertical Orbit','$JC$ - $JC_{L_2}$',18,'LaTex')
%     legend([p1, p2],'CR3BP','CR3BP w ZH')
%     
%     figure; hold all
%     p2 = plot(amp_SatEnc_ZH, POFamilyData_SatEnc_ZH_Vert(:, c_JC) - JC_L2_SatEnc_ZH, 'b', 'linewidth', 2);
%     p1 = plot(amp_SatEnc, POFamilyData_SatEnc_Vert(:, c_JC) - JC_L2_SatEnc, 'r', 'linewidth', 2);
%     PlotBoi2('Amplitude of Vertical Orbit','$JC$ - $JC_{L_2}$',18,'LaTex')
%     legend([p1, p2],'CR3BP','CR3BP w ZH')
    


%
    % -------------------------------------------------
    %%% Tp vs delta-amplitude for vertical orbits
    % -------------------------------------------------
    %%% Get bounds of all parameters for interpolation
    minAmp_JupEur = max([min(amp_JupEur), min(amp_JupEur_ZH)]);
    maxAmp_JupEur = min([max(amp_JupEur), min(amp_JupEur_ZH)]);
    
    minAmp_SatEnc = max([min(amp_SatEnc), min(amp_SatEnc_ZH)]);
    maxAmp_SatEnc = min([max(amp_SatEnc), min(amp_SatEnc_ZH)]);
    
    minTp_JupEur = max([min(POFamilyData_JupEur_Vert(:, c_Tp_n)), min(POFamilyData_JupEur_Vert_ZH(:, c_Tp_n))]);
    maxTp_JupEur = min([max(POFamilyData_JupEur_Vert(:, c_Tp_n)), max(POFamilyData_JupEur_Vert_ZH(:, c_Tp_n))]);
    
    minTp_SatEnc = max([min(POFamilyData_SatEnc_Vert(:, c_Tp_n)), min(POFamilyData_SatEnc_Vert_ZH(:, c_Tp_n))]);
    maxTp_SatEnc = min([max(POFamilyData_SatEnc_Vert(:, c_Tp_n)), max(POFamilyData_SatEnc_Vert_ZH(:, c_Tp_n))]);
    
    %%% Time-period ranges
    nInterpPoints = 100;
    TpRange_JupEur = linspace(minTp_JupEur, maxTp_JupEur, nInterpPoints);
    TpRange_SatEnc = linspace(minTp_SatEnc, maxTp_SatEnc, nInterpPoints);
    
    %%% Calculate interpolated data
    interpolatedAmps_JupEur    = interp1(POFamilyData_JupEur_Vert(:, c_Tp_n), amp_JupEur, TpRange_JupEur);
    interpolatedAmps_JupEur_ZH = interp1(POFamilyData_JupEur_Vert_ZH(:, c_Tp_n), amp_JupEur_ZH, TpRange_JupEur);
    
    interpolatedAmps_SatEnc    = interp1(POFamilyData_SatEnc_Vert(:, c_Tp_n), amp_SatEnc, TpRange_SatEnc);
    interpolatedAmps_SatEnc_ZH = interp1(POFamilyData_SatEnc_Vert_ZH(:, c_Tp_n), amp_SatEnc_ZH, TpRange_SatEnc);
    
    deltaAmp_JupEur = interpolatedAmps_JupEur_ZH - interpolatedAmps_JupEur;
    deltaAmp_SatEnc = interpolatedAmps_SatEnc_ZH - interpolatedAmps_SatEnc;
    
    figure; hold all
    plot(TpRange_JupEur, deltaAmp_JupEur, 'm', 'linewidth', 2)
    PlotBoi2('Normalized Time Period','$\Delta$Amplitude of Vertical Orbit',18,'LaTex')
    
    figure; hold all
    plot(TpRange_SatEnc, deltaAmp_SatEnc, 'm', 'linewidth', 2)
    PlotBoi2('Normalized Time Period','$\Delta$Amplitude of Vertical Orbit',18,'LaTex')
    
    
% % %     correspondingTps    = interp1(amp_SatEnc,    POFamilyData_SatEnc_Vert(:, c_Tp_n),    linspace(min(amp_SatEnc), max(amp_SatEnc), length(amp_SatEnc)));
% % %     correspondingTps_ZH = interp1(amp_SatEnc_ZH, POFamilyData_SatEnc_ZH_Vert(:, c_Tp_n), linspace(min(amp_SatEnc), max(amp_SatEnc), length(amp_SatEnc)));
% % %     
% % %     correspondingAmps   = interp1(POFamilyData_SatEnc_Vert(:, c_Tp_n),    amp_SatEnc,    linspace(min(POFamilyData_SatEnc_ZH_Vert(:, c_Tp_n)), max(POFamilyData_SatEnc_ZH_Vert(:, c_Tp_n)), length(POFamilyData_SatEnc_ZH_Vert(:, c_Tp_n))));
% % %     correspondingAmps_ZH = interp1(POFamilyData_SatEnc_ZH_Vert(:, c_Tp_n),    amp_SatEnc_ZH,    linspace(min(POFamilyData_SatEnc_ZH_Vert(:, c_Tp_n)), max(POFamilyData_SatEnc_ZH_Vert(:, c_Tp_n)), length(POFamilyData_SatEnc_ZH_Vert(:, c_Tp_n))));
% % % %     figure; hold all
% % % %     plot(linspace(min(amp_SatEnc), max(amp_SatEnc), length(amp_SatEnc)), correspondingTps_ZH - correspondingTps, 'm', 'linewidth', 2);
% % % %     PlotBoi2('Amplitude of Vertical Orbit','$\Delta$Normalized Time Period',18,'LaTex')
% % %     figure; hold all
% % %     plot(linspace(min(POFamilyData_SatEnc_ZH_Vert(:, c_Tp_n)), max(POFamilyData_SatEnc_ZH_Vert(:, c_Tp_n)), length(POFamilyData_SatEnc_ZH_Vert(:, c_Tp_n))), correspondingAmps_ZH - correspondingAmps, 'm', 'linewidth', 2);
% % %     PlotBoi2('Normalized Time Period','$\Delta$Amplitude of Vertical Orbit',18,'LaTex')%     legend([p1, p2],'CR3BP','CR3BP w ZH')
end
% ========================================================================
%%% Formatting Structures
% ========================================================================
% -------------------------------------------------
%%% 
% -------------------------------------------------

% --------------------------
%%% 
% --------------------------

% ========================================================================
%%% Closeout
% ========================================================================
tocWhole = toc(ticWhole);
fprintf('\n(Elapsed time: %1.4f seconds)\n',tocWhole)













