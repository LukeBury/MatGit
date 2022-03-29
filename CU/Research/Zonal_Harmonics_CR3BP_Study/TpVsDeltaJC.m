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
dataSetPath = '~/CU_Google_Drive/Documents/MatGit/CU/Research/Zonal_Harmonics_CR3BP_Study/ZH_POs_for_direct_comparison/';
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
run_ZAmplitude_vs_Tp_for_verticals = 1;

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
% %%% Jupiter-Europa
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
tNorm_JupEur = sqrt((rNorm_JupEur^3)/(bodies.constants.G*(jupiter.mass + europa.mass))); % n <-> sec

rNorm_SatEnc = enceladus.a;         % n <-> km
tNorm_SatEnc = sqrt((rNorm_SatEnc^3)/(bodies.constants.G*(saturn.mass + enceladus.mass))); % n <-> sec
% -------------------------------------------------
%%% Load data 
% -------------------------------------------------
POFamilyData_JupEur_Lyap     = dlmread([dataSetPath, famName_JupEur_Lyap, '.300.txt'],',',1,0);
POFamilyData_JupEur_Vert     = dlmread([dataSetPath, famName_JupEur_Vert, '.300.txt'],',',1,0);
POFamilyData_JupEur_SHalo    = dlmread([dataSetPath, famName_JupEur_SHalo, '.300.txt'],',',1,0);
POFamilyData_JupEur_Lyap_ZH  = dlmread([dataSetPath, famName_JupEur_ZH_Lyap, '.300.txt'],',',1,0);
POFamilyData_JupEur_Vert_ZH  = dlmread([dataSetPath, famName_JupEur_ZH_Vert, '.300.txt'],',',1,0);
POFamilyData_JupEur_SHalo_ZH = dlmread([dataSetPath, famName_JupEur_ZH_SHalo, '.300.txt'],',',1,0);

POFamilyData_SatEnc_Lyap     = dlmread([dataSetPath, famName_SatEnc_Lyap, '.300.txt'],',',1,0);
POFamilyData_SatEnc_Vert     = dlmread([dataSetPath, famName_SatEnc_Vert, '.300.txt'],',',1,0);
POFamilyData_SatEnc_SHalo    = dlmread([dataSetPath, famName_SatEnc_SHalo, '.300.txt'],',',1,0);
POFamilyData_SatEnc_Lyap_ZH  = dlmread([dataSetPath, famName_SatEnc_ZH_Lyap, '.300.txt'],',',1,0);
POFamilyData_SatEnc_Vert_ZH  = dlmread([dataSetPath, famName_SatEnc_ZH_Vert, '.300.txt'],',',1,0);
POFamilyData_SatEnc_SHalo_ZH = dlmread([dataSetPath, famName_SatEnc_ZH_SHalo, '.300.txt'],',',1,0);

% -------------------------------------------------
%%% Determine JC of L2 in each system
% -------------------------------------------------
prms_JupEur.u = europa.MR;
prms_SatEnc.u = enceladus.MR;

prms_JupEur.n = 1;
prms_SatEnc.n = 1;

rLp_JupEur = collinearEquilibriumPoints_ZH(prms_JupEur);
rLp_SatEnc = collinearEquilibriumPoints_ZH(prms_SatEnc);


prms_JupEur_ZH = prms_JupEur;
prms_SatEnc_ZH = prms_SatEnc;
prms_JupEur_ZH.J2p = jupiter.J2; prms_JupEur_ZH.J4p = jupiter.J4; prms_JupEur_ZH.J6p = jupiter.J6; prms_JupEur_ZH.J2s = europa.J2;    prms_JupEur_ZH.R1 = jupiter.R/rNorm_JupEur; prms_JupEur_ZH.R2 = europa.R_n;
prms_SatEnc_ZH.J2p = saturn.J2;  prms_SatEnc_ZH.J4p = saturn.J4;  prms_SatEnc_ZH.J6p = saturn.J6;  prms_SatEnc_ZH.J2s = enceladus.J2; prms_SatEnc_ZH.R1 = saturn.R/rNorm_SatEnc;  prms_SatEnc_ZH.R2 = enceladus.R_n;

prms_JupEur_ZH.n = europa.meanMot    * tNorm_JupEur;
prms_SatEnc_ZH.n = enceladus.meanMot * tNorm_SatEnc;

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
    %%% Differences in jacobi constant between each member of family and L2
    %%% of that system
    dJCL2_JupEur_Lyap     = POFamilyData_JupEur_Lyap(:,c_JC)     - JC_L2_JupEur;
    dJCL2_JupEur_Lyap_ZH  = POFamilyData_JupEur_Lyap_ZH(:,c_JC)  - JC_L2_JupEur_ZH;
    dJCL2_JupEur_Vert     = POFamilyData_JupEur_Vert(:,c_JC)     - JC_L2_JupEur;
    dJCL2_JupEur_Vert_ZH  = POFamilyData_JupEur_Vert_ZH(:,c_JC)  - JC_L2_JupEur_ZH;
    dJCL2_JupEur_SHalo    = POFamilyData_JupEur_SHalo(:,c_JC)    - JC_L2_JupEur;
    dJCL2_JupEur_SHalo_ZH = POFamilyData_JupEur_SHalo_ZH(:,c_JC) - JC_L2_JupEur_ZH;
    
    dJCL2_SatEnc_Lyap     = POFamilyData_SatEnc_Lyap(:,c_JC)     - JC_L2_SatEnc;
    dJCL2_SatEnc_Lyap_ZH  = POFamilyData_SatEnc_Lyap_ZH(:,c_JC)  - JC_L2_SatEnc_ZH;
    dJCL2_SatEnc_Vert     = POFamilyData_SatEnc_Vert(:,c_JC)     - JC_L2_SatEnc;
    dJCL2_SatEnc_Vert_ZH  = POFamilyData_SatEnc_Vert_ZH(:,c_JC)  - JC_L2_SatEnc_ZH;
    dJCL2_SatEnc_SHalo    = POFamilyData_SatEnc_SHalo(:,c_JC)    - JC_L2_SatEnc;
    dJCL2_SatEnc_SHalo_ZH = POFamilyData_SatEnc_SHalo_ZH(:,c_JC) - JC_L2_SatEnc_ZH;
    
    
    %%% Time period bounds that are shared between the classical and ZH
    %%% families
    TpRange_JupEur_Lyap  = POFamilyData_JupEur_Lyap(:,c_Tp_n);
    TpRange_JupEur_Vert  = POFamilyData_JupEur_Vert(:,c_Tp_n);
    TpRange_JupEur_SHalo = POFamilyData_JupEur_SHalo(:,c_Tp_n);
    TpRange_SatEnc_Lyap  = POFamilyData_SatEnc_Lyap(:,c_Tp_n);
    TpRange_SatEnc_Vert  = POFamilyData_SatEnc_Vert(:,c_Tp_n);
    TpRange_SatEnc_SHalo = POFamilyData_SatEnc_SHalo(:,c_Tp_n);    

    %%% Difference interpolated data
    delta_dJCs_JupEur_Lyap  = dJCL2_JupEur_Lyap_ZH  - dJCL2_JupEur_Lyap;
    delta_dJCs_JupEur_Vert  = dJCL2_JupEur_Vert_ZH  - dJCL2_JupEur_Vert;
    delta_dJCs_JupEur_SHalo = dJCL2_JupEur_SHalo_ZH - dJCL2_JupEur_SHalo;
    
    delta_dJCs_SatEnc_Lyap  = dJCL2_SatEnc_Lyap_ZH  - dJCL2_SatEnc_Lyap;
    delta_dJCs_SatEnc_Vert  = dJCL2_SatEnc_Vert_ZH  - dJCL2_SatEnc_Vert;
    delta_dJCs_SatEnc_SHalo = dJCL2_SatEnc_SHalo_ZH - dJCL2_SatEnc_SHalo;
    

    
    figure; hold all
    plot(TpRange_JupEur_Lyap, delta_dJCs_JupEur_Lyap, 'm', 'linewidth', 2)
    PlotBoi2('Normalized Time Period', '$\Delta$($JC_{Lyap}$ - $JC_{L_2}$)',20,'LaTex')
    figure; hold all
    plot(TpRange_JupEur_Vert, delta_dJCs_JupEur_Vert, 'm', 'linewidth', 2)
    PlotBoi2('Normalized Time Period', '$\Delta$($JC_{Vert}$ - $JC_{L_2}$)',20,'LaTex')
    figure; hold all
    plot(TpRange_JupEur_SHalo, delta_dJCs_JupEur_SHalo, 'm', 'linewidth', 2)
    PlotBoi2('Normalized Time Period', '$\Delta$($JC_{SHalo}$ - $JC_{L_2}$)',20,'LaTex')
    
    figure; hold all
    plot(TpRange_SatEnc_Lyap, delta_dJCs_SatEnc_Lyap, 'm', 'linewidth', 2)
    PlotBoi2('Normalized Time Period', '$\Delta$($JC_{Lyap}$ - $JC_{L_2}$)',20,'LaTex')
    figure; hold all
    plot(TpRange_SatEnc_Vert, delta_dJCs_SatEnc_Vert, 'm', 'linewidth', 2)
    PlotBoi2('Normalized Time Period', '$\Delta$($JC_{Vert}$ - $JC_{L_2}$)',20,'LaTex')
    figure; hold all
    plot(TpRange_SatEnc_SHalo, delta_dJCs_SatEnc_SHalo, 'm', 'linewidth', 2)
    PlotBoi2('Normalized Time Period', '$\Delta$($JC_{SHalo}$ - $JC_{L_2}$)',20,'LaTex')
    


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
        
        [~, X_JupEur] = ode113(@Int_CR3Bn, linspace(0, Tp_n, 10000), X0_n', options, prms_JupEur);
        
        amp_JupEur(kk) = max(X_JupEur(:,3));
    end
    
    amp_JupEur_ZH = NaN(size(POFamilyData_JupEur_Vert_ZH,1),1);
    for kk = 1:size(POFamilyData_JupEur_Vert_ZH,1)
        X0_n = POFamilyData_JupEur_Vert_ZH(kk, c_x0_n:c_zd0_n);
        Tp_n = POFamilyData_JupEur_Vert_ZH(kk, c_Tp_n);
        
        [~, X_JupEur_ZH] = ode113(@Int_CR3Bn_J2pJ4pJ6pJ2s, linspace(0, Tp_n, 10000), X0_n', options, prms_JupEur_ZH);
        
        amp_JupEur_ZH(kk) = max(X_JupEur_ZH(:,3));
    end
    
    amp_SatEnc = NaN(size(POFamilyData_SatEnc_Vert,1),1);
    for kk = 1:size(POFamilyData_SatEnc_Vert,1)
        X0_n = POFamilyData_SatEnc_Vert(kk, c_x0_n:c_zd0_n);
        Tp_n = POFamilyData_SatEnc_Vert(kk, c_Tp_n);
        
        [~, X_SatEnc] = ode113(@Int_CR3Bn, linspace(0, Tp_n, 10000), X0_n', options, prms_SatEnc);
        
        amp_SatEnc(kk) = max(X_SatEnc(:,3));
    end
    
    amp_SatEnc_ZH = NaN(size(POFamilyData_SatEnc_Vert_ZH,1),1);
    for kk = 1:size(POFamilyData_SatEnc_Vert_ZH,1)
        X0_n = POFamilyData_SatEnc_Vert_ZH(kk, c_x0_n:c_zd0_n);
        Tp_n = POFamilyData_SatEnc_Vert_ZH(kk, c_Tp_n);
        
        [~, X_SatEnc_ZH] = ode113(@Int_CR3Bn_J2pJ4pJ6pJ2s, linspace(0, Tp_n, 10000), X0_n', options, prms_SatEnc_ZH);
        
        amp_SatEnc_ZH(kk) = max(X_SatEnc_ZH(:,3));
    end
    
%     figure; hold all
%     p2 = plot(amp_JupEur_ZH, POFamilyData_JupEur_Vert_ZH(:, c_Tp_n), 'b', 'linewidth', 2);
%     p1 = plot(amp_JupEur, POFamilyData_JupEur_Vert(:, c_Tp_n), 'r', 'linewidth', 2);
%     PlotBoi2('Amplitude of Vertical Orbit','Normalized Time Period',20,'LaTex')
%     legend([p1, p2],'CR3BP','CR3BP w ZH')
%     
%     figure; hold all
%     p2 = plot(amp_SatEnc_ZH, POFamilyData_SatEnc_Vert_ZH(:, c_Tp_n), 'b', 'linewidth', 2);
%     p1 = plot(amp_SatEnc, POFamilyData_SatEnc_Vert(:, c_Tp_n), 'r', 'linewidth', 2);
%     PlotBoi2('Amplitude of Vertical Orbit','Normalized Time Period',20,'LaTex')
%     legend([p1, p2],'CR3BP','CR3BP w ZH')


    % -------------------------------------------------
    %%% Tp vs delta-amplitude for vertical orbits
    % -------------------------------------------------
    deltaAmp_JupEur = amp_JupEur_ZH - amp_JupEur;
    deltaAmp_SatEnc = amp_SatEnc_ZH - amp_SatEnc;
    
    figure; hold all
    plot(TpRange_JupEur_Vert, deltaAmp_JupEur, 'm', 'linewidth', 2)
    PlotBoi2('Normalized Time Period','$\Delta$Amplitude of Vertical Orbit',20,'LaTex')
    
    figure; hold all
    plot(TpRange_SatEnc_Vert, deltaAmp_SatEnc, 'm', 'linewidth', 2)
    PlotBoi2('Normalized Time Period','$\Delta$Amplitude of Vertical Orbit',20,'LaTex')
    
    
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













