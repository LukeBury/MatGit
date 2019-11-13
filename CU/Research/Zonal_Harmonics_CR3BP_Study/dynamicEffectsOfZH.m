% ========================================================================
%%% Description
% ========================================================================
% A script for showing how selected zonal harmonic terms affect zero
% velocity curves and locations of collinear equilibrium points

% Created: 10/02/19
% Author : Luke Bury, luke.bury@colorado.edu
% ========================================================================
%%% Initialization
% ========================================================================
clear
clc
% close all
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
run_collinearEquilibrium_analysis = 1;
run_ZVCurve_analysis              = 0;

% ========================================================================
%%% Analysis of collinear equilibrium points
% ========================================================================
if run_collinearEquilibrium_analysis == 1
% -------------------------------------------------
%%% Creating cell array of familes
% -------------------------------------------------
sysName{1} = 'Earth_Moon';
sysName{2} = 'Mars_Phobos';
sysName{3} = 'Jupiter_Io';
sysName{4} = 'Jupiter_Europa';
sysName{5} = 'Jupiter_Ganymede';
sysName{6} = 'Jupiter_Callisto';
sysName{7} = 'Saturn_Enceladus';
sysName{8} = 'Saturn_Titan';
sysName{9} = 'Neptune_Triton';

% -------------------------------------------------
%%% Loop through systems and determine movement of L1, L2, and L3
% -------------------------------------------------
for kk = 1:length(sysName)
    clear prms
    sysName_k = sysName{kk};
    
    %%% Set primary & secondary
    [primary, secondary] = assignPrimaryAndSecondary_CR3BP(sysName_k, bodies);
    
    %%% Factor for normalizing distances
    rNorm = secondary.a; % n <-> km
    
    %%% Setting parameters structure
    prms.u = secondary.MR;
    prms.R1 = primary.R / rNorm;
    prms.R2 = secondary.R_n;
    prms.J2p = primary.J2;
%     prms.J4p = primary.J4;
%     prms.J6p = primary.J6;
    prms.J2s = secondary.J2;
    
    %%% Equillibrium Points
    rLPs_n_van = EquilibriumPoints(prms.u, 1:3);
    rLPs_n_ZH  = collinearEquilibriumPoints_ZH(prms);
    
    dL123_n = rLPs_n_ZH(:,1) - rLPs_n_van(:,1);
    dL123 = dL123_n .* rNorm;
    
    fprintf('%s\n',sysName_k)
%     fprintf('dL1n: %1.3e,\tdL2n: %1.3e,\tdL3n: %1.3e\n', dL123_n(1), dL123_n(2), dL123_n(3))
    fprintf('dL1:  %1.3f,\t\tdL2:  %1.3f,\t\tdL3:  %1.3f km\n', dL123(1), dL123(2), dL123(3))
    fprintf('\n============================================================\n')

end % kk = 1:length(sysName)

% -------------------------------------------------
%%% Do it again, but just for Jupiter_Europa and Saturn_Enceladus, and bring in J4p and J6p
% -------------------------------------------------
fprintf('============================================================\n')
fprintf('============================================================\n')
clear prms
for kk = 1:2
    if kk     == 1
        sysName_k = 'Jupiter_Europa';
    elseif kk == 2
        sysName_k = 'Saturn_Enceladus';
    end
    
    %%% Set primary & secondary
    [primary, secondary] = assignPrimaryAndSecondary_CR3BP(sysName_k, bodies);
    
    %%% Factor for normalizing distances
    rNorm = secondary.a; % n <-> km
    
    %%% Setting parameters structure
    prms.u = secondary.MR;
    prms.R1 = primary.R / rNorm;
    prms.R2 = secondary.R_n;
    prms.J2p = primary.J2;
    prms.J4p = primary.J4;
    prms.J6p = primary.J6;
    prms.J2s = secondary.J2;
    
    %%% Equillibrium Points
    rLPs_n_van = EquilibriumPoints(prms.u, 1:3);
    rLPs_n_ZH  = collinearEquilibriumPoints_ZH(prms);
    
    dL123_n = rLPs_n_ZH(:,1) - rLPs_n_van(:,1);
    dL123 = dL123_n .* rNorm;
    
    fprintf('%s\n',sysName_k)
    fprintf('dL1n: %1.3e,\tdL2n: %1.3e,\tdL3n: %1.3e\n', dL123_n(1), dL123_n(2), dL123_n(3))
    fprintf('dL1:  %1.3f,\t\tdL2:  %1.3f,\t\tdL3:  %1.3f km\n', dL123(1), dL123(2), dL123(3))
    fprintf('\n============================================================\n')

end % kk = 1:length(sysName)

end % run_collinearEquilibrium_analysis == 1






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




















% 
% clc
% 
% for kk = 1:length(sysName)
%     clear prms
%     sysName_k = sysName{kk};
%     
%     %%% Set primary & secondary
%     [primary, secondary] = assignPrimaryAndSecondary_CR3BP(sysName_k, bodies);
%     
% %     rNorm = secondary.a; % n <-> km
% %     tNorm = 1/secondary.meanMot; % n <-> sec
% 
% %     fprintf('%s\n',sysName_k)
% %     fprintf('dL1n: %1.3e,\tdL2n: %1.3e,\tdL3n: %1.3e\n', dL123_n(1), dL123_n(2), dL123_n(3))
%     fprintf('%8s\t& %1.1f\t& %1.4e\t& %1.0f\t& %1.8e\t& %1.8f\\\\\n',secondary.title, secondary.R, secondary.MR, secondary.a, secondary.meanMot, secondary.J2)
% %     fprintf('\n============================================================\n')
% 
% end % kk = 1:length(sysName)
% 
% 
% 





