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

% ========================================================================
%%% Analysis of collinear equilibrium points
% ========================================================================
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
%%% Perturbed Method - Loop through systems and determine movement of L1, L2, and L3
% -------------------------------------------------
fprintf('\nPerturbed Method accounting for J2p and J2s in mean motion, but  with no perturbations in model\n')

for kk = 1:length(sysName)
    clear prms
    sysName_k = sysName{kk};
    
    %%% Set primary & secondary
    [primary, secondary] = assignPrimaryAndSecondary_CR3BP(sysName_k, bodies);
    
    %%% Factor for normalizing distances
    rNorm = secondary.a; % n <-> km
    
    %%% Getting normalized mean motion
%     tNorm = sqrt((rNorm^3)/(bodies.constants.G*(primary.mass + secondary.mass)));
%     prms.n = secondary.meanMot * tNorm;
    
    %%% Setting parameters structure
    prms.u = secondary.MR;
%     prms.R1 = primary.R / rNorm;
%     prms.R2 = secondary.R_n;
%     prms.J2p = primary.J2;
% %     prms.J4p = primary.J4;
% %     prms.J6p = primary.J6;
%     prms.J2s = secondary.J2;


%     prms.n = sqrt(1 + (3/2)*((1-prms.u)*(prms.R1^2)*prms.J2p + prms.u*(prms.R2^2)*prms.J2s));
    prms.n = sqrt(1 + (3/2)*((1-secondary.MR)*((primary.R/rNorm)^2)*primary.J2 + secondary.MR*(secondary.R_n^2)*secondary.J2));
    
    %%% Equillibrium Points
    rLPs_n_van = EquilibriumPoints(prms.u, 1, 1:3);
    rLPs_n_ZH  = collinearEquilibriumPoints_ZH(prms);
    
    dL123_n = rLPs_n_ZH(:,1) - rLPs_n_van(:,1);
    dL123 = dL123_n .* rNorm;
    
%     fprintf('%s\n',sysName_k)
%     fprintf('dL1n: %1.3e,\tdL2n: %1.3e,\tdL3n: %1.3e\n', dL123_n(1), dL123_n(2), dL123_n(3))
%     fprintf('dL1:  %1.3f,\t\tdL2:  %1.3f,\t\tdL3:  %1.3f km\n', dL123(1), dL123(2), dL123(3))
    fprintf('%s\t & %s\t & %1.9f\t & %1.4f\t & %1.4f\t & %1.4f \\\\\n',primary.title, secondary.title, prms.n, dL123(1), dL123(2), dL123(3))
%     fprintf('\n============================================================\n')

end % kk = 1:length(sysName)


% -------------------------------------------------
%%% Universal Method - Loop through systems and determine movement of L1, L2, and L3
% -------------------------------------------------
fprintf('\n\nUniversal Method with no perturbations in model\n')
for kk = 1:length(sysName)
    clear prms
    sysName_k = sysName{kk};
    
    %%% Set primary & secondary
    [primary, secondary] = assignPrimaryAndSecondary_CR3BP(sysName_k, bodies);
    
    %%% Factor for normalizing distances
    rNorm = secondary.a; % n <-> km
    
    %%% Getting normalized mean motion
    tNorm = sqrt((rNorm^3)/(bodies.constants.G*(primary.mass + secondary.mass)));
    prms.n = secondary.meanMot * tNorm;
    
    %%% Setting parameters structure
    prms.u = secondary.MR;
%     prms.R1 = primary.R / rNorm;
%     prms.R2 = secondary.R_n;
%     prms.J2p = primary.J2;
% %     prms.J4p = primary.J4;
% %     prms.J6p = primary.J6;
%     prms.J2s = secondary.J2;


%     prms.u = 0.1899309048e-6;
%     prms.R1 = 60268/238042;
%     prms.R2 = 252.1/238042;
%     prms.J2p = 1.6298e-2;
%     prms.J2s = 2.5e-3;
%     prms.n = sqrt(1 + (3/2)*((1-prms.u)*(prms.R1^2)*prms.J2p + prms.u*(prms.R2^2)*prms.J2s));
    
    %%% Equillibrium Points
    rLPs_n_van = EquilibriumPoints(prms.u, 1, 1:3);
    rLPs_n_ZH  = collinearEquilibriumPoints_ZH(prms);
    
    dL123_n = rLPs_n_ZH(:,1) - rLPs_n_van(:,1);
    dL123 = dL123_n .* rNorm;
    
%     fprintf('%s\n',sysName_k)
%     fprintf('dL1n: %1.3e,\tdL2n: %1.3e,\tdL3n: %1.3e\n', dL123_n(1), dL123_n(2), dL123_n(3))
%     fprintf('dL1:  %1.3f,\t\tdL2:  %1.3f,\t\tdL3:  %1.3f km\n', dL123(1), dL123(2), dL123(3))
    fprintf('%s\t & %s\t & %1.9f\t & %1.4f\t & %1.4f\t & %1.4f \\\\\n',primary.title, secondary.title, prms.n, dL123(1), dL123(2), dL123(3))
%     fprintf('\n============================================================\n')

end % kk = 1:length(sysName)


% -------------------------------------------------
%%% Universal Method - with perturbations
% -------------------------------------------------
fprintf('\n\nUniversal Method with J2p and J2s\n')
for kk = 1:length(sysName)
    clear prms
    sysName_k = sysName{kk};
    
    %%% Set primary & secondary
    [primary, secondary] = assignPrimaryAndSecondary_CR3BP(sysName_k, bodies);
    
    %%% Factor for normalizing distances
    rNorm = secondary.a; % n <-> km
    
    %%% Getting normalized mean motion
    tNorm = sqrt((rNorm^3)/(bodies.constants.G*(primary.mass + secondary.mass)));
    prms.n = secondary.meanMot * tNorm;
    
    %%% Setting parameters structure
    prms.u = secondary.MR;
    prms.R1 = primary.R / rNorm;
    prms.R2 = secondary.R_n;
    prms.J2p = primary.J2;
% %     prms.J4p = primary.J4;
% %     prms.J6p = primary.J6;
    prms.J2s = secondary.J2;


%     prms.u = 0.1899309048e-6;
%     prms.R1 = 60268/238042;
%     prms.R2 = 252.1/238042;
%     prms.J2p = 1.6298e-2;
%     prms.J2s = 2.5e-3;
%     prms.n = sqrt(1 + (3/2)*((1-prms.u)*(prms.R1^2)*prms.J2p + prms.u*(prms.R2^2)*prms.J2s));
    
    %%% Equillibrium Points
    rLPs_n_van = EquilibriumPoints(prms.u, 1, 1:3);
    rLPs_n_ZH  = collinearEquilibriumPoints_ZH(prms);
    
    dL123_n = rLPs_n_ZH(:,1) - rLPs_n_van(:,1);
    dL123 = dL123_n .* rNorm;
    
%     fprintf('%s\n',sysName_k)
%     fprintf('dL1n: %1.3e,\tdL2n: %1.3e,\tdL3n: %1.3e\n', dL123_n(1), dL123_n(2), dL123_n(3))
%     fprintf('dL1:  %1.3f,\t\tdL2:  %1.3f,\t\tdL3:  %1.3f km\n', dL123(1), dL123(2), dL123(3))
    fprintf('%s\t & %s\t & %1.9f\t & %1.4f\t & %1.4f\t & %1.4f \\\\\n',primary.title, secondary.title, prms.n, dL123(1), dL123(2), dL123(3))
%     fprintf('\n============================================================\n')

end % kk = 1:length(sysName)




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

















