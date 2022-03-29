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
    
    %%% Getting normalized mean motion
%     tNorm = sqrt((rNorm^3)/(bodies.constants.G*(primary.mass + secondary.mass)));
%     prms.n = secondary.meanMot * tNorm;
    
    %%% Setting parameters structure
    prms.u = secondary.MR;
    prms.R1 = primary.R / rNorm;
    prms.R2 = secondary.R_n;
    prms.J2p = primary.J2;
%     prms.J4p = primary.J4;
%     prms.J6p = primary.J6;
    prms.J2s = secondary.J2;


%     prms.u = 0.1899309048e-6;
%     prms.R1 = 60268/238042;
%     prms.R2 = 252.1/238042;
%     prms.J2p = 1.6298e-2;
%     prms.J2s = 2.5e-3;
    prms.n = sqrt(1 + (3/2)*((1-prms.u)*(prms.R1^2)*prms.J2p + prms.u*(prms.R2^2)*prms.J2s));
    
    %%% Equillibrium Points
    rLPs_n_van = EquilibriumPoints(prms.u, 1, 1:3);
    rLPs_n_ZH  = collinearEquilibriumPoints_ZH(prms);
    
    dL123_n = rLPs_n_ZH(:,1) - rLPs_n_van(:,1);
    dL123 = dL123_n .* rNorm;
    
    fprintf('%s\n',sysName_k)
    fprintf('dL1n: %1.3e,\tdL2n: %1.3e,\tdL3n: %1.3e\n', dL123_n(1), dL123_n(2), dL123_n(3))
    fprintf('dL1:  %1.3f,\t\tdL2:  %1.3f,\t\tdL3:  %1.3f km\n', dL123(1), dL123(2), dL123(3))
    fprintf('%1.3f      & %1.3f      & %1.3f \n', dL123(1), dL123(2), dL123(3))
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
    
    %%% Getting normalized mean motion
    tNorm = sqrt((rNorm^3)/(bodies.constants.G*(primary.mass + secondary.mass)));
    prms.n = secondary.meanMot * tNorm;
    
    %%% Setting parameters structure
    prms.u = secondary.MR;
    prms.R1 = primary.R / rNorm;
    prms.R2 = secondary.R_n;
    prms.J2p = primary.J2;
    prms.J4p = primary.J4;
    prms.J6p = primary.J6;
    prms.J2s = secondary.J2;
    
    %%% Equillibrium Points
    rLPs_n_van = EquilibriumPoints(prms.u, 1, 1:3);
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




















%%% Set primary & secondary
[primary, secondary] = assignPrimaryAndSecondary_CR3BP('Saturn_Enceladus', bodies);
% 
% prms.u = secondary.MR;
% prms.R1 = primary.R;
% prms.R2 = 252.1;
% secondary.a = 238042;
% rNorm = secondary.a;
% prms.J2p = 1.6298e-2;
% prms.J2s = 2.5e-3;

prms.u = 0.19e-6;
secondary.a = 238042;
rNorm = secondary.a;
prms.R1 = 60268/rNorm;
prms.R2 = 252.1/rNorm;
prms.J2p = 1.6298e-2;
prms.J2s = 2.5e-3;
% prms.n = sqrt(1 + (3/2)*((1-prms.u)*(prms.R1^2)*prms.J2p + prms.u*(prms.R2^2)*prms.J2s));
% prms.n = sqrt(1 + (3/2)*((prms.R1^2)*prms.J2p + (prms.R2^2)*prms.J2s));
prms.n = sqrt(1 + (3/2)*((prms.R1^2)*prms.J2p/(secondary.a^2) + (prms.R2^2)*prms.J2s/(secondary.a^2)));

%%% Equillibrium Points
rLPs_n_van = EquilibriumPoints(prms.u, 1, 1:3);
rLPs_n_ZH  = collinearEquilibriumPoints_ZH(prms);

dL123_n = rLPs_n_ZH(:,1) - rLPs_n_van(:,1);
dL123 = dL123_n .* rNorm;

% fprintf('%s\n',sysName_k)
fprintf('dL1n: %1.3e,\tdL2n: %1.3e,\tdL3n: %1.3e\n', dL123_n(1), dL123_n(2), dL123_n(3))
fprintf('dL1:  %1.3f,\t\tdL2:  %1.3f,\t\tdL3:  %1.3f km\n', dL123(1), dL123(2), dL123(3))
% fprintf('%1.3f      & %1.3f      & %1.3f \n', dL123(1), dL123(2), dL123(3))



% page 13

