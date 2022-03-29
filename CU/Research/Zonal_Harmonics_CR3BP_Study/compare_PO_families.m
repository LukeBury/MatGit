% ========================================================================
%%% Description
% ========================================================================
% 

% Created: 
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
run_JE_Lyap  = true;
run_JE_Vert  = true;
run_JE_SHalo = true;

run_SE_Lyap  = true;
run_SE_Vert  = true;
run_SE_SHalo = true;


plot_primary   = 0;
plot_secondary = 1;
plot_L1        = 0;
plot_L2        = 0;

plot_onlyOneOrbit    = 1;
    chosenOrbitIndex = 285;
plot_stabilityIndices = 1;

% ========================================================================
%%% Setup
% ========================================================================
% -------------------------------------------------
%%% Basic options
% -------------------------------------------------
%%% Run options
n_Solutions = 10;

%%% Multiple shooter options
n_Nodes = 12;
iterMax = 500;

%%% Plot options
lw = 2; % linewidth

%%% Integration options
tol = 1e-13;
options = odeset('RelTol',tol,'AbsTol',tol);

% -------------------------------------------------
%%% Set bounds of TPs for PO families
% -------------------------------------------------
targetVec_JE_Lyap  = linspace(3.08, 5.51, n_Solutions);
targetVec_JE_Vert  = linspace(3.19, 6.18, n_Solutions);
targetVec_JE_SHalo = linspace(1.81, 3.11, n_Solutions);
targetVec_SE_Lyap  = linspace(3.11, 4.29, n_Solutions);
targetVec_SE_Vert  = linspace(3.22, 6.27, n_Solutions);
targetVec_SE_SHalo = linspace(2.3, 3.09, n_Solutions);

% -------------------------------------------------
%%% Load data
% -------------------------------------------------
%%% Filenames
fam_JE_Lyap_cr3bp = 'Jupiter_Europa.CR3BP.L2_Lyapunov.txt';
fam_JE_Lyap_ZH    = 'Jupiter_Europa.CR3BP_J2pJ4pJ6pJ2s.L2_Lyapunov.txt';

fam_JE_Vert_cr3bp = 'Jupiter_Europa.CR3BP.L2_Vertical.txt';
fam_JE_Vert_ZH = 'Jupiter_Europa.CR3BP_J2pJ4pJ6pJ2s.L2_Vertical.txt';

fam_JE_SHalo_cr3bp = 'Jupiter_Europa.CR3BP.L2_SHalo.txt';
fam_JE_SHalo_ZH = 'Jupiter_Europa.CR3BP_J2pJ4pJ6pJ2s.L2_SHalo.txt';

fam_SE_Lyap_cr3bp = 'Saturn_Enceladus.CR3BP.L2_Lyapunov.txt';
fam_SE_Lyap_ZH = 'Saturn_Enceladus.CR3BP_J2pJ4pJ6pJ2s.L2_Lyapunov.txt';

fam_SE_Vert_cr3bp = 'Saturn_Enceladus.CR3BP.L2_Vertical.txt';
fam_SE_Vert_ZH = 'Saturn_Enceladus.CR3BP_J2pJ4pJ6pJ2s.L2_Vertical.txt';

fam_SE_SHalo_cr3bp = 'Saturn_Enceladus.CR3BP.L2_SHalo.txt';
fam_SE_SHalo_ZH = 'Saturn_Enceladus.CR3BP_J2pJ4pJ6pJ2s.L2_SHalo.txt';

%%% Path from mbin to data
dataPathFromMBin = '/Data/InitialConditions/PO_Families/';

%%% Filenames with paths
datafile_JE_Lyap_cr3bp = [mbinPath, dataPathFromMBin, fam_JE_Lyap_cr3bp];
datafile_JE_Lyap_ZH = [mbinPath, dataPathFromMBin, fam_JE_Lyap_ZH];

datafile_JE_Vert_cr3bp = [mbinPath, dataPathFromMBin, fam_JE_Vert_cr3bp];
datafile_JE_Vert_ZH = [mbinPath, dataPathFromMBin, fam_JE_Vert_ZH];

datafile_JE_SHalo_cr3bp = [mbinPath, dataPathFromMBin, fam_JE_SHalo_cr3bp];
datafile_JE_SHalo_ZH = [mbinPath, dataPathFromMBin, fam_JE_SHalo_ZH];

datafile_SE_Lyap_cr3bp = [mbinPath, dataPathFromMBin, fam_SE_Lyap_cr3bp];
datafile_SE_Lyap_ZH = [mbinPath, dataPathFromMBin, fam_SE_Lyap_ZH];

datafile_SE_Vert_cr3bp = [mbinPath, dataPathFromMBin, fam_SE_Vert_cr3bp];
datafile_SE_Vert_ZH = [mbinPath, dataPathFromMBin, fam_SE_Vert_ZH];

datafile_SE_SHalo_cr3bp = [mbinPath, dataPathFromMBin, fam_SE_SHalo_cr3bp];
datafile_SE_SHalo_ZH = [mbinPath, dataPathFromMBin, fam_SE_SHalo_ZH];


%%% Read data in matrices
if run_JE_Lyap
    FamData_JE_Lyap_cr3bp = dlmread(datafile_JE_Lyap_cr3bp,',',1,0);
    FamData_JE_Lyap_ZH = dlmread(datafile_JE_Lyap_ZH,',',1,0);
end
if run_JE_Vert
    FamData_JE_Vert_cr3bp = dlmread(datafile_JE_Vert_cr3bp,',',1,0);
    FamData_JE_Vert_ZH = dlmread(datafile_JE_Vert_ZH,',',1,0);
end
if run_JE_SHalo
    FamData_JE_SHalo_cr3bp = dlmread(datafile_JE_SHalo_cr3bp,',',1,0);
    FamData_JE_SHalo_ZH = dlmread(datafile_JE_SHalo_ZH,',',1,0);
end
if run_SE_Lyap
    FamData_SE_Lyap_cr3bp = dlmread(datafile_SE_Lyap_cr3bp,',',1,0);
    FamData_SE_Lyap_ZH = dlmread(datafile_SE_Lyap_ZH,',',1,0);
end
if run_SE_Vert
    FamData_SE_Vert_cr3bp = dlmread(datafile_SE_Vert_cr3bp,',',1,0);
    FamData_SE_Vert_ZH = dlmread(datafile_SE_Vert_ZH,',',1,0);
end
if run_SE_SHalo
    FamData_SE_SHalo_cr3bp = dlmread(datafile_SE_SHalo_cr3bp,',',1,0);
    FamData_SE_SHalo_ZH = dlmread(datafile_SE_SHalo_ZH,',',1,0);
end
% -------------------------------------------------
%%% Set up the system
% -------------------------------------------------
% --------------------------
% Set primary & secondary
% --------------------------
[jupiter, europa] = assignPrimaryAndSecondary_CR3BP('Jupiter_Europa', bodies);
[saturn, enceladus] = assignPrimaryAndSecondary_CR3BP('Saturn_Enceladus', bodies);

%%% Normalizing constants
rNorm_JE = europa.a;         % n <-> km
rNorm_SE = enceladus.a;         % n <-> km

tNorm_JE_ZH = sqrt((rNorm_JE^3)/(bodies.constants.G*(jupiter.mass + europa.mass)));
tNorm_SE_ZH = sqrt((rNorm_JE^3)/(bodies.constants.G*(saturn.mass + enceladus.mass)));


%%% Create prms structures
prms_JE_cr3bp.n  = 1;
prms_JE_cr3bp.u  = europa.MR;
prms_JE_cr3bp.R1 = jupiter.R / rNorm_JE;
prms_JE_cr3bp.R2 = europa.R_n;

prms_SE_cr3bp.n  = 1;
prms_SE_cr3bp.u  = enceladus.MR;
prms_SE_cr3bp.R1 = saturn.R / rNorm_SE;
prms_SE_cr3bp.R2 = enceladus.R_n;

prms_JE_ZH     = prms_JE_cr3bp;
prms_JE_ZH.n   = europa.meanMot * tNorm_JE_ZH;
prms_JE_ZH.J2p = jupiter.J2;
prms_JE_ZH.J4p = jupiter.J4;
prms_JE_ZH.J6p = jupiter.J6;
prms_JE_ZH.J2s = europa.J2;

prms_SE_ZH     = prms_SE_cr3bp;
prms_SE_ZH.n   = enceladus.meanMot * tNorm_SE_ZH;
prms_SE_ZH.J2p = saturn.J2;
prms_SE_ZH.J4p = saturn.J4;
prms_SE_ZH.J6p = saturn.J6;
prms_SE_ZH.J2s = enceladus.J2;


%%% Equillibrium Points
rLPs_JE_cr3bp_n = collinearEquilibriumPoints_ZH(prms_JE_cr3bp);
rLPs_SE_cr3bp_n = collinearEquilibriumPoints_ZH(prms_SE_cr3bp);

rLPs_JE_ZH_n = collinearEquilibriumPoints_ZH(prms_JE_ZH);
rLPs_SE_ZH_n = collinearEquilibriumPoints_ZH(prms_SE_ZH);













% ========================================================================
%%% Closeout
% ========================================================================
tocWhole = toc(ticWhole);
fprintf('Elapsed time: %1.4f seconds\n',tocWhole)






















