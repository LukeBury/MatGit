% ========================================================================
%%% Description
% ========================================================================
% For determining the minimum possible landing speeds on the surface of
% various moons

% Created: 11/14/19
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
%%% Systems
% -------------------------------------------------
systemNames = {'Earth_Moon.CR3BP',...
               'Jupiter_Europa.CR3BP',...
               'Jupiter_Ganymede.CR3BP',...
               'Saturn_Enceladus.CR3BP',...
               'Neptune_Triton.CR3BP'};
% systemNames = {'Jupiter_Europa.CR3BP'};
           
for sys_i = 1:length(systemNames)
    sysName = systemNames{sys_i};
    
    % --------------------------
    % Set primary & secondary
    % --------------------------
    [primary, secondary] = assignPrimaryAndSecondary_CR3BP(sysName, bodies);
    
    % --------------------------
    % Setup
    % --------------------------
    %%% Normalizing constants
    rNorm = secondary.a;         % n <-> km
    tNorm = 1/secondary.meanMot; % n <-> sec
    vNorm = rNorm / tNorm;       % n <-> km/sec
    
    %%% Setting parameter structure
    prms.u  = secondary.MR;
    prms.R1 = primary.R / rNorm;
    prms.R2 = secondary.R_n;
    
    %%% Equillibrium Points
    rLPs_n = collinearEquilibriumPoints_ZH(prms);
    
    % --------------------------
    % Calculations
    % --------------------------
    %%% Calculate Jacobi constant of L2
    if isequal(sysName, 'Earth_Moon.CR3BP')
        LP = 1;
    else
        LP = 2;
    end
    JC_LP = getJacobiConstant_ZH([rLPs_n(LP,:), 0, 0, 0], prms);
    
    %%% Calculat JC on surface
%     JC_spx = getJacobiConstant_ZH([1-prms.u+prms.R2, 0, 0, 0, 0, 0], prms);
%     JC_smx = getJacobiConstant_ZH([1-prms.u-prms.R2, 0, 0, 0, 0, 0], prms);
    JC_spy = getJacobiConstant_ZH([1-prms.u, prms.R2, 0, 0, 0, 0], prms);
%     JC_smy = getJacobiConstant_ZH([1-prms.u, -prms.R2, 0, 0, 0, 0], prms);
    
%     velocityDiff_spx_mps = sqrt(JC_spx - JC_L2) * vNorm * 1000
%     velocityDiff_smx_mps = sqrt(JC_smx - JC_L2) * vNorm * 1000
    velocityDiff_spy_mps = sqrt(JC_spy - JC_LP) * vNorm * 1000;
%     velocityDiff_smy_mps = sqrt(JC_smy - JC_L2) * vNorm * 1000

    fprintf('============================================\n')
    fprintf('%s\n',sysName)
    fprintf('Mass Ratio: %1.1e\n', prms.u)
    fprintf('Secondary Radius: %1.0f\n', secondary.R)
    fprintf('Minimum impact velocity from L%1d: %1.0f mps\n', LP, velocityDiff_spy_mps)
    
end
fprintf('============================================\n')

% ========================================================================
%%% Closeout
% ========================================================================
tocWhole = toc(ticWhole);
fprintf('Elapsed time: %1.4f seconds\n',tocWhole)















