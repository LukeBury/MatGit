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
systemNames = ['Earth_Moon.CR3BP',...
               'Jupiter_Europa.CR3BP',...
               'Saturn_Ganymede.CR3BP',...
               'Saturn_Enceladus.CR3BP',...
               'Uranus_Triton.CR3BP'];

for sys_i = 1:length(systemNames)
    sysName = systemNames(sys_i);
    
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
    [JC_L2] = JacobiConstantCalculator(prms.u, rLPs_n(2,:), [0,0,0]);
    
    %%% Calculat JC on surface
    JC_spx = getJacobiConstant_ZH([1-prms.u+prms.R2, 0, 0, 0, 0, 0], prms);
    JC_smx = getJacobiConstant_ZH([1-prms.u-prms.R2, 0, 0, 0, 0, 0], prms);
    JC_spy = getJacobiConstant_ZH([1-prms.u, prms.R2, 0, 0, 0, 0], prms);
    JC_smy = getJacobiConstant_ZH([1-prms.u, -prms.R2, 0, 0, 0, 0], prms);
    
    
    
    warning('finish this and also make a plot in plotShallowImpactAngles.m that shows the manifolds propagated backwards')
% % % %     %%% Calculate Jacobi constant of L2
% % % %     [JC_L2] = JacobiConstantCalculator(mu, L2_BCR, [0,0,0]);
% % % % 
% % % %     %%% Calculate difference between sc JC and JC_L2
% % % %     dJC_L2 = JC_L2 - JC_sc;
% % % % 
% % % %     %%% Convert the difference in JC to difference in velocity
% % % %     L2FlyoverVelocity_kps = sqrt(dJC_L2)*vNorm;
% % % %     L2FlyoverVelocity_mps = L2FlyoverVelocity_kps * 1000;

% % % %     %%% Calculate Jacobi constant of positive-x surface
% % % %     x_posXSurface = 1 - prms.u + prms.R2;
% % % %     [JC_posXSurface] = getJacobiConstant_ZH([x_posXSurface, 0, 0, 0, 0, 0], prms);
% % % % 
% % % %     %%% Calculate difference between sc JC and JC_L2
% % % %     dJC = JC_sc - JC_posXSurface;
% % % % 
% % % %     %%% Make sure the surface is accessible
% % % %     if dJC > 0
% % % %         warning('Surface not reachable')
% % % %         landingVelocity_mps = NaN;
% % % %         return
% % % %     end
% % % % 
% % % %     %%% Convert the difference in JC to difference in velocity
% % % %     landingVelocity_mps = sqrt(abs(dJC)) * vNorm * 1000;
end



% -------------------------------------------------
%%% System
% -------------------------------------------------












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
fprintf('Elapsed time: %1.4f seconds\n',tocWhole)
















