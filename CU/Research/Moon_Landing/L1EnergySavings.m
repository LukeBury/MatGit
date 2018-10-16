clear
clc
close all
mbinPath = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/mbin';
addpath(genpath(mbinPath))
tic

% ========================================================================
%%% Run/Plot Switches
% ========================================================================

% ========================================================================
%%% Importing Data
% ========================================================================
%%% General data on solar system bodies
bodies = getBodyData(mbinPath);

%%% Color options/schemes
colors = get_colors();

% ========================================================================
%%% Setup
% ========================================================================
% -------------------------------------------------
% Choose 3B systems 
% -------------------------------------------------
primaries   = {bodies.earth, bodies.jupiter, bodies.jupiter,  bodies.saturn,    bodies.saturn, bodies.neptune};
secondaries = {bodies.moon,  bodies.europa,  bodies.ganymede, bodies.enceladus, bodies.titan,  bodies.triton};

if length(primaries) ~= length(secondaries)
    warning('Incorrect system allignment')
    return
end

systems{length(primaries)}.primary = [];
for kk = 1:length(primaries)
    systems{kk}.primary = primaries{kk};
    systems{kk}.secondary = secondaries{kk};
    
    systems{kk}.rNorm = systems{kk}.secondary.a;         % n <-> km
    systems{kk}.tNorm = 1/systems{kk}.secondary.meanMot; % n <-> sec
    systems{kk}.vNorm = systems{kk}.rNorm / systems{kk}.tNorm;       % n <-> km/sec
    
end

dataMat = zeros(length(systems),2);
counter = 0;
for sys_i = systems
    counter = counter + 1;
    sys = sys_i{1};
    
    %%% Acquire Collinear Lagrange points
    L123 = EquilibriumPoints(sys.secondary.MR,1:3); % [3x3] of L1, L2, L3 normalized BCR coordinates
    L1 = L123(1,:);
    L2 = L123(2,:);
    
    %%% Jacobi constant of Lagrange points
    [JC_L1] = JacobiConstantCalculator(sys.secondary.MR,L1,[0,0,0]);
    [JC_L2] = JacobiConstantCalculator(sys.secondary.MR,L2,[0,0,0]);
    
    %%% Jacobi constant of surface points
%     [JC_sPX] = JacobiConstantCalculator(sys.secondary.MR,[1-sys.secondary.MR+sys.secondary.R_n,0,0],[0,0,0]);
    [JC_sNX] = JacobiConstantCalculator(sys.secondary.MR,[1-sys.secondary.MR-sys.secondary.R_n,0,0],[0,0,0]);
    [JC_sY] = JacobiConstantCalculator(sys.secondary.MR,[1-sys.secondary.MR,sys.secondary.R_n,0],[0,0,0]);
    
% % %     fprintf('-------------\n') % ** Difference is between 2-5 orders of magnitude
% % %     JC_L1 - JC_L2
% % %     JC_sNX - JC_sPX
% % %     JC_sY - JC_sPX
% % %     JC_sY - JC_sNY

%     %%% Caculating differences in JC
%     d_L1sNX = JC_sNX - JC_L1;
%     d_L1sY  = JC_sY  - JC_L1;
%     d_L2sNX = JC_sNX - JC_L2;
%     d_L2sY  = JC_sY  - JC_L2;
    
    %%% Calculating minimum necessary DV to stop on surface
    sys.dV_L1sNX_mps = sqrt(JC_sNX - JC_L1)*systems{kk}.vNorm*1000; % m/s
    sys.dV_L1sY_mps  = sqrt(JC_sY  - JC_L1)*systems{kk}.vNorm*1000; % m/s
    sys.dV_L2sNX_mps = sqrt(JC_sNX - JC_L2)*systems{kk}.vNorm*1000; % m/s
    sys.dV_L2sY_mps  = sqrt(JC_sY  - JC_L2)*systems{kk}.vNorm*1000; % m/s
    
    sys.dV_minLanding = sys.dV_L1sY_mps;
    
    %%% Calculating differences between L1 and L2
    sys.diff_sNX_mps = sys.dV_L2sNX_mps - sys.dV_L1sNX_mps;
    sys.diff_sY_mps  = sys.dV_L2sY_mps  - sys.dV_L1sY_mps;

    %%% Printing
    fprintf('------------------------------\n')
    fprintf('%s - %s\n',sys.primary.name,sys.secondary.name)
%     fprintf('dV - L1 to sNX: %2.3f (m/s)\n',sys.dV_L1sNX_mps)
%     fprintf('dV - L2 to sNX: %2.3f (m/s)\n\n',sys.dV_L2sNX_mps)
    fprintf('Minimum \x0394V to land:\t %2.1f (m/s)\n',sys.dV_minLanding)
    fprintf('L1 savings at (-x):\t %2.5f (m/s)\n',sys.diff_sNX_mps)
    fprintf('L1 savings at (+-y):\t %2.5f (m/s)\n\n',sys.diff_sY_mps)
    
    figure(1); hold all
    plot(sys.secondary.MR, sys.dV_L1sNX_mps,'r.','markersize',30)
    plot(sys.secondary.MR, sys.dV_L2sNX_mps,'b.','markersize',30)
    
    figure(2); hold all
    plot(sys.secondary.MR, sys.diff_sY_mps,'b.','markersize',30)
    
    %%% Storing values for matrix to be printed to LaTex
    dataMat(counter,1) = round(sys.dV_minLanding,2);
    dataMat(counter,2) = round(sys.diff_sY_mps,2);
end

figure(1)
set(gca,'xscale','log')
PlotBoi2('System Mass Ratio', 'Minimum $\Delta$V for Lander (m/s)', 14, 'LaTex')

figure(2)
set(gca,'xscale','log')
PlotBoi2('System Mass Ratio', 'L$_1$ $\Delta$V Savings at surface +-y (m/s)', 14, 'LaTex')


fprintf('================================================\n')
fprintf('Results correlated with MR, but also dependent on Rn\n')
fprintf('L1 will save you 15.6 m/s at the moon, but tops around 2 m/s for outer-moons\n')
fprintf('Low MR moons like Enceledus barely see savings\n')
fprintf('Easiest places to land are +-y, Hardest is -x, closely followed by +x\n\n')








