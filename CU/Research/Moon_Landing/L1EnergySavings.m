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
bodies = getBodyData();

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


for sys_i = systems
    sys = sys_i{1};
    
    %%% Acquire Collinear Lagrange points
    L123 = EquilibriumPoints(sys.secondary.MR,1:3); % [3x3] of L1, L2, L3 normalized BCR coordinates
    L1 = L123(1,:);
    L2 = L123(2,:);
    
    %%% Jacobi constant of Lagrange points
    [JC_L1] = JacobiConstantCalculator(sys.secondary.MR,L1,[0,0,0]);
    [JC_L2] = JacobiConstantCalculator(sys.secondary.MR,L2,[0,0,0]);
    
    %%% Jacobi constant of surface points
    [JC_sPX] = JacobiConstantCalculator(sys.secondary.MR,[1-sys.secondary.MR+sys.secondary.R_n,0,0],[0,0,0]);
    [JC_sNX] = JacobiConstantCalculator(sys.secondary.MR,[1-sys.secondary.MR-sys.secondary.R_n,0,0],[0,0,0]);
    
    %%% Caculating differences in JC
    d_L1sPX = abs(JC_L1 - JC_sPX);
    d_L1sNX = abs(JC_L1 - JC_sNX);
    d_L2sPX = abs(JC_L2 - JC_sPX);
    d_L2sNX = abs(JC_L2 - JC_sNX);
    
    %%% Calculating minimum necessary DV to stop on surface
    sys.dV_L1sPX_mps = sqrt(d_L1sPX)*systems{kk}.vNorm*1000;
    sys.dV_L1sNX_mps = sqrt(d_L1sNX)*systems{kk}.vNorm*1000;
    sys.dV_L2sPX_mps = sqrt(d_L2sPX)*systems{kk}.vNorm*1000;
    sys.dV_L2sNX_mps = sqrt(d_L2sNX)*systems{kk}.vNorm*1000;
    
    %%% Calculating differences between L1 and L2
    sys.diff = sys.dV_L2sPX_mps - sys.dV_L1sPX_mps;
    
    %%% Printing
    fprintf('------------------------------\n')
    fprintf('%s - %s\n',sys.primary.name,sys.secondary.name)
    fprintf('dV - L1 to sPX: %2.2f (m/s)\n',sys.dV_L1sPX_mps)
    fprintf('dV - L1 to sNX: %2.3f (m/s)\n',sys.dV_L1sNX_mps)
    fprintf('dV - L2 to sPX: %2.3f (m/s)\n',sys.dV_L2sPX_mps)
    fprintf('dV - L2 to sNX: %2.3f (m/s)\n\n',sys.dV_L2sNX_mps)
    fprintf('L1 savings: %2.3f (m/s)\n\n',sys.diff)
    
    figure(1); hold all
    plot(sys.secondary.MR, sys.dV_L1sPX_mps,'bx')
    plot(sys.secondary.MR, sys.dV_L2sPX_mps,'bx')
    
    figure(2); hold all
    plot(sys.secondary.MR, sys.diff,'bx')
    
    
end

figure(1)
set(gca,'xscale','log')
PlotBoi2('System Mass Ratio', 'Minimum $\Delta$V for Lander (m/s)', 14, 'LaTex')

figure(2)
set(gca,'xscale','log')
PlotBoi2('System Mass Ratio', 'L$_1$ $\Delta$V Savings (m/s)', 14, 'LaTex')











