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
%%% 3B system
primary   = bodies.jupiter;
secondary = bodies.europa;

%%% Normalizing constants
rNorm = secondary.a;         % n <-> km
tNorm = 1/secondary.meanMot; % n <-> sec
vNorm = rNorm / tNorm;       % n <-> km/sec

% -------------------------------------------------
% Work
% -------------------------------------------------
%%% Acquire Collinear Lagrange points
L123 = EquilibriumPoints(secondary.MR,1:3); % [3x3] of L1, L2, L3 normalized BCR coordinates
L1 = L123(1,:);
L2 = L123(2,:);

%%% Jacobi constant of Lagrange points
[JC_L1] = JacobiConstantCalculator(secondary.MR,L1,[0,0,0]);
[JC_L2] = JacobiConstantCalculator(secondary.MR,L2,[0,0,0]);
   
nAngles = 16;
r0s = zeros(nAngles,3);
JCs = zeros(nAngles,1);

angles = 0:(2*pi/nAngles):(2*pi-2*pi/nAngles);

for kk = 1:nAngles
    r0s(kk,:) = R3([-1,0,0],angles(kk)) .* secondary.R_n + [1-secondary.MR, 0, 0] ;
    
    [JC_ang] = JacobiConstantCalculator(secondary.MR,r0s(kk,:),[0,0,0]);
    JCs(kk) = JC_ang;
    
end

figure(1); hold all
plot(angles, JCs,'x','color',colors.std.black)
% Plo 

fprintf('Easiest places to land are +-y\n Hardest is -x, closely followed by +x\n')


