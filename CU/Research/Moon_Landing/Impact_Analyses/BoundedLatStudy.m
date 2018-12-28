%%% INFORMATION
% When grid searches are performed over intial conditions in the L1 or L2
% necks at low energies, we find that ballistic impact trajectories only 
% reach certain latitudes until a threshold energy is reached which allow
% impacts of all latitudes. Curiously, for sufficiently low energies, the
% bounds set on reachable latitudes don't come from the energy domain -
% these trajectories have a sufficient jacobi constant to reach any point
% on the surface. So what part of phase space describes this boundary? In
% this script, I plan to investigate angular velocity as a boundary, and
% will also study the maximum latitude reached before impact.
% ========================================================================
clear
clc
close all
mbinPath = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/mbin';
moonFuncsPath = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/CU/Research/Moon_Landing/Moon_Landing_funcs';
addpath(genpath(mbinPath))
addpath(genpath(moonFuncsPath))
tic
 
% ========================================================================
%%% Run Switches
% ========================================================================
plot_trajectories = 1;

% ========================================================================
% Data files
% ========================================================================
%%% File location
MatlabOutputsPath = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/';

%%% Low-angle landing file
lowImpactAngleFile = [MatlabOutputsPath,'F.iGS_eurL2_100mps_50km_149v0s_land.txt'];
 

% ========================================================================
%%% Importing Data and Setting up System
% ========================================================================
%%% General data on solar system bodies
bodies = getBodyData(mbinPath);

%%% Color options/schemes
colors = get_colors();

%%% Bodies
primary = bodies.jupiter;
secondary = bodies.europa;

%%% Colinear Lagrange Points
L123 = EquilibriumPoints(secondary.MR,1:3); % [3x3] of L1, L2, L3 normalized BCR coordinates

%%% Normalizing constants
rNorm = secondary.a;         % n <-> km
tNorm = 1/secondary.meanMot; % n <-> sec
vNorm = rNorm / tNorm;       % n <-> km/sec


% ========================================================================
%%% Studying Angular Momentum
% ========================================================================
% ---------------------------------------
% Grabbing initial conditions
% ---------------------------------------
%%% Load low impact data
lowImpactMat = dlmread(lowImpactAngleFile,',',1,0);
nLowImpactTraj = size(lowImpactMat,1);

% ---------------------------------------
% Setting up for integration
% ---------------------------------------
%%% Time
t_i = 0; % sec
t_f = 2*pi;
n_dt = 1000;
time0_n = linspace(t_i,t_f,n_dt);

%%% Choosing ode45 tolerance
tol = 2.22045e-14;

%%% Setting integrator options
options = odeset('RelTol',tol,'AbsTol',tol);
options_ImpactEscape = odeset('Events',@event_ImpactEscape_CR3Bn,'RelTol',tol,'AbsTol',tol);

%%% Setting extra parameters
prms.u = secondary.MR;
prms.R2_n = secondary.R_n;
prms.L1x = L123(1,1);
prms.L2x = L123(2,1);


%%% How many trajectories to loop through?
nTraj = size(lowImpactMat,1);
% nTraj = 50;

%%% Preallocating for storage of all trajectories
trajs(nTraj).azEl = [];
trajs(nTraj).rBCR = [];


%%% Looping through trajectories
for ii = 1:nTraj
    %%% Clearing / Preallocating
    X_BCR_n = [];
    
    %%% Initial conditions
    X0_n = lowImpactMat(ii,2:7)';

    % ---------------------------------------
    % Propagating trajectory
    % ---------------------------------------
    %%% Propagating trajectory
    [time_n, X_BCR_n, time_eventImpact, X_eventImpact, index_eventImpact] = ode113(@Int_CR3Bn,...
                time0_n, X0_n, options_ImpactEscape, prms);
    
    if plot_trajectories == 1
        %%% Storing trajectory
        trajs(ii).rBCR = X_BCR_n(:,1:3);
    end
    
    % ---------------------------------------
    % Angular momentum
    % ---------------------------------------
    trajs(ii).azEl = zeros(size(X_BCR_n,1),2);
%     angMoms(1).azEl = zeros(size(X_BCR_n,1),2);

    for kk = 1:size(X_BCR_n,1)
        %%% Determine azimuth and elevation at current time
        [Az_rad,El_rad,Rad] = cart2sph(X_BCR_n(kk,4), X_BCR_n(kk,5), X_BCR_n(kk,6));

        %%% Store azimuth and elevation
        trajs(ii).azEl(kk,:) = [Az_rad, El_rad];
%         angMoms(1).azEl_rad(kk,:) = [Az_rad, El_rad];
    end

end % ii = 1:size(lowImpactMat,1)

%%% Rearranging data
AzEls_rad = vertcat(trajs.azEl);
if plot_trajectories == 1
    rBCRs = vertcat(trajs.rBCR);
end

% ---------------------------------------
% Plotting angular momentum
% ---------------------------------------
figure(22); hold all
plot(AzEls_rad(:,1),AzEls_rad(:,2),'b.')
PlotBoi2('Azimuth, $rad$','Elevation, $rad$',14,'LaTex')
xlim([-pi pi])
ylim([-pi/2 pi/2])


% ---------------------------------------
% Plotting system
% ---------------------------------------
if plot_trajectories == 1
    figure; hold all
    plot3(rBCRs(:,1),rBCRs(:,2),rBCRs(:,3),'b','linewidth',1.5)
    plotBodyTexture3(secondary.R_n,[1-secondary.MR,0,0],secondary.img)
    PlotBoi3('$x_n$','$y_n$','$z_n$',16,'LaTex')
    axis equal
end




toc




