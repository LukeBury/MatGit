clear
clc
close all
mbinPath = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/mbin';
addpath(genpath(mbinPath))

% ========================================================================
%%% Importing Data
% ========================================================================
%%% General data on solar system bodies
bodies = getBodyData(mbinPath);

%%% Color options/schemes
colors = get_colors();

% ========================================================================
%%% Run Switches
% ========================================================================


% ========================================================================
%%% Load Trajectories
% ========================================================================
% ------------------------------------
%%% System
% ------------------------------------
primary   = bodies.jupiter;
secondary = bodies.europa;

% ------------------------------------
%%% Initial conditions and TOF
% ------------------------------------
%%% Low-angle-landing-file column specifiers
low_X0_c = 2:7;
low_lat_c = 9;
low_lon_c = 10;
low_tf_c  = 11;

%%% Choose data file
% myFile = dlmread('/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_eurL2_50mps_50km_149v0s_land.txt',',',1,0);
% myFile = dlmread('/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_eurL2_100mps_50km_149v0s_land.txt',',',1,0);
% myFile = dlmread('/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_eurL2_150mps_50km_149v0s_land.txt',',',1,0);
% myFile = dlmread('/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_eurL2_200mps_50km_149v0s_land.txt',',',1,0);
% myFile = dlmread('/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_eurL2_250mps_50km_149v0s_land.txt',',',1,0);
% myFile = dlmread('/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_eurL2_300mps_50km_149v0s_land.txt',',',1,0);
myFile = dlmread('/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/F.iGS_eurL2_350mps_50km_149v0s_land.txt',',',1,0);

%%% Choose trajectory
myIdx = find(abs(myFile(:,low_lat_c)) == max(abs(myFile(:,low_lat_c))));
%%% For Europa 350, these three are all high-latitude with short TOF
% myIdx = 64;
% myIdx = 269;
% myIdx = 782;

myTraj = myFile(myIdx,:);

X0_n = myTraj(low_X0_c)';
% b = sortrows(myFile,low_lat_c);
% b(1:5,:)
% X0_n = b(5,low_X0_c)';

%%% Look at landing time vs impact latitude of all trajectories
% % figure
% % plot(myFile(:,low_tf_c),myFile(:,low_lat_c),'.')
% % PlotBoi2('tf','Impact Latitude, deg',18,'LaTex')


% ------------------------------------
%%% Normalizing factors and equillibrium points
% ------------------------------------
%%% Normalizing factors
rNorm = secondary.a;         % n <-> km
tNorm = 1/secondary.meanMot; % n <-> sec
vNorm = rNorm / tNorm;       % n <-> km/sec

%%% Collinear equilibrium points
L123 = EquilibriumPoints(secondary.MR,1:3); % [3x3] of L1, L2, L3 normalized BCR coordinates

% ========================================================================
%%% Propagate and Plot
% ========================================================================
% ------------------------------------
%%% Time and integration options
% ------------------------------------
%%% Create initial normalized time vector
t0 = 0;
n_t = 10000;
tf = 4*pi;
time0_n = linspace(t0,tf,n_t);

%%% Choosing ode tolerance
tol = 1e-13;

%%% Setting integrator options
options = odeset('Event',@event_Impact_CR3Bn, 'RelTol',tol,'AbsTol',tol);

%%% Setting necessary parameters for integration
prms.u    = secondary.MR;
prms.R2_n = secondary.R_n;

% ------------------------------------
%%% Integrate
% ------------------------------------
[time_n, X_BCR_n] = ode113(@Int_CR3Bn, time0_n, X0_n, options, prms);

% ------------------------------------
%%% Plot
% ------------------------------------
figure
figure(1); hold all
plot3(X_BCR_n(:,1),X_BCR_n(:,2),X_BCR_n(:,3),'r','linewidth',1.5)
[JC_scInitial] = JacobiConstantCalculator(secondary.MR,X0_n(1:3)',X0_n(4:6)');
plotCR3BP_YZNeck( JC_scInitial, secondary.MR , 2, 0, prms, colors.std.black, 1.5)
plotCR3BP_YZNeck( JC_scInitial, secondary.MR , 1, 0, prms, colors.std.black, 1.5)
plotCR3BP_Neck(secondary,L123,JC_scInitial,600,200,colors.std.black,1.5)
plotBodyTexture3(secondary.R_n, [1-secondary.MR, 0, 0], secondary.img)
view(0,90)
axis equal
PlotBoi3('$X_n$','$Y_n$','$Z_n$',18,'LaTex')

% ------------------------------------
%%% Post processing
% ------------------------------------
[lat_deg, lon_deg] = BCR2latlon(X_BCR_n(end,1:3), 'secondary', secondary.MR);
fprintf('Impact Latitude:   %2.3f°\n',lat_deg)
fprintf('Impact Longitude:  %2.3f°\n',lon_deg)




























