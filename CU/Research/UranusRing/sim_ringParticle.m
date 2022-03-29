% ========================================================================
%%% Description
% ========================================================================
% 

% Created: 07/30/19
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
runManifolds = 1;
    manifoldsFromMonodromy = 1;
    manifoldsFromPert      = 0;

% ========================================================================
%%% Setup
% ========================================================================
% -------------------------------------------------
%%% System
% -------------------------------------------------
%%% Primary
uranus    = bodies.uranus;

%%% Secondaries
cordelia  = bodies.cordelia;
ophelia   = bodies.ophelia;

%%% Normalizing constants
rNorm_cor = cordelia.a;            % n <-> km
tNorm_cor = 1 / cordelia.meanMot;  % n <-> sec
vNorm_cor = rNorm_cor / tNorm_cor; % n <-> km/sec

rNorm_oph = ophelia.a;             % n <-> km
tNorm_oph = 1 / ophelia.meanMot;   % n <-> sec
vNorm_oph = rNorm_oph / tNorm_oph; % n <-> km/sec

%%% Equillibrium Points
rLPs_cor_n = EquilibriumPoints(cordelia.MR);
rLPs_oph_n  = EquilibriumPoints(ophelia.MR);

% ------------------------------------------------- 
%%% Integration Options
% -------------------------------------------------
%%% Tolerance
tol = 1e-13;

%%% Options
options = odeset('RelTol',tol,'AbsTol',tol);
options_yEqualsZero = odeset('Events',@event_yEqualsZero,'RelTol',tol,'AbsTol',tol);

%%% Create parameter struct for integration
prms_cor.u    = cordelia.MR;
prms_cor.R1_n = uranus.R/rNorm_cor;
prms_cor.R2_n = cordelia.R/rNorm_cor;
prms_cor.J21  = uranus.J2;
prms_cor.J22  = 0;

prms_oph.u  = ophelia.MR;
prms_oph.R1_n = uranus.R/rNorm_oph;
prms_oph.R2_n = ophelia.R/rNorm_oph;
prms_oph.J21  = uranus.J2;
prms_oph.J22  = 0;

% % ========================================================================
% %%% Perturb Lagrange points to compute manifolds
% % ========================================================================
% % -------------------------------------------------
% %%% Perturb the lagrange points and create new initial conditions
% % -------------------------------------------------
% %%% Set perturbation vector
% pert_n = [1,0].*1e-8;
% 
% %%% Add perturbation vector to L1 of Ophelia and L2 of Cordelia
% X0_cor_part_n = [rLPs_cor_n(2,1:2) + pert_n, zeros(1,4)]';
% X0_oph_part_n = [rLPs_oph_n(1,1:2) - pert_n, zeros(1,4)]';

% ========================================================================
%%% Create initial conditions for ring particles
% ========================================================================
% -------------------------------------------------
%%% Data from paper
% -------------------------------------------------
% Mean distance: 51149.3 km
% velocity: 10.650 km/s
% Period: 8.3823 hr
% width: 20 to 96 km
% e: 0.00794
% i: 0
a_ring_km        = 51149.3; % km
e_ring           = 0.00794; % 
minWidth_ring_km = 26;      % km
maxWidth_ring_km = 96;      % km
% -------------------------------------------------
%%% Checking minimum and maximum radius of ring particles
% -------------------------------------------------
%%% Calculating periapsis and apoapsis distances
rp_ring_km = a_ring_km * (1 - e_ring); % km
ra_ring_km = a_ring_km * (1 + e_ring); % km

%%% Calculating minimum and maximum theoretical radius for ring
minRadius_ring_km = rp_ring_km - maxWidth_ring_km; % km
maxRadius_ring_km = ra_ring_km + maxWidth_ring_km; % km


% -------------------------------------------------
%%% Create Fake Epsilon Ring
% -------------------------------------------------
nPoints = 10000;
innerEdge = [ones(nPoints,1), zeros(nPoints,2)].*minRadius_ring_km;
outerEdge = [ones(nPoints,1), zeros(nPoints,2)].*maxRadius_ring_km;

thetas = linspace(0,2*pi,nPoints);
for kk = 2:nPoints
    innerEdge(kk,:) = R3_rot([innerEdge(1,:)],thetas(kk));
    outerEdge(kk,:) = R3_rot([outerEdge(1,:)],thetas(kk));
end


% -------------------------------------------------
%%% Create an inertial initial position for a ring particle
% -------------------------------------------------
%%% Set other initial orbital elements
i_ring_deg    = 0;
raan_ring_deg = 0;
w_ring_deg    = 0;
ta_ring_deg   = 0;

%%% Convert orbital elements to an inertial state centered around primary
[r0_ring_inertial_dim, v0_ring_inertial_dim] = COE2RV(a_ring_km, e_ring, i_ring_deg, raan_ring_deg, w_ring_deg, ta_ring_deg, uranus.u);
r0_ring_inertial_dim = r0_ring_inertial_dim';
v0_ring_inertial_dim = v0_ring_inertial_dim';

%%% Getting conditions for maximum-energy ring particle
[r0_maxCircRing_inertial_dim, v0_maxCircRing_inertial_dim] = COE2RV(maxRadius_ring_km, 0, i_ring_deg, raan_ring_deg, w_ring_deg, ta_ring_deg, uranus.u);
r0_maxCircRing_inertial_dim = r0_maxCircRing_inertial_dim';
v0_maxCircRing_inertial_dim = v0_maxCircRing_inertial_dim';

%%% Getting conditions for minimum-energy ring particle
[r0_minCircRing_inertial_dim, v0_minCircRing_inertial_dim] = COE2RV(minRadius_ring_km, 0, i_ring_deg, raan_ring_deg, w_ring_deg, ta_ring_deg, uranus.u);
r0_minCircRing_inertial_dim = r0_minCircRing_inertial_dim';
v0_minCircRing_inertial_dim = v0_minCircRing_inertial_dim';

% -------------------------------------------------
%%% Transform initial state of particle to the rotating frames
%%% of both Cordelia and Ophelia and normalize to respective dimensionless
%%% units
% -------------------------------------------------
%%% Quick warning
warning('Not accounting for barycenter ... less than a meter difference though')

%%% Convert cordelia states of test-particle from ring
X0_ring_cor_rot = X_BCI2BCR([r0_ring_inertial_dim, v0_ring_inertial_dim], 0, cordelia.meanMot);
X0_ring_cor_rot_n = [X0_ring_cor_rot(:,1:3)./rNorm_cor, X0_ring_cor_rot(:,4:6)./vNorm_cor];

%%% Convert cordelia states of theoretical maximum circular ring particle
X0_maxCircRing_cor_rot = X_BCI2BCR([r0_maxCircRing_inertial_dim, v0_maxCircRing_inertial_dim], 0, cordelia.meanMot);
X0_maxCircRing_cor_rot_n = [X0_maxCircRing_cor_rot(:,1:3)./rNorm_cor, X0_maxCircRing_cor_rot(:,4:6)./vNorm_cor];

%%% Convert cordelia states of theoretical minimum circular ring particle
X0_minCircRing_cor_rot = X_BCI2BCR([r0_minCircRing_inertial_dim, v0_minCircRing_inertial_dim], 0, cordelia.meanMot);
X0_minCircRing_cor_rot_n = [X0_minCircRing_cor_rot(:,1:3)./rNorm_cor, X0_minCircRing_cor_rot(:,4:6)./vNorm_cor];

%%% Convert ophelia states of test-particle from ring
X0_ring_oph_rot = X_BCI2BCR([r0_ring_inertial_dim, v0_ring_inertial_dim], 0, ophelia.meanMot);
X0_ring_oph_rot_n = [X0_ring_oph_rot(:,1:3)./rNorm_oph, X0_ring_oph_rot(:,4:6)./vNorm_oph];

%%% Convert ophelia states of theoretical maximum circular ring particle
X0_maxCircRing_oph_rot = X_BCI2BCR([r0_maxCircRing_inertial_dim, v0_maxCircRing_inertial_dim], 0, ophelia.meanMot);
X0_maxCircRing_oph_rot_n = [X0_maxCircRing_oph_rot(:,1:3)./rNorm_oph, X0_maxCircRing_oph_rot(:,4:6)./vNorm_oph];

%%% Convert ophelia states of theoretical minimum circular ring particle
X0_minCircRing_oph_rot = X_BCI2BCR([r0_minCircRing_inertial_dim, v0_minCircRing_inertial_dim], 0, ophelia.meanMot);
X0_minCircRing_oph_rot_n = [X0_minCircRing_oph_rot(:,1:3)./rNorm_oph, X0_minCircRing_oph_rot(:,4:6)./vNorm_oph];
% ========================================================================
%%% Integrating particle in each system (cordelia and ophelia)
% ========================================================================
% ------------------------------------------------- 
%%% Final prep for integration
% -------------------------------------------------
%%% Specify time bounds
t0 = 0;
% tf_ringParticle = 4*pi;
% tf_ringParticle = 2000*pi;
tf_ringParticle = 16000*pi;
% tf_ringParticle = 400*pi;
% tf_ringParticle = 4000*pi;

%%% Create time vector
time0_n = [t0,tf_ringParticle];
time0_n_vec = linspace(time0_n(1),time0_n(end),2000000);

% ------------------------------------------------- 
%%% Integrate and store results
% -------------------------------------------------
%%% Integrate 
[t_ring_cor_n, X_ring_cor_rot_n] = ode113(@Int_CR3Bn, time0_n_vec, X0_ring_cor_rot_n, options, prms_cor);
[t_ring_oph_n, X_ring_oph_rot_n] = ode113(@Int_CR3Bn, time0_n_vec, X0_ring_oph_rot_n, options, prms_oph);

% ------------------------------------------------- 
%%% Computing Jacobi Constants
% -------------------------------------------------
%%% Test particles from ring in each system
[JCs_cor_ring_vec] = JacobiConstantCalculator(cordelia.MR, X_ring_cor_rot_n(:,1:3), X_ring_cor_rot_n(:,4:6));
[JCs_oph_ring_vec] = JacobiConstantCalculator(ophelia.MR, X_ring_oph_rot_n(:,1:3), X_ring_oph_rot_n(:,4:6));

%%% Maximum and minimum theoretical ring energies in each system
[JCs_cor_maxCircRing] = JacobiConstantCalculator(cordelia.MR, X0_maxCircRing_cor_rot_n(1:3), X0_maxCircRing_cor_rot_n(4:6));
[JCs_cor_minCircRing] = JacobiConstantCalculator(cordelia.MR, X0_minCircRing_cor_rot_n(1:3), X0_minCircRing_cor_rot_n(4:6));
[JCs_oph_maxCircRing] = JacobiConstantCalculator(ophelia.MR, X0_maxCircRing_oph_rot_n(1:3), X0_maxCircRing_oph_rot_n(4:6));
[JCs_oph_minCircRing] = JacobiConstantCalculator(ophelia.MR, X0_minCircRing_oph_rot_n(1:3), X0_minCircRing_oph_rot_n(4:6));

%%% Important Lagrange points
[JC_L2_cor] = JacobiConstantCalculator(cordelia.MR,rLPs_cor_n(2,1:3),[0, 0, 0]);
[JC_L1_oph] = JacobiConstantCalculator(ophelia.MR,rLPs_oph_n(1,1:3),[0, 0, 0]);

% ========================================================================
%%% Transform rings to Inertial
% ========================================================================
%%% Transforming states to normalized inertial
X_ring_cor_inertial_n = X_BCR2BCI(X_ring_cor_rot_n, t_ring_cor_n, 1 );
X_ring_oph_inertial_n = X_BCR2BCI(X_ring_oph_rot_n, t_ring_oph_n, 1 );

%%% Dimensionalizing inertial states
X_ring_cor_inertial = [X_ring_cor_inertial_n(:,1:3).*rNorm_cor, X_ring_cor_inertial_n(:,4:6).*vNorm_cor];
X_ring_oph_inertial = [X_ring_oph_inertial_n(:,1:3).*rNorm_oph, X_ring_oph_inertial_n(:,4:6).*vNorm_oph];

% ========================================================================
%%% Plot Results
% ========================================================================
figure
subplot(2,2,1); hold all
title('Cordelia, Rotating')
plot(X_ring_cor_rot_n(:,1),X_ring_cor_rot_n(:,2),'b.')
plot(X_ring_cor_rot_n(1,1),X_ring_cor_rot_n(1,2),'bo')
plot(X_ring_cor_rot_n(end,1),X_ring_cor_rot_n(end,2),'bx')
plot(rLPs_cor_n(1:2,1), rLPs_cor_n(1:2,2), '^', 'markeredgecolor', colors.std.blue, 'markerfacecolor', colors.std.ltblue, 'markersize',8)
plotBody2(cordelia.R/rNorm_cor, [1-cordelia.MR,0], colors.std.black, colors.std.black, 1, 1)
plotBody2(uranus.R/rNorm_cor,   [-cordelia.MR,0],  colors.std.cyan,  colors.std.black, 1, 0.5)
axis equal
PlotBoi2('$x_{rot}$','$y_{rot}$',18,'LaTex')

subplot(2,2,2); hold all
title('Ophelia, Rotating')
plot(X_ring_oph_rot_n(:,1),X_ring_oph_rot_n(:,2),'m.')
plot(X_ring_oph_rot_n(1,1),X_ring_oph_rot_n(1,2),'mo')
plot(X_ring_oph_rot_n(end,1),X_ring_oph_rot_n(end,2),'mx')
plot(rLPs_oph_n(1:2,1), rLPs_oph_n(1:2,2), '^', 'markeredgecolor', colors.std.blue, 'markerfacecolor', colors.std.ltblue, 'markersize',8)
plotBody2(ophelia.R/rNorm_oph, [1-ophelia.MR,0], colors.std.black, colors.std.black, 1, 1)
plotBody2(uranus.R/rNorm_oph,  [-ophelia.MR,0],  colors.std.cyan,  colors.std.black, 1, 0.5)
axis equal
PlotBoi2('$x_{rot}$','$y_{rot}$',18,'LaTex')

subplot(2,2,3); hold all
title('Cordelia, Inertial')
plot(X_ring_cor_inertial_n(:,1),X_ring_cor_inertial_n(:,2),'b.')
plot(X_ring_cor_inertial_n(1,1),X_ring_cor_inertial_n(1,2),'bo')
plot(X_ring_cor_inertial_n(end,1),X_ring_cor_inertial_n(end,2),'bx')
plotBody2(uranus.R/rNorm_cor,   [0,0],  colors.std.cyan,  colors.std.black, 1, 0.5)
axis equal
PlotBoi2('$x$','$y$',18,'LaTex')

subplot(2,2,4); hold all
title('Ophelia, Inertial')
plot(X_ring_oph_inertial_n(:,1),X_ring_oph_inertial_n(:,2),'m.')
plot(X_ring_oph_inertial_n(1,1),X_ring_oph_inertial_n(1,2),'mo')
plot(X_ring_oph_inertial_n(end,1),X_ring_oph_inertial_n(end,2),'mx')
plotBody2(uranus.R/rNorm_oph,   [0,0],  colors.std.cyan,  colors.std.black, 1, 0.5)
axis equal
PlotBoi2('$x$','$y$',18,'LaTex')


%% =======================================================================
%%% Look at manifolds of Lyapunov orbits from Cordelia L2 and Ophelia L1
% ========================================================================
if runManifolds == 1
% -------------------------------------------------
%%% General Manifold Setup
% -------------------------------------------------
%%% Number of points along PO to observe manifold from
n_manifolds = 50;

%%% Integration time for unstable manifolds
tf_man = 8000*2*pi;

%%% Initialize STM
stm0        = eye(6);
stm0_colVec = reshape(stm0,36,1);

run_cor_manifold = 1;
run_oph_manifold = 1;
% -------------------------------------------------
%%% Set IC for Lyapunov orbit and integrate the full PO
% -------------------------------------------------
% %%% Largest(ish) Lyapunov at Cordelia L2 (smallest JC)
PO_cor_L2_Lyap = [1.000745138331199;
                     0.000000013561826;
                     0.000000000000000;
                     0.000000014923329;
                     -0.001468288082842;
                     0.000000000000000];
Tp_cor_L2_Lyap = 4.019928182484143;

%%% Smallest(ish) Lyapunov at Cordelia L2 (Largest JC)
% PO_cor_L2_Lyap = [1.000554317117410;
%                  -0.000000000000000;
%                  0.000000000000000;
%                  -0.000000000014150;
%                  -0.000009259839589;
%                  0.000000000000000];
% Tp_cor_L2_Lyap = 3.034220417738762;

% %%% Largest(ish) Lyapunov at Ophelia L1 (smallest JC)
PO_oph_L1_Lyap = [0.999791280335056;
                 -0.000000000000000;
                 0.000000000000000;
                 -0.000000000008347;
                 -0.002209798953327;
                 0.000000000000000];
Tp_oph_L1_Lyap = 3.909793685141573;


%%% Smallest(ish) Lyapunov at Ophelia L1 (Largest JC)
% PO_oph_L1_Lyap = [0.999412754405857;
%                  -0.000000000000000;
%                  0.000000000000000;
%                  0.000000000000900;
%                  -0.000005607138979;
%                  0.000000000000000];
% Tp_oph_L1_Lyap = 3.031762317014002;


%%% Integrating Lyapunov orbit to get full PO
if run_cor_manifold == 1
[~, X_cor_L2_Lyap_full] = ode113(@Int_CR3Bn, [0, Tp_cor_L2_Lyap], PO_cor_L2_Lyap(1:6), options, prms_cor);
end
if run_oph_manifold == 1
    [~, X_oph_L1_Lyap_full] = ode113(@Int_CR3Bn, [0, Tp_oph_L1_Lyap], PO_oph_L1_Lyap(1:6), options, prms_oph);
end


%%% Plot PO
if run_cor_manifold == 1
    figure(100); hold all
    plot3(X_cor_L2_Lyap_full(:,1),X_cor_L2_Lyap_full(:,2),X_cor_L2_Lyap_full(:,3),'k','linewidth',2)
    axis equal
end
if run_oph_manifold == 1
    figure(200); hold all
    plot3(X_oph_L1_Lyap_full(:,1),X_oph_L1_Lyap_full(:,2),X_oph_L1_Lyap_full(:,3),'k','linewidth',2)
    axis equal
end

%%% Integrating L2 Lyapunov orbit at cordelia to get "n_manifolds" nodes
% if run_cor_manifold == 1
    [~, X_cor_L2_Lyap_nodes] = ode113(@Int_CR3Bn, linspace(0, Tp_cor_L2_Lyap, n_manifolds+1), PO_cor_L2_Lyap(1:6), options, prms_cor);
% end
% if run_oph_manifold == 1
    [~, X_oph_L1_Lyap_nodes] = ode113(@Int_CR3Bn, linspace(0, Tp_oph_L1_Lyap, n_manifolds+1), PO_oph_L1_Lyap(1:6), options, prms_oph);
% end
    % -------------------------------------------------
    %%% If finding manifolds from the unstable eigenvector of the monodromy
    %%% matrix
    % -------------------------------------------------
    if manifoldsFromMonodromy == 1
        %%% Scaling for perturbing the IC in the direction of the unstable
        %%% eigenvector
        pertScale = 1e-8;

% %         %%% Set options to try and cut out unstable manifolds that go to
% %         %%% interior region
% %         prms_cor.xRegion = 0.9998;
% %         prms_cor.yRegion = 1e-3;
% %         options_xyRegionEscape = odeset('Events',@event_xyRegionEscape_CR3Bn,'RelTol',tol,'AbsTol',tol);
        
        %%% Preallocating for parfor loop
        trajData_cor = cell(n_manifolds,1);
        trajData_oph = cell(n_manifolds,1);
        
        %%% Integrate full PO from each node to get monodromy matrix and
        %%% corresponding unstable eigenvector to perturb the PO and reveal the
        %%% manifold with
        parfor kk = 1:n_manifolds
            if run_cor_manifold == 1
                %%% Current X0 is the 'kk'th node of the PO
                X0_cor_L2_Lyap_m_i = X_cor_L2_Lyap_nodes(kk,:)';

                %%% Integrate the current node with an STM
                [~, Xstm_cor_L2_Lyap_i] = ode113(@Int_CR3BnSTM, [0, Tp_cor_L2_Lyap], [X0_cor_L2_Lyap_m_i; stm0_colVec], options, prms_cor);

                %%% Get monodromy matrix and determine its stable and unstable
                %%% eigenvectors
                stm_tf_t0 = reshape(Xstm_cor_L2_Lyap_i(end,7:42),6,6);
                monodromy = stm_tf_t0;
                [eVec_stable, eVec_unstable] = getStableAndUnstableEigenvectors(monodromy);

                %%% Create the perturbed initial condition and integrate the manifold
                X0_cor_L2_Lyap_stable_i = X0_cor_L2_Lyap_m_i + eVec_stable.*pertScale;

                %%% Integrate the manifold
                [~, Xstm_cor_L2_Lyap_stable] = ode113(@Int_CR3Bn, [tf_man, 0], X0_cor_L2_Lyap_stable_i, options, prms_cor);

                %%% Adding a 'NaN' so trajectories don't combine later on the
                %%% plot
                Xstm_cor_L2_Lyap_stable = [Xstm_cor_L2_Lyap_stable; NaN(1,6)];

                %%% Storing data
                trajData_cor{kk}.Xstm_cor_L2_Lyap_stable = Xstm_cor_L2_Lyap_stable(:,1:3);
                trajData_cor{kk}.X0_cor_L2_Lyap_stable = Xstm_cor_L2_Lyap_stable(1,1:6);
            end
            
            if run_oph_manifold == 1
                %%% Current X0 is the 'kk'th node of the PO
                X0_oph_L1_Lyap_m_i = X_oph_L1_Lyap_nodes(kk,:)';

                %%% Integrate the current node with an STM
                [~, Xstm_oph_L1_Lyap_i] = ode113(@Int_CR3BnSTM, [0, Tp_oph_L1_Lyap], [X0_oph_L1_Lyap_m_i; stm0_colVec], options, prms_oph);

                %%% Get monodromy matrix and determine its stable and unstable
                %%% eigenvectors
                stm_tf_t0 = reshape(Xstm_oph_L1_Lyap_i(end,7:42),6,6);
                monodromy = stm_tf_t0;
                [eVec_stable, eVec_unstable] = getStableAndUnstableEigenvectors(monodromy);

                %%% Create the perturbed initial condition and integrate the manifold
                X0_oph_L1_Lyap_stable_i = X0_oph_L1_Lyap_m_i + eVec_stable.*pertScale;

                %%% Integrate the manifold
                [~, Xstm_oph_L1_Lyap_stable] = ode113(@Int_CR3Bn, [tf_man, 0], X0_oph_L1_Lyap_stable_i, options, prms_oph);

                %%% Adding a 'NaN' so trajectories don't combine later on the
                %%% plot
                Xstm_oph_L1_Lyap_stable = [Xstm_oph_L1_Lyap_stable; NaN(1,6)];

                %%% Storing data
                trajData_oph{kk}.Xstm_oph_L1_Lyap_stable = Xstm_oph_L1_Lyap_stable(:,1:3);
                trajData_oph{kk}.X0_oph_L1_Lyap_stable = Xstm_oph_L1_Lyap_stable(1,1:6);
            end
        end % parfor
        
        if run_cor_manifold == 1
            %%% Restructuring trajectory data from parfor - now one large
            %%% struct array
            trajData_cor_structArray = [trajData_cor{:}];

            %%% Grabbing specific data sets and turning them into single
            %%% matrices
            X_cor_L2_Lyap_stable_n = [];
            X0s_cor_L2_Lyap_stable_n = [];

            for kk = 1:n_manifolds
                X_cor_L2_Lyap_stable_n = [X_cor_L2_Lyap_stable_n; [trajData_cor_structArray(kk).Xstm_cor_L2_Lyap_stable]];
                X0s_cor_L2_Lyap_stable_n = [X0s_cor_L2_Lyap_stable_n; [trajData_cor_structArray(kk).X0_cor_L2_Lyap_stable]];
            end

            %%% Plotting manifolds
            figure(100)
            plot3(X_cor_L2_Lyap_stable_n(:,1),X_cor_L2_Lyap_stable_n(:,2),X_cor_L2_Lyap_stable_n(:,3),'g')
        end
        if run_oph_manifold
            %%% Restructuring trajectory data from parfor - now one large
            %%% struct array
            trajData_oph_structArray = [trajData_oph{:}];

            %%% Grabbing specific data sets and turning them into single
            %%% matrices
            X_oph_L1_Lyap_stable_n = [];
            X0s_oph_L1_Lyap_stable_n = [];

            for kk = 1:n_manifolds
                X_oph_L1_Lyap_stable_n = [X_oph_L1_Lyap_stable_n; [trajData_oph_structArray(kk).Xstm_oph_L1_Lyap_stable]];
                X0s_oph_L1_Lyap_stable_n = [X0s_oph_L1_Lyap_stable_n; trajData_oph_structArray(kk).X0_oph_L1_Lyap_stable];
            end

            %%% Plotting manifolds
            figure(200)
            plot3(X_oph_L1_Lyap_stable_n(:,1),X_oph_L1_Lyap_stable_n(:,2),X_oph_L1_Lyap_stable_n(:,3),'g')
        end
        

        
    end % manifoldsFromMonodromy
    
% % % % %     % -------------------------------------------------
% % % % %     %%% If finding manifolds from small, constant perturbation
% % % % %     % -------------------------------------------------
% % % % %     if manifoldsFromPert == 1
% % % % %         %%% Choose perturbation vector
% % % % %         XPert_cor = [1; 0; 0; 0; 0; 0].*(1e-1)/rNorm_cor;
% % % % %         
% % % % %         %%% Preallocating for parfor loop
% % % % %         trajData_cor = cell(n_manifolds,1);
% % % % %         
% % % % %         parfor kk = 1:n_manifolds
% % % % %             %%% Grab the initial state of the current 
% % % % %             X0_cor_L2_Lyap_m_i = X_cor_L2_Lyap_nodes(kk,:)'
% % % % %             
% % % % %             %%% Perturb in the "plus" and "minus" directions
% % % % %             X0_cor_L2_Lyap_p_i = X0_cor_L2_Lyap_m_i + XPert_cor;
% % % % %             X0_cor_L2_Lyap_m_i = X0_cor_L2_Lyap_m_i - XPert_cor;
% % % % %             
% % % % %             %%% Integrate each case backward in time (stable manifolds)
% % % % %             [~, X_cor_L2_Lyap_p_stable_i] = ode113(@Int_CR3Bn, [tf_cor_man, 0], X0_cor_L2_Lyap_p_i, options, prms_cor);
% % % % %             [~, X_cor_L2_Lyap_m_stable_i] = ode113(@Int_CR3Bn, [tf_cor_man, 0], X0_cor_L2_Lyap_m_i, options, prms_cor);
% % % % %             
% % % % %             %%% Adding a 'NaN' so trajectories don't combine later on the
% % % % %             %%% plot
% % % % %             X_cor_L2_Lyap_p_stable_i = [X_cor_L2_Lyap_p_stable_i; NaN(1,6)];
% % % % %             X_cor_L2_Lyap_m_stable_i = [X_cor_L2_Lyap_m_stable_i; NaN(1,6)];
% % % % %             
% % % % %             %%% Store data
% % % % %             trajData_cor{kk}.X_cor_L2_Lyap_p_stable = X_cor_L2_Lyap_p_stable_i(:,1:3);
% % % % %             trajData_cor{kk}.X_cor_L2_Lyap_m_stable = X_cor_L2_Lyap_m_stable_i(:,1:3);
% % % % %             
% % % % %         end % parfor
% % % % %         
% % % % %         %%% Restructuring trajectory data from parfor - now one large
% % % % %         %%% struct array
% % % % %         trajData_cor_structArray = [trajData_cor{:}];
% % % % %         
% % % % %         %%% Grabbing specific data sets and turning them into single
% % % % %         %%% matrices
% % % % %         X_cor_L2_Lyap_p_stable = [];
% % % % %         X_cor_L2_Lyap_m_stable = [];
% % % % %         
% % % % %         for kk = 1:n_manifolds
% % % % %             X_cor_L2_Lyap_p_stable = [X_cor_L2_Lyap_p_stable; [trajData_cor_structArray(kk).X_cor_L2_Lyap_p_stable]];
% % % % %             X_cor_L2_Lyap_m_stable = [X_cor_L2_Lyap_m_stable; [trajData_cor_structArray(kk).X_cor_L2_Lyap_m_stable]];
% % % % %         end
% % % % %         
% % % % %         %%% Plotting manifolds
% % % % %         plot3(X_cor_L2_Lyap_p_stable(:,1),X_cor_L2_Lyap_p_stable(:,2),X_cor_L2_Lyap_p_stable(:,3),'g')
% % % % %         plot3(X_cor_L2_Lyap_m_stable(:,1),X_cor_L2_Lyap_m_stable(:,2),X_cor_L2_Lyap_m_stable(:,3),'g')
% % % % %     end












end % runManifolds
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

%%
figure('position',[283 171 852 627]); hold all
axis equal
PlotBoi2('$X$, $km$','$Y$, $km$',20,'LaTex')

r = minRadius_ring_km;
R = maxRadius_ring_km;
xf = 0;
Xf = 0;
yf = 0;
Yf = 0;
% Creates a annulus patch object and returns the handle.  Input arguments 
% are the inner radius, outer radius, inner x offset, outer x offset, inner
% y offset and outer Y offset.  Changes to the edgecolor and linestyle are
% allowed, and will preserve the correct look of the annulus
t = linspace(0,2*pi,200);
x = xf + r*cos(t);
y = yf + r*sin(t);
X = Xf + R*cos(t);
Y = Yf + R*sin(t);
P = patch([x X],[y Y],colors.std.grey,'linestyle','non','facealph',.5);
L(1) = line(x,y,'color',colors.std.white);
L(2) = line(X,Y,'color',colors.std.white);
axis equal
% plistener(P,'edgecolor',@edgec) % listeners for changes to props.
% plistener(P,'linestyle',@lnstl)


r_cor_L2_Lyap_stable = X_cor_L2_Lyap_stable_n(:,1:3).*rNorm_cor;
r_oph_L2_Lyap_stable = X_oph_L1_Lyap_stable_n(:,1:3).*rNorm_oph;
 
r_ring_cor_rot = X_ring_cor_rot_n(:,1:3).*rNorm_cor;
r_ring_oph_rot = X_ring_oph_rot_n(:,1:3).*rNorm_oph;


% 
% fill(outerEdge(:,1),outerEdge(:,2),colors.std.ltgrey)
% fill(innerEdge(:,1),innerEdge(:,2),colors.std.white)

plot3(r_cor_L2_Lyap_stable(:,1),r_cor_L2_Lyap_stable(:,2),r_cor_L2_Lyap_stable(:,3),'g')
plot3(r_oph_L2_Lyap_stable(:,1),r_oph_L2_Lyap_stable(:,2),r_oph_L2_Lyap_stable(:,3),'g')

plot3(r_ring_cor_rot(:,1),r_ring_cor_rot(:,2),r_ring_cor_rot(:,3),'k')
plot3(r_ring_oph_rot(:,1),r_ring_oph_rot(:,2),r_ring_oph_rot(:,3),'r')














function [] = plistener(ax,prp,func)
  % Sets the properties. From proplistener by Yair Altman.
  psetact = 'PropertyPostSet';
  hC = handle(ax);
  hSrc = hC.findprop(prp);
  hl = handle.listener(hC, hSrc, psetact, {func,ax});
  p = findprop(hC, 'Listeners__');
  if isempty(p)
      p = schema.prop(hC, 'Listeners__', 'handle vector');
      set(p,'AccessFlags.Serialize', 'off', ...
          'AccessFlags.Copy', 'off',...
          'FactoryValue', [], 'Visible', 'off');
  end
  hC.Listeners__ = hC.Listeners__(ishandle(hC.Listeners__));
  hC.Listeners__ = [hC.Listeners__; hl];
end
function [] = edgec(varargin)
  St = get(varargin{3},'edgecolor');
  set(L,'color',St)
end
function [] = lnstl(varargin)
  St = get(varargin{3},'linestyle');
  set(varargin{3},'linestyle','none')
  set(L,'linestyle',St)
end






