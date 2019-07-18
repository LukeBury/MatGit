% ========================================================================
%%% Description
% ========================================================================
% 

% Created: 
% Author : Luke Bury, luke.bury@colorado.edu
% ========================================================================
%%% Initialization
% ========================================================================
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

% ========================================================================
%%% Run Switches
% ========================================================================


% ========================================================================
%%% Setup
% ========================================================================
% -------------------------------------------------
%%% System
% -------------------------------------------------
%%% Choose bodies
primary = bodies.jupiter;   secondary = bodies.europa;
% primary = bodies.jupiter;   secondary = bodies.ganymede;
% primary = bodies.jupiter;   secondary = bodies.callisto;
% primary = bodies.saturn;    secondary = bodies.enceladus;
% primary = bodies.saturn;    secondary = bodies.titan;
% primary = bodies.neptune;   secondary = bodies.triton;

%%% Normalizing constants
rNorm = secondary.a;         % n <-> km
tNorm = 1/secondary.meanMot; % n <-> sec
vNorm = rNorm / tNorm;       % n <-> km/sec

%%% Shortcut variables
mu   = secondary.MR;
R2_n = secondary.R_n;

%%% Equillibrium Points
rLPs_n = EquilibriumPoints(mu);

% ------------------------------------------------- 
%%% Integration Options
% -------------------------------------------------
tol = 1e-13;
options = odeset('RelTol',tol,'AbsTol',tol);
prms.u    = mu;
prms.R2_n = R2_n;
t0 = 0;

% ------------------------------------------------- 
%%% Initial conditions
% -------------------------------------------------
% --------------------------
% X0 and Tp corresponding to L2 Lyapunov biffurcation at Europa
% --------------------------
X0_n = [1.017026611070739; 0.000000014441540; 0;...
        0.000000002451836; 0.020022207280852; 0];
T0_PO_n = 3.124263922352081;


% ========================================================================
%%% Find new family
% ========================================================================
% ------------------------------------------------- 
%%% Perturb X0 in z direction to find biffurcated family
% -------------------------------------------------
X0_guess_n = X0_n + [0; 0; 1e-6; 0; 0; 1e-6];

% ------------------------------------------------- 
%%% Integrate first guess
% -------------------------------------------------
%%% Initial stm is identity (t0 to t0)
stm0 = reshape(eye(6),36,1);

%%% Set full initial state that includes error
X0_guess_n = [X0_guess_n; stm0];

% X0_n = [X0_n; stm0];

% ------------------------------------------------- 
%%% Shooter Options
% -------------------------------------------------
n_iterations = 30;
stateTarget = X0_guess_n(1:6);
dX0_scale = 1;

% ------------------------------------
%%% Initial integration
% ------------------------------------
%%% Integrate trajectory
time0_n = [0, T0_PO_n];
[time_n, X_BCR_n_guess] = ode113(@Int_CR3BnSTM, time0_n, X0_guess_n, options, prms);

% ------------------------------------
%%% Setting up for shooting
% ------------------------------------
%%% Initialize new variables to iterate over
X0_new = X0_guess_n;
X_BCR_n_new = X_BCR_n_guess;


%%% Iterate
for kk = 1:n_iterations
    % ------------------------------------
    %%% Correct state
    % ------------------------------------
    %%% Set error in final state
%     dXf = X_BCR_n_new(end,1:6)' - stateTarget;
    dXf = X_BCR_n_new(end,1:6)' - X_BCR_n_new(1,1:6)';

    %%% Grab STM(tf, t0) and invert it to get STM(t0, tf)
    stm_tf_t0 = reshape(X_BCR_n_new(end,7:42),6,6);

    stm_t0_tf = inv(stm_tf_t0); % state

    %%% Map dXf to t0
    dX0 = stm_t0_tf * dXf;

    %%% Correct X0 with dX0
    X0_new = [X0_new(1:6) - dX0.*dX0_scale; stm0]; % full state

    % ------------------------------------
    %%% Integrate new correction
    % ------------------------------------
    %%% Integrate new state
    [time_n, X_BCR_n_new] = ode113(@Int_CR3BnSTM, time0_n, X0_new, options, prms);

    %%% Print results from current iteration
    fprintf('-------------------------------------------- %1d\n',kk)
    fprintf('Norm of dXf:  \t%1.4e\n',norm(dXf))

end

figure; hold all
plot3(X_BCR_n_new(:,1),X_BCR_n_new(:,2),X_BCR_n_new(:,3))



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
toc(ticWhole)
















