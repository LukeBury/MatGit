% ========================================================================
%%% Description
% ========================================================================
% Finding POs by also correcting time vector

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
% primary = bodies.earth;     secondary = bodies.moon;
% primary = bodies.jupiter;   secondary = bodies.europa;
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
X0_guess_n = X0_n - [0; 0; 1e-5; 0; 0; 1e-5];

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
F_new = [X0_new(1:6); T0_PO_n];


%%% Iterate
for kk = 1:n_iterations
    % ------------------------------------
    %%% Correct state
    % ------------------------------------
    %%% Set constraint vector
    c = X_BCR_n_new(end,1:6)' - X_BCR_n_new(1,1:6)';

    %%% Grab STM(tf, t0)
    stm_tf_t0 = reshape(X_BCR_n_new(end,7:42),6,6);
    
    %%% Grab EOM at Xf
    Xdot_f = Int_CR3Bn(0,X_BCR_n_new(end,1:6)',prms);
    
    %%% Create DF matrix
    DF = [stm_tf_t0 - eye(6), Xdot_f];
    
    %%% Correct Free Variables
    F_new = F_new - (DF'/(DF*(DF')))*c;
    X0_new = [F_new(1:6); stm0];
    T_PO_new = F_new(7);

    % ------------------------------------
    %%% Integrate new correction
    % ------------------------------------
    %%% Integrate new state
    [time_n, X_BCR_n_new] = ode113(@Int_CR3BnSTM, [0, T_PO_new], X0_new, options, prms);

    %%% Print results from current iteration
    fprintf('-------------------------------------------- %1d\n',kk)
    fprintf('Norm of constraint vector:  \t%1.4e\n',norm(c))

end

figure; hold all
plot3(X_BCR_n_new(:,1),X_BCR_n_new(:,2),X_BCR_n_new(:,3))




% -------------------------------------------------
%%% Integrating solution to x-y plane crossing for state simplicity
% -------------------------------------------------
options_XZStop = odeset('Event',@event_yEqualsZero,'RelTol',tol,'AbsTol',tol);
[time_n_XY, X_BCR_n_XY] = ode113(@Int_CR3Bn, [0, T_PO_new], X0_new(1:6), options_XZStop, prms);









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

















