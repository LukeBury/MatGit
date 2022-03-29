% ========================================================================
%%% Description
% ========================================================================
% 

% Created: 
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
%%% Jupiter
bodies.jupiter.name  = 'jupiter';
bodies.jupiter.mass  = 1.89819e27; % kg
bodies.jupiter.u     = 1.268e8; % km^3 * s^-2
bodies.jupiter.color = [235,214,173]./255;
bodies.jupiter.a     = 778.57e6; % km, https://nssdc.gsfc.nasa.gov/planetary/factsheet/jupiterfact.html
bodies.jupiter.R     = 69911; % km

%%% Europa
bodies.europa.name    = 'europa';
bodies.europa.color   = [0, 1, 1];
bodies.europa.mass    = 4.799e22; % kg
bodies.europa.u       = 3203.413216; % km^3 / s^2
bodies.europa.a       = 671100; % km
bodies.europa.R       = 1560.8; % km 
bodies.europa.R_n     = bodies.europa.R / bodies.europa.a;
bodies.europa.meanMot = 2*pi/(3.551181*86400); % rad/s
bodies.europa.Tp      = 2*pi/bodies.europa.meanMot; % sec
bodies.europa.MR      = bodies.europa.mass / (bodies.europa.mass + bodies.jupiter.mass); % Mass ratio w/ primary

%%% Choose bodies
primary = bodies.jupiter;   secondary = bodies.europa;


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
options_yStop = odeset('Event',@event_yEqualsZero,'RelTol',tol,'AbsTol',tol);
prms.u    = mu;
prms.R2_n = R2_n;
t0 = 0;

% ------------------------------------------------- 
%%% Initial conditions
% -------------------------------------------------
% --------------------------
% X0 and Tp corresponding to L2 Lyapunov biffurcation at Europa
% --------------------------
X0_guess_n = [rLPs_n(2,1) - 0.003; 0; 0; 0; 0.018; 0];
% X0_guess_n = [rLPs_n(1,1) - 0.003; 0; 0; 0; 0.023; 0];
Tp_n = 3.2;
        

% ------------------------------------------------- 
%%% Integrate first guess
% -------------------------------------------------
%%% Initial stm is identity (t0 to t0)
stm0 = reshape(eye(6),36,1);

%%% Set full initial state that includes error
X0_guess_n = [X0_guess_n; stm0];

% ------------------------------------------------- 
%%% Shooter Options
% -------------------------------------------------
n_iterations = 10;

% ------------------------------------
%%% Initial integration
% ------------------------------------
%%% Integrate trajectory
[time_n, X_BCR_n_guess] = ode113(@Int_CR3BnSTM, [0, Tp_n], X0_guess_n, options_yStop, prms);

figure; hold all
plot3(X_BCR_n_guess(:,1),X_BCR_n_guess(:,2),X_BCR_n_guess(:,3))
plot3(rLPs_n(1,1), rLPs_n(1,2), rLPs_n(1,3), '^','markeredgecolor',colors.black,'markerfacecolor',colors.ltgrn)
plot3(rLPs_n(2,1), rLPs_n(2,2), rLPs_n(2,3), '^','markeredgecolor',colors.black,'markerfacecolor',colors.ltgrn)

% ------------------------------------
%%% Setting up for shooting
% ------------------------------------
%%% Initialize new variables to iterate over
X0_new = X0_guess_n;
X_BCR_n_new = X_BCR_n_guess;
F_new = [X0_new(1:6); Tp_n];


%%% Iterate
for kk = 1:n_iterations
    dXf = X_BCR_n_new(end,1:6)' - X_BCR_n_new(1,1:6)';
    
    %%% Grab STM(tf, t0)
    stm_tf_t0 = reshape(X_BCR_n_new(end,7:42),6,6);
    stm_t0_tf = inv(stm_tf_t0);
        
    %%% Correct ydot0 with dxf
    dydot0_new = stm_tf_t0(1,5) \ dXf(1);
    
    X0_new(5) = X0_new(5) - dydot0_new;
    % ------------------------------------
    %%% Integrate new correction
    % ------------------------------------
    %%% Integrate new state
    [time_n, X_BCR_n_new] = ode113(@Int_CR3BnSTM, [0, Tp_n], X0_new, options_yStop, prms);

    %%% Print results from current iteration
    fprintf('-------------------------------------------- %1d\n',kk)
    fprintf('Norm of error vector:  \t%1.4e\n',dXf(1))

end

figure; hold all
plot3(X_BCR_n_new(:,1),X_BCR_n_new(:,2),X_BCR_n_new(:,3),'b','linewidth',2)
plot3(rLPs_n(1,1), rLPs_n(1,2), rLPs_n(1,3), '^','markeredgecolor',colors.black,'markerfacecolor',colors.ltgrn)
plot3(rLPs_n(2,1), rLPs_n(2,2), rLPs_n(2,3), '^','markeredgecolor',colors.black,'markerfacecolor',colors.ltgrn)
plotBody2(secondary.R_n, [1-secondary.MR,0,0],colors.ltblue,colors.black,1.5,1)
view(0,90)
PlotBoi3_CR3Bn(22)
axis equal





















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
















