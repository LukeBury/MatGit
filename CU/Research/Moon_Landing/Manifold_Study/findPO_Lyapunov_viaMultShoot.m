% ========================================================================
%%% Description
% ========================================================================
% Function to find Lypunov PO initial conditions by only changing the
% y-dot-0 of the first node, and full states of the remaining nodes. First
% node, or guess, should be on x-axis with only a y-velocity

% Created: 06/20/19
% Author : Luke Bury, luke.bury@colorado.edu
% ========================================================================
%%% Initialization
% ========================================================================
clear
clc
close all
mbinPath    = '~/CU_Google_Drive/Documents/MatGit/mbin';
PO_SavePath = '~/CU_Google_Drive/Documents/MatGit/CU/Research/Moon_Landing/Manifold_Study/PO_X0s';
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
run_CR3BP              = 1;
run_CR3BP_J2pJ4pJ6pJ2s = 0;


if (run_CR3BP == 1) && (run_CR3BP_J2pJ4pJ6pJ2s == 1)
    warning('These are mutually exclusive options')
    return
end

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
primary = bodies.saturn;    secondary = bodies.enceladus;
% primary = bodies.saturn;    secondary = bodies.titan;
% primary = bodies.neptune;   secondary = bodies.triton;
% primary = bodies.uranus;    secondary = bodies.cordelia;
% primary = bodies.uranus;    secondary = bodies.ophelia;

%%% Choose equillibrium
% LP = 1;
LP = 2;

%%% Normalizing constants
rNorm = secondary.a;         % n <-> km
tNorm = 1/secondary.meanMot; % n <-> sec
vNorm = rNorm / tNorm;       % n <-> km/sec

%%% Shortcut variables
mu   = secondary.MR;
R2_n = secondary.R_n;

%%% System parameters
prms.u   = mu;
prms.R1  = primary.R / rNorm;
prms.R2  = secondary.R_n;
% prms.J2p = primary.J2;
% prms.J4p = primary.J4;
% prms.J6p = primary.J6;
% prms.J2s = secondary.J2;

%%% Equillibrium Points
if run_CR3BP == 1
    rLPs_n = EquilibriumPoints(mu);
elseif run_CR3BP_J2pJ4pJ6pJ2s == 1
    rLPs_n = collinearEquilibriumPoints_ZH(prms);
end

% ------------------------------------------------- 
%%% Periodic Orbit options
% -------------------------------------------------
%%% Error tolerance for constraint vector
error_tol = 1e-13;

%%% Number of nodes for multiple shooter
n_Nodes = 20; 

%%% Maximum number of iterations for multiple shooter
iterMax = 500;

%%% Step size for multiple shooter (usually 1, sometimes 2 or maybe even 3)
stepSize = 1;

% ------------------------------------------------- 
%%% Integration Options
% -------------------------------------------------
tol      = 1e-13;
options  = odeset('RelTol',tol,'AbsTol',tol);

% ------------------------------------------------- 
%%% Guess at initial state
% -------------------------------------------------
% x0         = rLPs_n(LP,1) - 5/rNorm;
% r0         = [x0; 0; 0];
% v0_guess_n = [0; 0.0001275; 0];
% X0_guess_n = [r0; v0_guess_n];
x0         = 1.0040019600000001;
r0         = [x0; 0; 0];
v0_guess_n = [0; -0.0000696586977203; 0];
X0_guess_n = [r0; v0_guess_n];

T_guess_n = 3.0415795008451285;



% ------------------------------------------------- 
%%% Integrate and plot guess
% -------------------------------------------------
%%% Setup
stm0 = eye(6);
stm0_vec = reshape(stm0,36,1);

%%% Integrate guess
if run_CR3BP == 1
    [~, X_guess_n] = ode113(@Int_CR3Bn, linspace(0, T_guess_n,2), X0_guess_n, options, prms);
elseif run_CR3BP_J2pJ4pJ6pJ2s == 1
    [~, X_guess_n] = ode113(@Int_CR3BnSTM_J2pJ4pJ6pJ2s, linspace(0, T_guess_n,2), [X0_guess_n; stm0_vec], options, prms);
end

%%% Plot Guess
figure; hold all
plot3(X_guess_n(:,1),X_guess_n(:,2),X_guess_n(:,3),'r','linewidth',1.5)
PlotBoi3('$X_n$','$Y_n$','$Z_n$',20,'LaTex')
plot3(rLPs_n(LP,1), rLPs_n(LP,2), rLPs_n(LP,3), '^', 'markersize', 8, 'markerfacecolor', colors.std.ltgrn, 'markeredgecolor', colors.std.blue)
view(0,90)

% 989
% return
% ========================================================================
%%% Create nodes and send the guess to a multiple shooter
% ========================================================================
% ------------------------------------------------- 
%%% Integrate guess and divide into nodes for multiple shooter
% -------------------------------------------------
%%% Integrate again to discretize into nodes evenly spaced in time
if run_CR3BP == 1
    [T_nodes, X_nodes] = ode113(@Int_CR3Bn, linspace(0, T_guess_n,n_Nodes+1), X0_guess_n, options, prms);
elseif run_CR3BP_J2pJ4pJ6pJ2s == 1
    [T_nodes, X_nodes] = ode113(@Int_CR3BnSTM_J2pJ4pJ6pJ2s, linspace(0, T_guess_n,n_Nodes+1), [X0_guess_n; stm0_vec], options, prms);
end

%%% Setting initial nodes (getting rid of repeat end state)
X_nodes = X_nodes(1:n_Nodes,1:6);

%%% Create first free-variable column vector
F_new   = X_nodes';
F_new   = F_new(:);
F_new   = [F_new(5); F_new(7:end); T_guess_n];

% -------------------------------------------------
%%% Pre-While-Loop work
% -------------------------------------------------
%%% Set while-loop/convergence variables
iteration_counter = 0;
constraint_error  = 100;

% -------------------------------------------------
% Enter while loop - attempt to converge on next PO in family
% -------------------------------------------------
while (constraint_error > error_tol) && (iteration_counter < iterMax)

    %%% Count iteration
    iteration_counter = iteration_counter + 1;

    % -------------------------------------------------
    %%% Using multiple shooting - loop through nodes, integrate, and 
    %%% populate constraint vector and create DF matrix
    % -------------------------------------------------  
    %%% Integrate again to discretize into nodes evenly spaced in time
    if run_CR3BP == 1
        [DF_mat, constraints] = multShooter_stateContinuity_ydFree(n_Nodes, x0, F_new, @Int_CR3BnSTM, options, prms);
    elseif run_CR3BP_J2pJ4pJ6pJ2s == 1
        [DF_mat, constraints] = multShooter_stateContinuity_ydFree(n_Nodes, x0, F_new, @Int_CR3BnSTM_J2pJ4pJ6pJ2s, options, prms);
    end

    % -------------------------------------------------
    %%% Determine current error and create next F if necessary
    % -------------------------------------------------
    %%% Compute error
    constraint_error = norm(constraints);

    %%% Compute new free-variable vector if error is not converged
    if (constraint_error > error_tol)
        warning('off','MATLAB:nearlySingularMatrix')
        F_new = F_new - stepSize * DF_mat'*((DF_mat*(DF_mat'))\constraints);
        warning('on','MATLAB:nearlySingularMatrix')
    end

end % while constraint_error > error_tol

if constraint_error > error_tol
    warning('Solution not converged on')
end



% ========================================================================
%%% Integrate solution, plot and print result
% ========================================================================
%%% Integrate Solution
if run_CR3BP == 1
    [T_sol_n, X_sol_n] = ode113(@Int_CR3Bn, [0, F_new(end)], [x0; 0; 0; 0; F_new(1); 0], options, prms);
elseif run_CR3BP_J2pJ4pJ6pJ2s == 1
    [T_sol_n, X_sol_n] = ode113(@Int_CR3BnSTM_J2pJ4pJ6pJ2s, [0, F_new(end)], [x0; 0; 0; 0; F_new(1); 0; stm0_vec], options, prms);
end

%%% Plot Solution
figure; hold all
plot3(X_sol_n(:,1),X_sol_n(:,2),X_sol_n(:,3),'r','linewidth',1.5)
PlotBoi3('$X_n$','$Y_n$','$Z_n$',20,'LaTex')
axis equal

plot3(rLPs_n(LP,1), rLPs_n(LP,2), rLPs_n(LP,3), '^', 'markersize', 8, 'markerfacecolor', colors.std.ltgrn, 'markeredgecolor', colors.std.blue)

% plotSecondary(secondary)

view(0,90)

%%% Print Solution
fprintf('Iterations: %1.0d\n',iteration_counter)
fprintf('Error:      %1.1e\n\n',constraint_error)
prettyColVec([x0; 0; 0; 0; F_new(1); 0; F_new(end)])




% ========================================================================
%%% Closeout
% ========================================================================
toc(ticWhole);















