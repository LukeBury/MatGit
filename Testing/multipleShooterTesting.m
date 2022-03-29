% ========================================================================
%%% Description
% ========================================================================
% Testground for using a multiple shooter for state continuity

% Created: 02/03/20
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

bodies = getBodyData(mbinPath);

PO_ICs = get_PO_ICs();
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
%%% Integration Options and parameters
% -------------------------------------------------
tol = 1e-13;
options = odeset('RelTol',tol,'AbsTol',tol);

% prms.u    = mu;
% prms.R2_n = R2_n;

prms.u = 0.0121505856;
X0_guess_n = [0.8233977783758328; 0; 0; 0; 0.1262548581832116; 0];
Tp_guess_n = 2.7429956388426318;
% ------------------------------------------------- 
%%% Initial conditions
% -------------------------------------------------
% --------------------------
% Choose the X0 and Tp guess
% --------------------------
% X0_guess_n = [rLPs_n(2,1) - 0.003; 0; 0; 0; 0.018; 0];
% % X0_guess_n = [rLPs_n(1,1) - 0.003; 0; 0; 0; 0.023; 0];
% Tp_guess_n = 3.2;
        

% % X0_guess_n = [1.180888549172871; -0.000000000000002; 0.000137978545959; 0.000000000000128; -0.155846477971307; 0];
% X0_guess_n = [1.181; -0.000000000000002; 0.00014; 0.000000000000128; -0.156; 0];

% % X0_guess_n = [1.180888549172871; -0.000000000000002; -0.000137978545959; 0.000000000000128; -0.155846477971307; 0];
% X0_guess_n = [1.180888; -0.000000000000002; -0.000137978545959; 0.000000000000128; -0.155846477971307; 0];
% Tp_guess_n = 3.415510771287679;

% ------------------------------------------------- 
%%% Other options
% -------------------------------------------------
error_tol = 1e-12; 
% error_tol = 1e-10;
% error_tol = 1e-9;

%%% Number of nodes for multiple shooter
n_Nodes = 3; % 3

%%% Maximum number of multiple-shooter iterations for any given family
%%% member before kicking out of the loop and adjusting the tuning
%%% parameters
iterMax = 350;

%%% Specify a stepsize if you're a risk taker (should be between 0 and 1)
stepsize = 1;

% ------------------------------------------------- 
%%% Integrate first guess to get reference trajectory and to create nodes
%%% from it
% -------------------------------------------------
%%% Initialize STM
stm0        = eye(6);
stm0_colVec = reshape(stm0,36,1);

%%% Integrate reference trajectory (initial guess for full period)
[Tref_n, Xref_n] = ode113(@Int_CR3Bn, [0, Tp_guess_n], X0_guess_n, options, prms);

%%% Cut up the reference trajectory into N nodes
[T_nodes, X_nodes] = ode113(@Int_CR3BnSTM, linspace(0, Tp_guess_n,n_Nodes+1), [X0_guess_n; stm0_colVec], options, prms);

%%% Setting initial nodes and node-time
X_nodes = X_nodes(1:n_Nodes,1:6);
T_PO_i  = Tp_guess_n;

% ------------------------------------------------- 
%%% Create the first free-variable column vector
% -------------------------------------------------
F_new   = X_nodes';
F_new   = F_new(:);
F_new   = [F_new; T_PO_i];

% ------------------------------------------------- 
%%% Initializing things
% -------------------------------------------------
iteration_counter = 0;
constraint_error = 1000000;

% ========================================================================
%%% Multiple shooter
% ========================================================================
while (constraint_error > error_tol) && (iteration_counter < iterMax)

    %%% Count iteration
    iteration_counter = iteration_counter + 1;
    
    % --------------------------
    % Using multiple shooting - loop through nodes, integrate, and 
    % populate constraint vector and create DF matrix
    % --------------------------
    [DF_mat, constraints] = multShooter_stateContinuity(n_Nodes, F_new, @Int_CR3BnSTM, options, prms);

    %%% Compute error of constraint vector for current F
    constraint_error = norm(constraints);
    
    F_new
    %%% Compute new free-variable vector if error is not converged    
    if (constraint_error > error_tol)
        warning('off','MATLAB:nearlySingularMatrix')
        F_new = F_new - stepsize* DF_mat'*((DF_mat*(DF_mat'))\constraints);
        warning('on','MATLAB:nearlySingularMatrix')
    end
    
    fprintf('Iteration: %1d\n', iteration_counter)
    fprintf('Error:     %1.2e\n',constraint_error)
    fprintf('============================\n')
end % while constraint_error > error_tol

if (iteration_counter >= iterMax)
    warning('Multiple shooter didn''t converge')
    return
end

%%% Integrate converged solution
[Tsol_n, Xsol_n] = ode113(@Int_CR3Bn, [0, F_new(end)], F_new(1:6), options, prms);

%%% Plot reference (guess) trajectory 
figure; hold all
plot3(Xref_n(:,1),Xref_n(:,2),Xref_n(:,3))
PlotBoi3_CR3Bn(20)
% plotSecondary(secondary)

%%% Plot converged solution
figure; hold all
plot3(Xsol_n(:,1),Xsol_n(:,2),Xsol_n(:,3))
PlotBoi3_CR3Bn(20)
% plotSecondary(secondary)







% ========================================================================
%%% Closeout
% ========================================================================
tocWhole = toc(ticWhole);
fprintf('Elapsed time: %1.4f seconds\n',tocWhole)
















