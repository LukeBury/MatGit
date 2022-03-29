% ========================================================================
%%% Description
% ========================================================================
% For Ploltting the stable and unstable manifolds of periodic orbits in the
% CR3BP

% Created: 11/04/19
% Author : Luke Bury, luke.bury@colorado.edu
% ========================================================================
%%% Initialization
% ========================================================================
clear
clc
close all
mbinPath = '~/CU_Google_Drive/Documents/MatGit/mbin';
POfamilyPath = '~/CU_Google_Drive/Documents/MatGit/mbin/Data/InitialConditions/PO_Families/';
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

% ========================================================================
%%% Setup
% ========================================================================
% -------------------------------------------------
%%% Options
% -------------------------------------------------
%%% Number of manifolds per PO
n_nodes = 5;

%%% Scale the perturbation of the manifold node in the unstable direction
pertScale = 1e-8;

%%% Set propagation time for unstable manifolds
Tf_manifolds_n = 4*pi;



% -------------------------------------------------
%%% Integrator options
% -------------------------------------------------
tol                      = 1e-13;
options                  = odeset('RelTol',tol,'AbsTol',tol);



% -------------------------------------------------
%%% Loop through chosen indicies and compute manifolds at each
% -------------------------------------------------




% --------------------------
%%% Propagate initial conditions to get manifold node points
% --------------------------
%%% Grab IC
% X0_PO_i_n = [1.04574850772215;-6.75971996578334e-20;-0.194717827575194;1.00478436335094e-16;-0.149445161640467;4.36221967835965e-20];
% Tp_PO_i_n = 1.83212355593254;
% prms.u = 1.215060379322131e-02;

X0_PO_i_n = [1.0467099745706516,-0.0000000000000011,-0.0009488120693468,-0.0000000000240595,-0.0374261978743848,0.0000000000005953]';
Tp_PO_i_n = 3.1674713005139745;
prms.u = bodies.triton.MR;


[~, X_ref] = ode113(@Int_CR3Bn, [0, Tp_PO_i_n], X0_PO_i_n, options, prms);

% %%% Integrate to get manifold nodes
% [T_nodes, X_nodes] = ode113(@Int_CR3Bn, linspace(0, Tp_PO_i_n, n_nodes+1), X0_PO_i_n, options, prms);
% 
% %%% Get rid of the repeat final/first point
% X_nodes = X_nodes(1:end-1,:);
% T_nodes = T_nodes(1:end-1);

% --------------------------
%%% Get initial monodromy
% --------------------------
%%% STM vector
stm0        = eye(6);
stm0_colVec = reshape(stm0,36,1);
[T_nodes, XSTM_nodes] = ode113(@Int_CR3BnSTM, linspace(0, Tp_PO_i_n, n_nodes+1), [X0_PO_i_n; stm0_colVec], options, prms);

stm_tf_t0 = reshape(XSTM_nodes(end,7:42),6,6);
monodromy_N0 = stm_tf_t0;
[eVec_stable_N0, eVec_unstable_N0] = getStableAndUnstableEigenvectors(monodromy_N0);
    

X_N3 = XSTM_nodes(3,:);
T_N3 = T_nodes(3);

stm_N3_t0 = reshape(XSTM_nodes(3,7:42),6,6);

eVec_unstable_N3_from_N0 = stm_N3_t0 * eVec_unstable_N0;


%%% compute monodromy from N3
[~, XSTM_N3] = ode113(@Int_CR3BnSTM, [0, Tp_PO_i_n], [X_N3(1,1:6)'; stm0_colVec], options, prms);

stm_tf_t0_N3 = reshape(XSTM_N3(end,7:42),6,6);
monodromy_N3 = stm_tf_t0_N3;
[eVec_stable_N3, eVec_unstable_N3] = getStableAndUnstableEigenvectors(monodromy_N3);

%%% Normalize and compare
eVec_unstable_N3_from_N0 = eVec_unstable_N3_from_N0./norm(eVec_unstable_N3_from_N0);
eVec_unstable_N3 = eVec_unstable_N3./norm(eVec_unstable_N3);

eVec_unstable_N3_from_N0 - eVec_unstable_N3

% for manifold_i = 1:n_nodes
%     % --------------------------
%     %%% Pick node and integrate PO with STM from it
%     % --------------------------        
%     %%% Set IC of current manifold node
%     X0_manNode_i_n = X_nodes(manifold_i,:)';
% 
%     %%% Integrate manifold node to get monodromy matrix of PO at this
%     %%% node
%     [~, Xstm_man_n] = ode113(@Int_CR3BnSTM, [0, Tp_PO_i_n], [X0_manNode_i_n; stm0_colVec], options, prms);
% 
%     % --------------------------
%     %%% Perturb node X0 in directions of manifolds and integrate
%     % --------------------------  
%     %%% Get monodromy matrix and determine its stable and unstable
%     %%% eigenvectors
%     stm_tf_t0 = reshape(Xstm_man_n(end,7:42),6,6);
%     monodromy = stm_tf_t0;
%     [eVec_stable, eVec_unstable] = getStableAndUnstableEigenvectors(monodromy);
% 
%     if manifold_i == 1
%         firstUnstableEigenvector = eVec_unstable;
%     end
% 
%     if manifold_i == 30
%         [~, X30stm_man_n] = ode113(@Int_CR3BnSTM, [0, T_nodes(30)], [X_nodes(1,:)'; stm0_colVec], options, prms);
%          stm30_tf_t0 = reshape(X30stm_man_n(end,7:42),6,6);
%          thirtiethUnstableEigenvector = stm30_tf_t0 * firstUnstableEigenvector;
% 
%          eVec_unstable./norm(eVec_unstable) - thirtiethUnstableEigenvector./norm(thirtiethUnstableEigenvector)
%          989
% 
%     end
% 
%     %%% Create the perturbed initial condition and integrate the manifold
%     X0_man_unstable_p_i = X0_manNode_i_n + eVec_unstable.*pertScale;
%     X0_man_unstable_m_i = X0_manNode_i_n - eVec_unstable.*pertScale;
%     X0_man_stable_p_i   = X0_manNode_i_n + eVec_stable.*pertScale;
%     X0_man_stable_m_i   = X0_manNode_i_n - eVec_stable.*pertScale;
% 
%     %%% Integrate the manifolds
%     [~, X_man_unstable_p] = ode113(@Int_CR3Bn, [0, Tf_manifolds_n], X0_man_unstable_p_i, options, prms);
%     [~, X_man_unstable_m] = ode113(@Int_CR3Bn, [0, Tf_manifolds_n], X0_man_unstable_m_i, options, prms);
% 
%     [~, X_man_stable_p] = ode113(@Int_CR3Bn, [Tf_manifolds_n, 0], X0_man_stable_p_i, options, prms);
%     [~, X_man_stable_m] = ode113(@Int_CR3Bn, [Tf_manifolds_n, 0], X0_man_stable_m_i, options, prms);
%     % --------------------------
%     %%% Store manifolds for this PO
%     % --------------------------
%     PO_from_nodes{manifold_i}      = Xstm_man_n(:,1:3)
%     man_unstable_p_POi{manifold_i} = X_man_unstable_p(:,1:3);
%     man_unstable_m_POi{manifold_i} = X_man_unstable_m(:,1:3);
%     man_stable_p_POi{manifold_i}   = X_man_stable_p(:,1:3);
%     man_stable_m_POi{manifold_i}   = X_man_stable_m(:,1:3);
% end
%     
























% ========================================================================
%%% Formatting Structures
% ========================================================================
% -------------------------------------------------
%%% 
% -------------------------------------------------

% --------------------------
%%% 
% --------------------------

% ========================================================================
%%% Closeout
% ========================================================================
tocWhole = toc(ticWhole);
fprintf('\n(Elapsed time: %1.4f seconds)\n',tocWhole)
















