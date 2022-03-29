% ========================================================================
%%% Description
% ========================================================================
% After having identified a PO from which a new family bifurcates, this
% script will help identify the new family

% Created: 12/30/20
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


% ========================================================================
%%% Setup
% ========================================================================
% -------------------------------------------------
%%% Inputs / ICs
% -------------------------------------------------
%%% Primary_Secondary system
% systemName = 'Jupiter_Europa';
% systemName = 'Saturn_Titan';
systemName = 'Saturn_Enceladus';
% 989
% systemName = 'Earth_Moon';




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% 
% 
% 3P4 - 568
PO_0 = [0.9990984652798989;
 0.0000000000000000;
 -0.0027774454190227;
 -0.0000000000000019;
 -0.0073273546776289;
 0.0000000000000011;
 6.0800228609930249];

% -------------------------------------------------
%%% Bifurcation guess setup
% -------------------------------------------------
%%% Scaling for bump-vector
bumpVec_scalar =1e-3;

% -------------------------------------------------
%%% Shooter setup
% -------------------------------------------------
%%% Number of nodes to cut the initial guess into. Choose the one that
%%% creates the smallest s_important
% n_Nodes = 1; 
% n_Nodes = 2;
% n_Nodes = 3;
% n_Nodes = 4;
% n_Nodes = 5;
% n_Nodes = 6; 
% n_Nodes = 7; 
n_Nodes = 8; 

% -------------------------------------------------
%%% Parameter setup
% -------------------------------------------------
% --------------------------
%%% Set primary & secondary
% --------------------------
[primary, secondary] = assignPrimaryAndSecondary_CR3BP(systemName, bodies);

% --------------------------
%%% System
% --------------------------
%%% Normalizing constants
% rNorm:  n <-> km
% tNorm:  n <-> sec
% vNorm:  n <-> km/sec
[rNorm, tNorm, vNorm] = cr3bp_norms(primary, secondary, bodies.constants.G);

%%% Setting parameters
prms.u  = secondary.MR;
prms.n  = 1;
prms.R2 = secondary.R_n;


% ------------------------------------------------- 
%%% Integration Options
% -------------------------------------------------
tol = 1e-13;
options = odeset('RelTol',tol,'AbsTol',tol);

%%% Initialize STM as column vector to be added with state
stm0_colVec = reshape(eye(6),36,1);

% ========================================================================
%%% Find the tangent PO
% ========================================================================

% % 
% %%% Integrate the POs on either side of the detected bifurcation
[T_PO, X_PO] = ode113(@Int_CR3BnSTM, linspace(0, PO_0(7), 10000), [PO_0(1:6); stm0_colVec], options, prms);
% [T_PO2, X_PO2] = ode113(@Int_CR3BnSTM, linspace(0, otherPO2(7), 10000), [otherPO2(1:6); stm0_colVec], options, prms);
% % [T_PO3, X_PO3] = ode113(@Int_CR3BnSTM, linspace(0, otherPO3(7), 10000), [otherPO3(1:6); stm0_colVec], options, prms);
% 
% %%% Find alpha and beta parameters for each PO
% stm_tf_t0                           = reshape(X_PO(end,7:42),6,6);
% monodromy                           = stm_tf_t0;
% [eigenVectors_new, eigenValues_new] = eig(monodromy);
% [alpha, beta] = getBrouckeStabilityParameters(diag(eigenValues_new), monodromy);
% 
% stm_tf_t02                           = reshape(X_PO2(end,7:42),6,6);
% monodromy2                           = stm_tf_t02;
% [eigenVectors_new2, eigenValues_new2] = eig(monodromy2);
% [alpha2, beta2] = getBrouckeStabilityParameters(diag(eigenValues_new2), monodromy2);
% 
% % stm_tf_t03                           = reshape(X_PO3(end,7:42),6,6);
% % monodromy3                           = stm_tf_t03;
% % [eigenVectors_new3, eigenValues_new3] = eig(monodromy3);
% % [alpha3, beta3] = getBrouckeStabilityParameters(diag(eigenValues_new3), monodromy3);
% % 
% %%% Check for tangent bifurcation
% % plot_BrouckeStabilityDiagram([alpha; alpha2; alpha3], [beta; beta2; beta3],true)
% plot_BrouckeStabilityDiagram([alpha; alpha2], [beta; beta2], false)

% ========================================================================
%%% Find the tangent PO
% ========================================================================

% --------------------------
%%% Proceed to identify new guess
% --------------------------
%%% Create nodes along the period-multiplied PO
[~, X_nodes] = get_nodes([PO_0(1:6); stm0_colVec], [0, PO_0(7)], n_Nodes+1, @Int_CR3BnSTM, options, prms);
X_nodes = X_nodes(1:n_Nodes,1:6);

%%% If y0=0, use the correct multiple-shooting process for keeping it that
%%% way
if PO_0(2)==0
    %%% Build the free-variable vector
    F_new   = reshape(X_nodes',6*n_Nodes,1);
    F_new   = [F_new; PO_0(7)];
    F_new   = [F_new(1); F_new(3:end)];

    %%% Based on free-variable vector, build constraints vector and DF
    %%% matrix
    [DF_mat, constraints] = multShooter_stateContinuity_y0Fixed(n_Nodes, F_new, @Int_CR3BnSTM, options, prms, 0);
end

%%% Use single-value decomposition to find tangent direction of new family
s = svd(DF_mat); % s is equal to diag(S)
s_important = s(end-1)
[U,S,V] = svd(DF_mat);
V_vec = V(:,n_Nodes*6-1);
bumpVec = [V_vec(1); 0; V_vec(2:5); V_vec(end)]
if PO_0(2) ~= 0
    warning('Problem here')
end

%%% Add the bump to the initial PO to try and find the tangent family
newPO_guess = PO_0(1:7) + bumpVec.*bumpVec_scalar;

%%% Integrate the new guess
[T_newPO, X_newPO] = ode113(@Int_CR3BnSTM, linspace(0, newPO_guess(7), 10000), [newPO_guess(1:6); stm0_colVec], options, prms);

%%% Plot the results
figure; hold all
p1 = plot3(X_PO(:,1),X_PO(:,2),X_PO(:,3),'m');
p2 = plot3(X_newPO(:,1),X_newPO(:,2),X_newPO(:,3),'b');
PlotBoi3_CR3Bn(26)
legend([p1, p2], 'Old PO', 'New PO Guess')

missDistance_norm = norm(X_newPO(end,1:3) - X_newPO(1,1:3))

prettyColVec(newPO_guess)


%%% Find stability info
stm_tf_t0                           = reshape(X_newPO(end,7:42),6,6);
monodromy                           = stm_tf_t0;
[eigenVectors_new, eigenValues_new] = eig(monodromy);
[S1, S2]                            = getStabilityIndices(diag(eigenValues_new));
fprintf('Stability Indices: [%1.1f, %1.1f]\n\n', S1, S2)





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
fprintf('\nElapsed time: %1.4f seconds\n',tocWhole)
















