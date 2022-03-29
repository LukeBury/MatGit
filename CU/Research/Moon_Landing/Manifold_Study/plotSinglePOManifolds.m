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
plot_POs_alone = 0;

plot_unstable_manifolds = true;
plot_stable_manifolds   = true;

% ========================================================================
%%% Setup
% ========================================================================


PO_0 = [0.9880134301736027;
 0.0000000000000000;
 -0.0168339215530371;
 0.0000000000000000;
 -0.0040733760074276;
 0.0000000000000000;
 4.2714646720470677];



%%% Manifold Colors
% color_stable   = colors.red;
% color_unstable = colors.grn;

color_stable   = colors.blue2;
color_unstable = colors.red2;

% -------------------------------------------------
%%% Options
% -------------------------------------------------
%%% Number of manifolds per PO
n_nodes = 60;

%%% Scale the perturbation of the manifold node in the unstable direction
pertScale = 1e-5;

%%% Set propagation time for unstable manifolds
% Tf_manifolds_n = 0.7*pi;
Tf_manifolds_n = pi*4;
% Tf_manifolds_n = 24*pi;

%%% 3B System
famName = 'Jupiter_Europa';

% -------------------------------------------------
%%% Set up System 
% -------------------------------------------------
%%% Set primary & secondary
[primary, secondary] = assignPrimaryAndSecondary_CR3BP(famName, bodies);


%%% Normalizing constants
% rNorm:  n <-> km
% tNorm:  n <-> sec
% vNorm:  n <-> km/sec
[rNorm, tNorm, vNorm] = cr3bp_norms(primary, secondary, bodies.constants.G);


%%% Setting parameter structure
prms.u     = secondary.MR;
prms.rNorm = rNorm;
prms.R1    = primary.R / rNorm;
prms.R2    = secondary.R_n;
prms.n     = 1;

if contains(famName,'.CR3BP_J2pJ4pJ6pJ2s.')
    prms.J2p = primary.J2; prms.J4p = primary.J4; prms.J6p = primary.J6; prms.J2s = secondary.J2;
    warning('Add something for prms.n')
    return
end

%%% Equillibrium Points
rLPs_n = collinearEquilibriumPoints_ZH(prms);

% -------------------------------------------------
%%% Integrator options
% -------------------------------------------------
tol                      = 1e-13;
options                  = odeset('RelTol',tol,'AbsTol',tol);
options_XZStop           = odeset('Event',@event_yEqualsZeroPastL2,'RelTol',tol,'AbsTol',tol);
options_impactOrL1Escape = odeset('Event',@event_ImpactorL1Escape_CR3Bn,'RelTol',tol,'AbsTol',tol);
options_impact           = odeset('Event',@event_Impact_CR3Bn,'RelTol',tol,'AbsTol',tol);

% ========================================================================
%%% Find the manifolds and propagate them
% ========================================================================
% -------------------------------------------------
%%% Setup
% -------------------------------------------------
%%% STM vector
stm0        = eye(6);
stm0_colVec = reshape(stm0,36,1);

% -------------------------------------------------
%%% Loop through chosen indicies and compute manifolds at each
% -------------------------------------------------
% --------------------------
%%% Propagate initial conditions to get manifold node points
% --------------------------
%%% Integrate to get manifold nodes
[T_nodes, XSTM_nodes] = ode113(@Int_CR3BnSTM, linspace(0, PO_0(7), n_nodes+1), [PO_0(1:6); stm0_colVec], options, prms);

%%% Get the unstable eigenvector at the first node (initial condition)
monodromy_N0 = reshape(XSTM_nodes(end,7:42),6,6);
[eVec_stable_N0, eVec_unstable_N0] = getStableAndUnstableEigenvectors(monodromy_N0);

%%% Use state transition matrix to get stable and unstable eigenvectors 
%%% at each node
unstableEigenvectors = NaN(6, n_nodes);
stableEigenvectors   = NaN(6,n_nodes);
for node_i = 1:n_nodes
    %%% Grab STM from T0 to current node
    stm_Ni_T0                      = reshape(XSTM_nodes(node_i, 7:42),6,6);

    %%% Use STM to propagate the unstable eigenvector from N0 to the
    %%% current node
    unstableEigenvectors(:,node_i) = stm_Ni_T0 * eVec_unstable_N0;
    stableEigenvectors(:,node_i)   = stm_Ni_T0 * eVec_stable_N0;

    %%% Turn unstable eigenvector into a unit vector
    unstableEigenvectors(:,node_i) = unstableEigenvectors(:,node_i) ./ norm(unstableEigenvectors(:,node_i));
    stableEigenvectors(:,node_i)   = stableEigenvectors(:,node_i) ./ norm(stableEigenvectors(:,node_i));
end 


%%% Get rid of the repeat final/first point
X_nodes = XSTM_nodes(1:end-1,1:6);
T_nodes = T_nodes(1:end-1);

% --------------------------
%%% Preallocating for parfor
% --------------------------
man_unstable_p_POi = cell(n_nodes,1);
man_unstable_m_POi = cell(n_nodes,1);
man_stable_p_POi   = cell(n_nodes,1);
man_stable_m_POi   = cell(n_nodes,1);

% --------------------------
%%% Looping through manifolds ICs in parallel, computing the manifold
%%% direction, and propagating forward in time
% --------------------------
parfor node_i = 1:n_nodes
    % --------------------------
    %%% Pick node and integrate PO with STM from it
    % --------------------------        
    %%% Set IC of current manifold node
    X0_manNode_i_n = X_nodes(node_i,:)';

    %%% Create the perturbed initial condition and integrate the manifold
    X0_man_unstable_p_i = X0_manNode_i_n + unstableEigenvectors(:,node_i).*pertScale;
    X0_man_unstable_m_i = X0_manNode_i_n - unstableEigenvectors(:,node_i).*pertScale;
    X0_man_stable_p_i   = X0_manNode_i_n + stableEigenvectors(:,node_i).*pertScale;
    X0_man_stable_m_i   = X0_manNode_i_n - stableEigenvectors(:,node_i).*pertScale;

    %%% Integrate the manifolds
    [~, X_man_unstable_p] = ode113(@Int_CR3Bn, [0, Tf_manifolds_n], X0_man_unstable_p_i, options, prms);
    [~, X_man_unstable_m] = ode113(@Int_CR3Bn, [0, Tf_manifolds_n], X0_man_unstable_m_i, options, prms);

    [~, X_man_stable_p] = ode113(@Int_CR3Bn, [Tf_manifolds_n, 0], X0_man_stable_p_i, options, prms);
    [~, X_man_stable_m] = ode113(@Int_CR3Bn, [Tf_manifolds_n, 0], X0_man_stable_m_i, options, prms);
    % --------------------------
    %%% Store manifolds for this PO
    % --------------------------
    man_unstable_p_POi{node_i} = X_man_unstable_p(:,1:3);
    man_unstable_m_POi{node_i} = X_man_unstable_m(:,1:3);
    man_stable_p_POi{node_i}   = X_man_stable_p(:,1:3);
    man_stable_m_POi{node_i}   = X_man_stable_m(:,1:3);
end


    
    
    
    
    
    
figure; hold all
axis equal
plotSecondary(secondary)
PlotBoi3_CR3Bn(26)
view(0,90)
for man_i = 1:n_nodes
    if plot_unstable_manifolds
        p_u = plot3(man_unstable_p_POi{man_i}(:,1), man_unstable_p_POi{man_i}(:,2), man_unstable_p_POi{man_i}(:,3),'color', color_unstable,'linewidth',1);
        plot3(man_unstable_m_POi{man_i}(:,1), man_unstable_m_POi{man_i}(:,2), man_unstable_m_POi{man_i}(:,3),'color', color_unstable,'linewidth',1);
    end
    if plot_stable_manifolds
        p_s = plot3(man_stable_p_POi{man_i}(:,1), man_stable_p_POi{man_i}(:,2), man_stable_p_POi{man_i}(:,3),'color', color_stable,'linewidth',1);
        plot3(man_stable_m_POi{man_i}(:,1), man_stable_m_POi{man_i}(:,2), man_stable_m_POi{man_i}(:,3),'color', color_stable,'linewidth',1);
    end
end

if (plot_unstable_manifolds) && (~plot_stable_manifolds)
    legend(p_u, 'Unstable manifold','FontSize',14, 'location', 'best')
elseif (~plot_unstable_manifolds) && (plot_stable_manifolds)
    legend(p_s, 'Stable manifold','FontSize',14, 'location', 'best')
elseif (plot_unstable_manifolds) && (plot_stable_manifolds)
    legend([p_s, p_u], 'Stable manifold', 'Unstable manifold','FontSize',14, 'location', 'best')
end























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
















