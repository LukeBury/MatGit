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
plot_POs_alone = 1;

plot_unstable_manifolds = 1;
plot_stable_manifolds   = 1;

% ========================================================================
%%% Setup
% ========================================================================
% -------------------------------------------------
%%% Options
% -------------------------------------------------
%%% Number of manifolds per PO
n_nodes = 50;

%%% Scale the perturbation of the manifold node in the unstable direction
pertScale = 1e-8;

%%% Set propagation time for unstable manifolds
% Tf_manifolds_n = 4*pi;
989
Tf_manifolds_n = 3*pi;

%%% 3B System
famName_bodies = 'Earth_Moon';
% famName_bodies = 'Jupiter_Europa';
% famName_bodies = 'Jupiter_Ganymede';
% famName_bodies = 'Saturn_Enceladus';
% famName_bodies = 'Neptune_Triton';

%%% PO Family
famName_PO_Family = 'L1_Lyapunov';
% famName_PO_Family = 'L1_Vertical';
% famName_PO_Family = 'L1_SHalo';
% famName_PO_Family = 'L2_Lyapunov';
% famName_PO_Family = 'L2_Vertical';
% famName_PO_Family = 'L2_SHalo';

% -------------------------------------------------
%%% Set up System 
% -------------------------------------------------
%%% Create family name
famName = [famName_bodies, '.CR3BP.', famName_PO_Family];

%%% Set primary & secondary
[primary, secondary] = assignPrimaryAndSecondary_CR3BP(famName, bodies);

%%% Normalizing constants
rNorm = secondary.a;         % n <-> km
tNorm = 1/secondary.meanMot; % n <-> sec
vNorm = rNorm / tNorm;       % n <-> km/sec

%%% Setting parameter structure
prms.u     = secondary.MR;
prms.rNorm = rNorm;
prms.R1    = primary.R / rNorm;
prms.R2    = secondary.R_n;

if contains(famName,'.CR3BP_J2pJ4pJ6pJ2s.')
    prms.J2p = primary.J2; prms.J4p = primary.J4; prms.J6p = primary.J6; prms.J2s = secondary.J2;
end

%%% Equillibrium Points
rLPs_n = collinearEquilibriumPoints_ZH(prms);

% -------------------------------------------------
%%% Load PO family data
% -------------------------------------------------
%%% Load PO family
FamilyDataFile = [POfamilyPath, famName, '.txt'];
POFamilyData = dlmread(FamilyDataFile,',',1,0);

%%% Column specifiers
c_x0_n                = 1;
c_y0_n                = 2;
c_z0_n                = 3;
c_xd0_n               = 4;
c_yd0_n               = 5;
c_zd0_n               = 6;
c_Tp_n                = 7;
c_JC                  = 8;
c_L2FlythroughVel_mps = 9;
c_landingVelocity_mps = 10;
c_stabilityIndex1     = 11;
c_stabilityIndex2     = 12;

% -------------------------------------------------
%%% Integrator options
% -------------------------------------------------
tol                      = 1e-13;
options                  = odeset('RelTol',tol,'AbsTol',tol);
options_XZStop           = odeset('Event',@event_yEqualsZeroPastL2,'RelTol',tol,'AbsTol',tol);
options_impactOrL1Escape = odeset('Event',@event_ImpactorL1Escape_CR3Bn,'RelTol',tol,'AbsTol',tol);
options_impact           = odeset('Event',@event_Impact_CR3Bn,'RelTol',tol,'AbsTol',tol);

% -------------------------------------------------
%%% Set PO family indices to study manifolds of
% -------------------------------------------------
if     isequal(famName, 'Earth_Moon.CR3BP.L1_Lyapunov') % max - 507
%     PO_indicies = [27, 35, 52, 70, 95, 120, 150, 300, 500];
    PO_indicies = [180];
elseif isequal(famName, 'Earth_Moon.CR3BP.L1_Vertical') % max - 512
    %%% don't see impacts until vertical orbits are 0.1415 tall, or 54,440 km
    %%% 160 is just a bunch of close misses
    %%% ~183 is first hit
    %%% ~304 is last hit
    PO_indicies = [183, 187, 190, 195, 200, 300, 302, 304];
%     PO_indicies = [300];    
elseif isequal(famName, 'Earth_Moon.CR3BP.L1_SHalo') % max - 513
    %%% 46 has an example of the double impact that I have to handle for
    %%% 230 has some cool implications for earth
    %%% ~305 is last impact (not really, but effectively ..some around 370)
    PO_indicies = [2, 7, 18, 27, 40, 42, 48, 70, 110, 155, 170, 210, 230, 300];
%     PO_indicies = [490];
elseif isequal(famName, 'Jupiter_Europa.CR3BP.L2_Lyapunov')
%     PO_indicies = [1, 5, 15, 17, 19, 21, 24, 29, 35, 42, 49, 59, 75, 100, 150, 220, 340, 510];
    PO_indicies = [19];
elseif isequal(famName, 'Jupiter_Europa.CR3BP.L2_Vertical') % max - 511
%     PO_indicies = [2, 3, 5, 10, 14, 15, 17, 20, 25, 30, 50, 100]; % interesting things between 15 and 18 .. the manifold first grazes off the moon
%     *** Should be integrating for 6 pi to get extra set of intersections
%     *** Around 100 can be used for shallow landings at poles
    PO_indicies = [101];
elseif isequal(famName, 'Jupiter_Europa.CR3BP.L2_SHalo') % max - 525
%     PO_indicies = [5, 10, 15, 20, 25, 30, 35, 45, 60, 100];
%     PO_indicies = [60, 70, 80, 90, 95, 100, 105, 110];
    PO_indicies = [90, 91, 92, 93, 94, 95];
    
elseif isequal(famName, 'Jupiter_Ganymede.CR3BP.L2_Lyapunov') % max - 583
    %%% around 350 has large regions of shallow landing angles (4 different areas)
    PO_indicies = [1, 74, 79, 85, 86, 92, 100, 110, 150, 220, 350, 500];
%     PO_indicies = [550];
elseif isequal(famName, 'Jupiter_Ganymede.CR3BP.L2_Vertical') % max - 513
    %%% ~90 has great polar landers ... so does 109
    PO_indicies = [2, 10, 17, 18, 19, 28, 30, 54, 57, 64, 70, 77, 90, 109];
%     PO_indicies = [450];
elseif isequal(famName, 'Jupiter_Ganymede.CR3BP.L2_SHalo') % max - 541
    PO_indicies = [30, 45, 55, 65, 75, 85, 95, 110, 120, 135, 150, 185];
%     PO_indicies = [250];
elseif isequal(famName, 'Saturn_Enceladus.CR3BP.L2_Lyapunov')
%     PO_indicies = [1,];
    PO_indicies = [1];
elseif isequal(famName, 'Saturn_Enceladus.CR3BP.L2_Vertical')
%     PO_indicies = [2, 8, 15, 25, 45, 60, 61, 63, 65, 67, 70, 80, 83];
    %%% Around 80 is good for polar landings
    %%% Around 83, polar landing solutions pass over the landing site once
    %%%     before landing!
%     PO_indicies = [63, 65, 70, 80, 81, 82, 83];
    PO_indicies = [81];
elseif isequal(famName, 'Saturn_Enceladus.CR3BP.L2_SHalo')
    PO_indicies = [1, 10];
    
elseif isequal(famName, 'Neptune_Triton.CR3BP.L2_Lyapunov')
    PO_indicies = [1, 10];
elseif isequal(famName, 'Neptune_Triton.CR3BP.L2_Vertical')
    PO_indicies = [1, 10];
elseif isequal(famName, 'Neptune_Triton.CR3BP.L2_SHalo')
    PO_indicies = [60];
end

% ========================================================================
%%% (optional) Plot the POs whos indices were chosen
% ========================================================================
if plot_POs_alone
    p = figure; hold all
    p.Position = [112 385 560 420];
    axis equal
    plotSecondary(secondary)
    PlotBoi3_CR3Bn(20)
    view(0,90)
    
    for PO_index = PO_indicies
        %%% Grab initial conditions
        X0_PO_i_n = POFamilyData(PO_index, c_x0_n:c_zd0_n)';
        Tp_PO_i_n = POFamilyData(PO_index, c_Tp_n);
        
        %%% Integrate
        [~, X_ref_PO] = ode113(@Int_CR3Bn, linspace(0, Tp_PO_i_n*2, 100000), X0_PO_i_n, options, prms);
        
        %%% Plot
        plot3(X_ref_PO(:,1),X_ref_PO(:,2),X_ref_PO(:,3),'b','linewidth',1.5)
    end
    
end % plot_POs_from_chosen_indices

% ========================================================================
%%% Loop through chosen PO indicies, find the unstable manifolds, and
%%% propagate them
% ========================================================================
% -------------------------------------------------
%%% Pre-Loop setup
% -------------------------------------------------
%%% STM vector
stm0        = eye(6);
stm0_colVec = reshape(stm0,36,1);

%%% Preallocate outputs
manifolds_unstable_p_n = cell(length(PO_indicies),1);
manifolds_unstable_m_n = cell(length(PO_indicies),1);
manifolds_stable_p_n   = cell(length(PO_indicies),1);
manifolds_stable_m_n   = cell(length(PO_indicies),1);

% -------------------------------------------------
%%% Loop through chosen indicies and compute manifolds at each
% -------------------------------------------------
PO_i = 0;
for PO_index = PO_indicies
    PO_i = PO_i + 1;
    % --------------------------
    %%% Propagate initial conditions to get manifold node points
    % --------------------------
    %%% Grab IC
    X0_PO_i_n = POFamilyData(PO_index, c_x0_n:c_zd0_n)';
    Tp_PO_i_n = POFamilyData(PO_index, c_Tp_n);
    
    %%% Integrate to get manifold nodes
    [T_nodes, XSTM_nodes] = ode113(@Int_CR3BnSTM, linspace(0, Tp_PO_i_n, n_nodes+1), [X0_PO_i_n; stm0_colVec], options, prms);
    
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
    
    % --------------------------
    %%% Store manifolds for this PO into all-PO structure
    % --------------------------
    manifolds_unstable_p_n{PO_i} = man_unstable_p_POi;
    manifolds_unstable_m_n{PO_i} = man_unstable_m_POi;
    manifolds_stable_p_n{PO_i}   = man_stable_p_POi;
    manifolds_stable_m_n{PO_i}   = man_stable_m_POi;
end


for PO_i = 1:length(PO_indicies)
    
    figure; hold all
    axis equal
    plotSecondary(secondary)
    PlotBoi3_CR3Bn(20)
    view(0,90)
    for man_i = 1:n_nodes
        if plot_unstable_manifolds
            plot3(manifolds_unstable_p_n{PO_i}{man_i}(:,1), manifolds_unstable_p_n{PO_i}{man_i}(:,2), manifolds_unstable_p_n{PO_i}{man_i}(:,3),'r','linewidth',1)
            plot3(manifolds_unstable_m_n{PO_i}{man_i}(:,1), manifolds_unstable_m_n{PO_i}{man_i}(:,2), manifolds_unstable_m_n{PO_i}{man_i}(:,3),'r','linewidth',1)
        end
        if plot_stable_manifolds
            plot3(manifolds_stable_p_n{PO_i}{man_i}(:,1), manifolds_stable_p_n{PO_i}{man_i}(:,2), manifolds_stable_p_n{PO_i}{man_i}(:,3),'g','linewidth',1)
            plot3(manifolds_stable_m_n{PO_i}{man_i}(:,1), manifolds_stable_m_n{PO_i}{man_i}(:,2), manifolds_stable_m_n{PO_i}{man_i}(:,3),'g','linewidth',1)
        end
    end
    
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
















