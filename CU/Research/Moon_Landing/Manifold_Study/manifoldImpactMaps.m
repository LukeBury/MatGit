% ========================================================================
%%% Description
% ========================================================================
% 

% Created: 10/22/19
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
plot_POs_from_chosen_indices = 1;

store_and_plot_impact_trajectories    = 1;
store_and_plot_manifold_trajectories  = 1;


% ========================================================================
%%% Setup
% ========================================================================
% -------------------------------------------------
%%% Options
% -------------------------------------------------
%%% Number of manifolds per PO
% n_nodes = 200;
989
n_nodes = 20;

%%% Scale the perturbation of the manifold node in the unstable direction
pertScale = 1e-8;

%%% Set propagation time for unstable manifolds
% Tf_manifolds_n = 4*pi;
989
Tf_manifolds_n = 2*pi;

%%% Bins and colors for impact angles
impactAngleBins_deg = [0, 3, 10, 20, 45, 90];

n_impactAngleBins = length(impactAngleBins_deg) - 1;
impactAngleColors = colorScale([colors.std.cyan; colors.std.mag], n_impactAngleBins);

%%% 3B System
% famName_bodies = 'Earth_Moon';
famName_bodies = 'Jupiter_Europa';
% famName_bodies = 'Jupiter_Ganymede';
% famName_bodies = 'Saturn_Enceladus';
% famName_bodies = 'Neptune_Triton';

%%% PO Family
% famName_PO_Family = 'L1_Lyapunov';
% famName_PO_Family = 'L1_Vertical';
% famName_PO_Family = 'L1_SHalo';
famName_PO_Family = 'L2_Lyapunov';
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
    PO_indicies = [27, 35, 52, 70, 95, 120, 150, 300, 500];
%     PO_indicies = [507];
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
    PO_indicies = [70, 500];
elseif isequal(famName, 'Jupiter_Europa.CR3BP.L2_Lyapunov')
    %%% Around 340, can impact 4 regions at shallow angle
%     PO_indicies = [1, 5, 15, 17, 19, 21, 24, 29, 35, 42, 49, 59, 75, 100, 150, 220, 340, 510];
    PO_indicies = [10];
elseif isequal(famName, 'Jupiter_Europa.CR3BP.L2_Vertical') % max - 511
%     PO_indicies = [2, 3, 5, 10, 14, 15, 17, 20, 25, 30, 50, 100]; % interesting things between 15 and 18 .. the manifold first grazes off the moon
%     *** Should be integrating for 6 pi to get extra set of intersections
%     *** Around 100 can be used for shallow landings at poles
    PO_indicies = [101];
elseif isequal(famName, 'Jupiter_Europa.CR3BP.L2_SHalo') % max - 525
    %%% Impacts appear to stop after 226
    PO_indicies = [5, 10, 15, 20, 25, 30, 35, 45, 60, 92, 100, 140, 155, 175, 210, 220, 226];
    PO_indicies = [45, 500];
    
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
elseif isequal(famName, 'Saturn_Enceladus.CR3BP.L2_Lyapunov') % max - 509
    %%% Around 80, we see some manifolds interact with L1 Lyapunovs 
    %%% PO intersects Enceladus past ~325
    PO_indicies = [1, 2, 13, 15, 16, 22, 29, 42, 160, 250, 320];
%     PO_indicies = [320];
elseif isequal(famName, 'Saturn_Enceladus.CR3BP.L2_Vertical') % max - 612
    %%% Around 80 is good for polar landings
    %%% Around 83, polar landing solutions pass over the landing site once
    %%%     before landing!-
    %%% Family impacts with the anti-saturn side of Enceladus by the end
%     PO_indicies = [2, 8, 15, 25, 45, 60, 61, 63, 65, 67, 70, 77, 80, 82, 83];
    PO_indicies = [600];
elseif isequal(famName, 'Saturn_Enceladus.CR3BP.L2_SHalo') % max - 578
    %%% Around 180/190 we get some results whos symmetric NHalo results might
    %%% be great for the tiger strips
%     PO_indicies = [20, 30, 45, 75, 90, 100, 110, 120, 130, 140, 150, 160, 180, 190, 210, 220, 240, 260, 270, 280];
    PO_indicies = [290];
elseif isequal(famName, 'Neptune_Triton.CR3BP.L2_Lyapunov') % max - 720
    %%% Around 570, we see the same behavior as with Ganymede and Europa 
    %%%     where 4 different regions are accessible via shallow impact
    PO_indicies = [1, 50, 85, 95, 140, 190, 210, 270, 360, 400, 430, 470, 570];
%     PO_indicies = [570];
elseif isequal(famName, 'Neptune_Triton.CR3BP.L2_Vertical') % max - 633
    %%% Stops impacting after 220
    PO_indicies = [25, 130, 137, 138, 160, 170, 175, 179, 187, 217, 220];
%     PO_indicies = [222];

elseif isequal(famName, 'Neptune_Triton.CR3BP.L2_SHalo') % max - 616
    %%% 215 has some manifolds that would maybe be great polar impacters
    %%%     with the addition of a burn
    %%% No more impacts after 320 (found thus far)
%     PO_indicies = [70, 90, 105, 135, 155, 175, 190, 202, 215, 240, 290, 320];
    PO_indicies = [215];
end

% ========================================================================
%%% Propagate reference POs and plot the POs whos indices were chosen
% ========================================================================
if plot_POs_from_chosen_indices
    p = figure; hold all
    p.Position = [112 385 560 420];
    axis equal
    plotSecondary(secondary)
    PlotBoi3_CR3Bn(20)
    view(0,0)
end % plot_POs_from_chosen_indices
    
%%% Preallocate for storing reference POs
X_POs = cell(length(PO_indicies), 1);

PO_i = 0;
for PO_index = PO_indicies
    PO_i = PO_i + 1;
    
    %%% Grab initial conditions
    X0_PO_i_n = POFamilyData(PO_index, c_x0_n:c_zd0_n)';
    Tp_PO_i_n = POFamilyData(PO_index, c_Tp_n);

    %%% Integrate
    [~, X_ref_n] = ode113(@Int_CR3Bn, [0, Tp_PO_i_n], X0_PO_i_n, options, prms);
    X_POs{PO_i} = X_ref_n(:,1:3);
    
    %%% Plot
    if plot_POs_from_chosen_indices
        plot3(X_ref_n(:,1),X_ref_n(:,2),X_ref_n(:,3),'b','linewidth',1.5)
    end
end

% ========================================================================
%%% Loop through chosen PO indicies, find the unstable manifolds, and
%%% propagate them until impact (or escape)
% ========================================================================
% -------------------------------------------------
%%% Pre-Loop setup
% -------------------------------------------------
%%% STM vector
stm0        = eye(6);
stm0_colVec = reshape(stm0,36,1);

%%% Preallocate outputs
X0s_n             = cell(length(PO_indicies),1);
latLons_deg       = cell(length(PO_indicies),1);
impactAngles_deg  = cell(length(PO_indicies),1);
latLonHeadingHats = cell(length(PO_indicies),1);
impactColors      = cell(length(PO_indicies),1);
Ximpacts_n        = cell(length(PO_indicies),1);
impactBins        = cell(length(PO_indicies),1);

if store_and_plot_manifold_trajectories
    trajectories_n = cell(length(PO_indicies),1);
end
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
    
    %%% Integrate to get manifold nodes and state transition matrices
    %%% corresponding to each
    [~, XSTM_nodes] = ode113(@Int_CR3BnSTM, linspace(0, Tp_PO_i_n, n_nodes+1), [X0_PO_i_n; stm0_colVec], options, prms);
    
    %%% Get the unstable eigenvector at the first node (initial condition)
    monodromy_N0 = reshape(XSTM_nodes(end,7:42),6,6);
    [eVec_stable_N0, eVec_unstable_N0] = getStableAndUnstableEigenvectors(monodromy_N0);
    
    %%% Use state transition matrix to get stable and unstable eigenvectors 
    %%% at each node
    unstableEigenvectors = NaN(6, n_nodes);
    for node_i = 1:n_nodes
        %%% Grab STM from T0 to current node
        stm_Ni_T0 = reshape(XSTM_nodes(node_i, 7:42),6,6);
        
        %%% Use STM to propagate the unstable eigenvector from N0 to the
        %%% current node
        unstableEigenvectors(:,node_i) = stm_Ni_T0 * eVec_unstable_N0;
        
        %%% Turn unstable eigenvector into a unit vector
        unstableEigenvectors(:,node_i) = unstableEigenvectors(:,node_i) ./ norm(unstableEigenvectors(:,node_i));
    end 
    
    %%% Get rid of the repeat final/first point
    X_nodes = XSTM_nodes(1:end-1,1:6);
    
    % --------------------------
    %%% Preallocating for parfor
    % --------------------------
    traj_POi = cell(n_nodes,1);
    
    % --------------------------
    %%% Looping through manifolds ICs in parallel, computing the manifold
    %%% direction, and propagating forward in time
    % --------------------------
    parfor node_i = 1:n_nodes
        % --------------------------
        %%% Preallocate - will only be populated if there's an impact
        % --------------------------
        X0_impact          = [];
        X_impact_n         = [];
        latLon_p_deg       = [];   latLon_m_deg       = [];
        impactAngle_p_deg  = [];   impactAngle_m_deg  = [];
        latLonHeadingHat_p = [];   latLonHeadingHat_m = [];
        impactColor_p      = [];   impactColor_m      = [];
        impactBin_p        = [];   impactBin_m        = [];
        conditionCounter   = 1;
        
        %%% Set IC of current manifold node
        X0_manNode_i_n = X_nodes(node_i,:)';
        
        %%% Create the perturbed initial condition and integrate the manifold
        X0_man_unstable_p_i = X0_manNode_i_n + unstableEigenvectors(:,node_i).*pertScale;
        X0_man_unstable_m_i = X0_manNode_i_n - unstableEigenvectors(:,node_i).*pertScale;
        
        %%% Integrate the manifolds
        [~, X_man_unstable_p, t_p_event, X_p_event, index_p_event] = ode113(@Int_CR3Bn, [0, Tf_manifolds_n], X0_man_unstable_p_i, options_impact, prms);
        [~, X_man_unstable_m, t_m_event, X_m_event, index_m_event] = ode113(@Int_CR3Bn, [0, Tf_manifolds_n], X0_man_unstable_m_i, options_impact, prms);
        
        %%% If desired, store trajectory
        if store_and_plot_manifold_trajectories
            traj_POi{node_i}.r_n = [X_man_unstable_p(:,1:3); NaN(1,3); X_man_unstable_m(:,1:3)];
        end
        
        % --------------------------
        %%% Get impact conditions and information for each manifold set
        % --------------------------
        doubleImpactCounter = 0;
        if ~isempty(X_p_event)
            doubleImpactCounter = doubleImpactCounter + 1;
            [latLon_p_deg, impactAngle_p_deg, latLonHeadingHat_p, impactColor_p] = getImpactConditions(X_p_event, prms, impactAngleBins_deg, impactAngleColors);
            
            X0_impact  = X0_man_unstable_p_i;
            if store_and_plot_impact_trajectories
                X_impact_n = X_man_unstable_p(:,1:3);
            end
            
            impactBin_p = discretize(impactAngle_p_deg, impactAngleBins_deg);
        end
        
        if ~isempty(X_m_event)
            doubleImpactCounter = doubleImpactCounter + 1;
            [latLon_m_deg, impactAngle_m_deg, latLonHeadingHat_m, impactColor_m] = getImpactConditions(X_m_event, prms, impactAngleBins_deg, impactAngleColors);
            
            X0_impact  = X0_man_unstable_m_i;
            if store_and_plot_impact_trajectories
                X_impact_n = X_man_unstable_m(:,1:3);
            end
            
            impactBin_m = discretize(impactAngle_m_deg, impactAngleBins_deg);
        end
        
        if doubleImpactCounter == 2
            conditionCounter = 2;
        end
        
        
        % --------------------------
        %%% Store Impact Data
        % --------------------------
        traj_POi{node_i}.latLon_deg       = [latLon_p_deg; latLon_m_deg];
        traj_POi{node_i}.impactAngle_deg  = [impactAngle_p_deg; impactAngle_m_deg];
        traj_POi{node_i}.latLonHeadingHat = [latLonHeadingHat_p; latLonHeadingHat_m];
        traj_POi{node_i}.impactColor      = [impactColor_p; impactColor_m];
        traj_POi{node_i}.conditionCounter = conditionCounter;
        traj_POi{node_i}.impactBin        = [impactBin_p; impactBin_m];
            
        if doubleImpactCounter ~=2
            traj_POi{node_i}.X0_n             = X0_impact;
            if store_and_plot_impact_trajectories
                traj_POi{node_i}.X_impact_n   = X_impact_n;
            end
        elseif doubleImpactCounter == 2
            traj_POi{node_i}.X0_n             = [X0_man_unstable_p_i, X0_man_unstable_m_i];
            if store_and_plot_impact_trajectories
                traj_POi{node_i}.X_impact_n   = {X_man_unstable_p(:,1:3), X_man_unstable_m(:,1:3)};
            end
        end
    end
    
    % -------------------------------------------------
    %%% Rearrange data from parfor structure
    % -------------------------------------------------
    %%% Restructuring trajectory data from parfor - now one large struct
    %%% array
    traj_POi_structArray = [traj_POi{:}];
    
    %%% Counting conditions ... this value would equal n_manifolds, but in
    %%% some specific instances, manifolds in both the 'positive' and
    %%% 'negative' unstable direction end up impacting, so you get two
    %%% impact results from one point along the periodic orbit. Therefore,
    %%% this value will be greater than or equal to n_manifolds.
    n_conditions = sum([traj_POi_structArray(:).conditionCounter]);
    
    %%% Grabbing specific data sets and turning them into single matrices
    X0s_POi_n             = NaN(n_conditions,6);
    latLons_POi_deg       = NaN(n_conditions,2);
    impactAngles_POi_deg  = NaN(n_conditions,1);
    latLonHeadingHats_POi = NaN(n_conditions,2);
    impactColors_POi      = NaN(n_conditions,3);
    impactBins_POi        = NaN(n_conditions,1);
    if store_and_plot_impact_trajectories
        Ximpacts_n_POi    = {};
        impactTrajCounter = 0;
    end
    
    condition_i = 0;
    for kk = 1:n_nodes
        condition_i = condition_i + 1;
        if ~isempty(traj_POi_structArray(kk).latLon_deg)
            X0s_POi_n(condition_i,:)             = traj_POi_structArray(kk).X0_n(:,1)';
            latLons_POi_deg(condition_i,:)       = traj_POi_structArray(kk).latLon_deg(1,:);
            impactAngles_POi_deg(condition_i,:)  = traj_POi_structArray(kk).impactAngle_deg(1);
            latLonHeadingHats_POi(condition_i,:) = traj_POi_structArray(kk).latLonHeadingHat(1,:);
            impactColors_POi(condition_i,:)      = traj_POi_structArray(kk).impactColor(1,:);
            impactBins_POi(condition_i,:)        = traj_POi_structArray(kk).impactBin(1);
            
            if store_and_plot_impact_trajectories
                impactTrajCounter = impactTrajCounter + 1;
                if size(traj_POi_structArray(kk).latLon_deg,1) == 1
                    Ximpacts_n_POi{impactTrajCounter} = traj_POi_structArray(kk).X_impact_n;
                elseif size(traj_POi_structArray(kk).latLon_deg,1) == 2
                    Ximpacts_n_POi{impactTrajCounter} = traj_POi_structArray(kk).X_impact_n{1};
                end
            end
            
            %%% For manifold ICs that impacted in both the 'plus' and
            %%% 'minus' directions. This section stores the results of the
            %%% 2nd case.
            if size(traj_POi_structArray(kk).latLon_deg,1) == 2
                condition_i = condition_i + 1;
                
                X0s_POi_n(condition_i,:)             = traj_POi_structArray(kk).X0_n(:,2)';
                latLons_POi_deg(condition_i,:)       = traj_POi_structArray(kk).latLon_deg(2,:);
                impactAngles_POi_deg(condition_i,:)  = traj_POi_structArray(kk).impactAngle_deg(2);
                latLonHeadingHats_POi(condition_i,:) = traj_POi_structArray(kk).latLonHeadingHat(2,:);
                impactColors_POi(condition_i,:)      = traj_POi_structArray(kk).impactColor(2,:);
                impactBins_POi(condition_i,:)        = traj_POi_structArray(kk).impactBin(2);

                if store_and_plot_impact_trajectories
                    impactTrajCounter = impactTrajCounter + 1;
                    Ximpacts_n_POi{impactTrajCounter} = traj_POi_structArray(kk).X_impact_n{2};
                end
            end
        end
    end
    
    % -------------------------------------------------
    %%% Store data from this periodic orbit in larger cell array that
    %%% contains data from all periodic orbits studied
    % -------------------------------------------------
    X0s_n{PO_i}             = X0s_POi_n;
    latLons_deg{PO_i}       = latLons_POi_deg;
    impactAngles_deg{PO_i}  = impactAngles_POi_deg;
    latLonHeadingHats{PO_i} = latLonHeadingHats_POi;
    impactColors{PO_i}      = impactColors_POi;
    impactBins{PO_i}        = impactBins_POi;
    if store_and_plot_impact_trajectories
        Ximpacts_n{PO_i}    = Ximpacts_n_POi;
    end
    
    if store_and_plot_manifold_trajectories
        trajectories_POi_n = cell(n_nodes,1);
        for kk = 1:n_nodes
            trajectories_POi_n{kk} = traj_POi_structArray(kk).r_n;
        end
        
        trajectories_n{PO_i} = trajectories_POi_n;
    end
end

% ========================================================================
%%% Plot impact map and trajectories
% ========================================================================
% --------------------------
%%% Creating a color matrix corresponding to ticks of impact angle bins
% --------------------------
colorMatrixForColorbar = NaN(impactAngleBins_deg(end),3);
for kk = 1:impactAngleBins_deg(end)
    for jj = 2:length(impactAngleBins_deg)
        if (kk - impactAngleBins_deg(jj)) <= 0
            colorMatrixForColorbar(kk,:) = impactAngleColors(jj-1,:);
            break
        end
    end
end

% --------------------------
%%% Binning lat/lon impact data
% --------------------------
PObins_latLons_deg   = cell(length(PO_indicies),1);
PObins_latLonHeading = cell(length(PO_indicies),1);

for PO_i = 1:length(PO_indicies)
    latLon_plotData        = cell(n_impactAngleBins,1);
    latLonHeading_plotData = cell(n_impactAngleBins,1);
    for kk = 1:size(latLons_deg{PO_i},1)
%        if impactBins{PO_i}(kk) 
        if isnan(impactBins{PO_i}(kk)) == 0
            latLon_plotData{impactBins{PO_i}(kk)} = [latLon_plotData{impactBins{PO_i}(kk)}; latLons_deg{PO_i}(kk,:)];
            
            latLonHeading_plotData{impactBins{PO_i}(kk)} = [latLonHeading_plotData{impactBins{PO_i}(kk)}; latLonHeadingHats{PO_i}(kk,:)];
        end
    end
    PObins_latLons_deg{PO_i}   = latLon_plotData;
    PObins_latLonHeading{PO_i} = latLonHeading_plotData;
end

%%
%%% Make plots!
for PO_i = 1:length(PO_indicies)
    p = figure(100+PO_i); hold all
    p.Position = [772 385 560 420];
%     for kk = 1:size(latLons_deg{PO_i},1)
%         if ~isnan(latLons_deg{PO_i}(kk,2))
% 
%             plot(latLons_deg{PO_i}(kk,2), latLons_deg{PO_i}(kk,1),...
%                 'o','markersize',6,'markeredgecolor',colors.std.black,'markerfacecolor',impactColors{PO_i}(kk,:))
%             
%             quiver(latLons_deg{PO_i}(kk,2), latLons_deg{PO_i}(kk,1), ...
%                 latLonHeadingHats{PO_i}(kk,2), latLonHeadingHats{PO_i}(kk,1),...
%                 10, 'color', impactColors{PO_i}(kk,:),'linewidth',3)
%         end
%     end
    for bin_i = n_impactAngleBins:-1:1
        if ~isempty(PObins_latLons_deg{PO_i}{bin_i})
            plot(PObins_latLons_deg{PO_i}{bin_i}(:,2), PObins_latLons_deg{PO_i}{bin_i}(:,1),...
                '.','markersize',17,'markeredgecolor',colors.std.black,'markerfacecolor',impactAngleColors(bin_i,:))

            for quiver_i = 1:size(PObins_latLons_deg{PO_i}{bin_i},1)
                quiver(PObins_latLons_deg{PO_i}{bin_i}(quiver_i,2), PObins_latLons_deg{PO_i}{bin_i}(quiver_i,1), ...
                    PObins_latLonHeading{PO_i}{bin_i}(quiver_i,2), PObins_latLonHeading{PO_i}{bin_i}(quiver_i,1),...
                    10, 'color', impactAngleColors(bin_i,:),'linewidth',3)
            end
        end
    end
    
    xlim([-1, 1] .* 185)
    ylim([-1, 1] .* 95)
    PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 18, 'LaTex')
    cbar1 = colorbar;
    cbar1.FontName     = 'Arial';
    cbar1.FontSize     = 10;
    cbar1.Ticks        = impactAngleBins_deg./90;
    cbar1.TickLabels   = num2cell(impactAngleBins_deg);
    cbar1.Label.String = {'Impact Angle, deg'};
    cbar1.Label.Rotation = 0; % 0 = horizontal, 90 = vertical
    cbar1.Label.Position = [0.7, 1.05, 0];
    colormap(colorMatrixForColorbar)
    
    
    if store_and_plot_impact_trajectories
        p = figure(200+PO_i); hold all
        p.Position = [129 34 560 420];
        axis equal
        plotSecondary(secondary)
        PlotBoi3_CR3Bn(20)
        for man_i = 1:size(Ximpacts_n{PO_i},2)
            plot3(Ximpacts_n{PO_i}{man_i}(:,1), Ximpacts_n{PO_i}{man_i}(:,2), Ximpacts_n{PO_i}{man_i}(:,3), 'r', 'linewidth',0.5)
        end
        
        plot3(X_POs{PO_i}(:,1), X_POs{PO_i}(:,2), X_POs{PO_i}(:,3), 'b', 'linewidth',1)
    end
        
        
    if store_and_plot_manifold_trajectories
        p = figure(300+PO_i); hold all
        p.Position = [788 44 560 420];
        axis equal
        plotSecondary(secondary)
        PlotBoi3_CR3Bn(20)
        view(0,90)
        
        for man_i = 1:n_nodes
            plot3(trajectories_n{PO_i}{man_i}(:,1), trajectories_n{PO_i}{man_i}(:,2), trajectories_n{PO_i}{man_i}(:,3), 'r', 'linewidth',0.5)
        end
        
        plot3(X_POs{PO_i}(:,1), X_POs{PO_i}(:,2), X_POs{PO_i}(:,3), 'b', 'linewidth',1)
    end
end


if n_conditions > n_nodes
    warning('Some manifolds had double impacts')
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











% For Jay

%%% Europa - vertical manifold landing at north pole
% 500 manifolds
% Index 101, orbit 2 ... Ximpacts_n{1}{2}(:,1:3)




