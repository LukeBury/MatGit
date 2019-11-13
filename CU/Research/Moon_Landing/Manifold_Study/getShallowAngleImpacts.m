% ========================================================================
%%% Description
% ========================================================================
% This script is meant to propagate the unstable manifolds from families 
% of periodic orbits and track those that impact the surface of the
% secondary body. Of those that impact, those that arrive at an angle less
% than a set minimum value (originally 3 degrees) will have their
% information stored in a csv. This script is meant to run this process in
% parallel and write all outputs to files so that it may be run on remote
% machines.

% Created: 11/12/19
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

%%% Set paths based on computer
if isequal(computer,'MACI64')      % Mac
    mbinPath = '~/CU_Google_Drive/Documents/MatGit/mbin';
    POfamilyPath = '~/CU_Google_Drive/Documents/MatGit/mbin/Data/InitialConditions/PO_Families/';
    savepath = '~/CU_Google_Drive/Documents/MatGit/MatlabOutputs';
    computerTag = 'M';
    
elseif isequal(computer,'GLNXA64') % Fortuna
    mbinPath = '/home/lubu8198/MatGit/mbin';
    POfamilyPath = '/home/lubu8198/MatGit/mbin/Data/InitialConditions/PO_Families/';
    savepath = '/orc_raid/lubu8198/MatlabOutputs';
    computerTag = 'F';
%   N = maxNumCompThreads returns the current maximum number of computational threads N.
%   LASTN = maxNumCompThreads(N) sets the maximum number of computational threads to N, 
%        and returns the previous maximum number of computational threads, LASTN.
%   LASTN = maxNumCompThreads('automatic') sets the maximum number of computational threads
%        using what the MATLAB® software determines to be the most desirable. It additionally
%        returns the previous maximum number of computational threads, LASTN.

else 
    warning('This computer will explode in 5 seconds')
end

addpath(genpath(mbinPath))

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
save_data = 1;

% ========================================================================
%%% Setup
% ========================================================================
% -------------------------------------------------
%%% Options
% -------------------------------------------------
%%% Number of manifolds per PO
n_nodes = 2000;

%%% Scale the perturbation of the manifold node in the unstable direction
pertScale = 1e-8;

%%% Set propagation time for unstable manifolds
Tf_manifolds_n = 4*pi;

%%% Bins and colors for impact angles
impactAngleBins_deg = [0, 3, 10, 20, 45, 90];

n_impactAngleBins = length(impactAngleBins_deg) - 1;
impactAngleColors = colorScale([colors.std.cyan; colors.std.mag], n_impactAngleBins);

%%% 3B System
% famName_bodies = 'Earth_Moon';
% famName_bodies = 'Jupiter_Europa';
% famName_bodies = 'Jupiter_Ganymede';
famName_bodies = 'Saturn_Enceladus';
% famName_bodies = 'Neptune_Triton';

%%% PO Family
% famName_PO_Family = 'L1_Lyapunov';
% famName_PO_Family = 'L1_Vertical';
% famName_PO_Family = 'L1_SHalo';
% famName_PO_Family = 'L2_Lyapunov';
famName_PO_Family = 'L2_Vertical';
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
    PO_indicies = [340];
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
    PO_indicies = linspace(1,size(POFamilyData,1),size(POFamilyData,1));
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
Xfs_n             = cell(length(PO_indicies),1);
latLons_deg       = cell(length(PO_indicies),1);
impactAngles_deg  = cell(length(PO_indicies),1);
latLonHeadingHats = cell(length(PO_indicies),1);
impactBins        = cell(length(PO_indicies),1);

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
    if isempty(eVec_stable_N0)
        continue
    end
    
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
        Xf_impact          = [];
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
        
        % --------------------------
        %%% Get impact conditions and information for each manifold set
        % --------------------------
        doubleImpactCounter = 0;
        if ~isempty(X_p_event)
            doubleImpactCounter = doubleImpactCounter + 1;
            [latLon_p_deg, impactAngle_p_deg, latLonHeadingHat_p, impactColor_p] = getImpactConditions(X_p_event, prms, impactAngleBins_deg, impactAngleColors);
            
            X0_impact  = X0_man_unstable_p_i;
            Xf_impact  = X_p_event(end,:)';
            
            impactBin_p = discretize(impactAngle_p_deg, impactAngleBins_deg);
        end
        
        if ~isempty(X_m_event)
            doubleImpactCounter = doubleImpactCounter + 1;
            [latLon_m_deg, impactAngle_m_deg, latLonHeadingHat_m, impactColor_m] = getImpactConditions(X_m_event, prms, impactAngleBins_deg, impactAngleColors);
            
            X0_impact  = X0_man_unstable_m_i;
            Xf_impact  = X_m_event(end,:)';
          
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
            traj_POi{node_i}.Xf_n             = Xf_impact;

        elseif doubleImpactCounter == 2
            traj_POi{node_i}.X0_n             = [X0_man_unstable_p_i, X0_man_unstable_m_i];
            traj_POi{node_i}.Xf_n             = [X_p_event(end,:)', X_m_event(end,:)'];
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
    Xfs_POi_n             = NaN(n_conditions,6);
    latLons_POi_deg       = NaN(n_conditions,2);
    impactAngles_POi_deg  = NaN(n_conditions,1);
    latLonHeadingHats_POi = NaN(n_conditions,2);
    impactBins_POi        = NaN(n_conditions,1);
    
    condition_i = 0;
    for kk = 1:n_nodes
        condition_i = condition_i + 1;
        if ~isempty(traj_POi_structArray(kk).latLon_deg)
            X0s_POi_n(condition_i,:)             = traj_POi_structArray(kk).X0_n(:,1)';
            Xfs_POi_n(condition_i,:)             = traj_POi_structArray(kk).Xf_n(:,1)';
            latLons_POi_deg(condition_i,:)       = traj_POi_structArray(kk).latLon_deg(1,:);
            impactAngles_POi_deg(condition_i,:)  = traj_POi_structArray(kk).impactAngle_deg(1);
            latLonHeadingHats_POi(condition_i,:) = traj_POi_structArray(kk).latLonHeadingHat(1,:);
            impactBins_POi(condition_i,:)        = traj_POi_structArray(kk).impactBin(1);
            
            %%% For manifold ICs that impacted in both the 'plus' and
            %%% 'minus' directions. This section stores the results of the
            %%% 2nd case.
            if size(traj_POi_structArray(kk).latLon_deg,1) == 2
                condition_i = condition_i + 1;
                
                X0s_POi_n(condition_i,:)             = traj_POi_structArray(kk).X0_n(:,2)';
                Xfs_POi_n(condition_i,:)             = traj_POi_structArray(kk).Xf_n(:,2)';
                latLons_POi_deg(condition_i,:)       = traj_POi_structArray(kk).latLon_deg(2,:);
                impactAngles_POi_deg(condition_i,:)  = traj_POi_structArray(kk).impactAngle_deg(2);
                latLonHeadingHats_POi(condition_i,:) = traj_POi_structArray(kk).latLonHeadingHat(2,:);
                impactBins_POi(condition_i,:)        = traj_POi_structArray(kk).impactBin(2);

            end
        end
    end
    
    % -------------------------------------------------
    %%% Store data from this periodic orbit in larger cell array that
    %%% contains data from all periodic orbits studied
    % -------------------------------------------------
    X0s_n{PO_i}             = X0s_POi_n;
    Xfs_n{PO_i}             = Xfs_POi_n;
    latLons_deg{PO_i}       = latLons_POi_deg;
    impactAngles_deg{PO_i}  = impactAngles_POi_deg;
    latLonHeadingHats{PO_i} = latLonHeadingHats_POi;
    impactBins{PO_i}        = impactBins_POi;
end


% --------------------------
%%% Binning lat/lon impact data
% --------------------------
PObins_latLons_deg   = cell(length(PO_indicies),1);
PObins_latLonHeading = cell(length(PO_indicies),1);

allImpactData = [];


for PO_i = 1:length(PO_indicies)
%     latLon_plotData        = cell(n_impactAngleBins,1);
%     latLonHeading_plotData = cell(n_impactAngleBins,1);
    for kk = 1:size(latLons_deg{PO_i},1)
% %        if impactBins{PO_i}(kk) 
%         if ~isnan(impactBins{PO_i}(kk))
%             latLon_plotData{impactBins{PO_i}(kk)} = [latLon_plotData{impactBins{PO_i}(kk)}; latLons_deg{PO_i}(kk,:)];
%             
%             latLonHeading_plotData{impactBins{PO_i}(kk)} = [latLonHeading_plotData{impactBins{PO_i}(kk)}; latLonHeadingHats{PO_i}(kk,:)];
%         end
        
        if impactBins{PO_i}(kk) == 1
            JC_current     = getJacobiConstant_ZH(X0s_n{PO_i}(kk,:), prms);
            landingVel_mps = norm(Xfs_n{PO_i}(kk,4:6)) * vNorm * 1000;
            
            % x0, y0, z0, xd0, yd0, zd0, lat, lon, impactAngle, heading_x, heading_y, JC, landingVel_mps
            allImpactData  = [allImpactData;...
                            X0s_n{PO_i}(kk,:), latLons_deg{PO_i}(kk,:), impactAngles_deg{PO_i}(kk,:), latLonHeadingHats{PO_i}(kk,:), JC_current, landingVel_mps];
                        
        end
        
    end
%     PObins_latLons_deg{PO_i}   = latLon_plotData;
%     PObins_latLonHeading{PO_i} = latLonHeading_plotData;
end

% -------------------------------------------------
% Writing data to CSV
% -------------------------------------------------
if save_data
    %%% Create File Name
    filename_shallowImpact = fullfile(savepath, sprintf('shallowImpacts.%s.%s.nodes%1.0f.txt', computerTag, famName, n_nodes));

    %%% Opening files and writing header
    f_shallowImpacts = fopen(filename_shallowImpact, 'wt');
    fprintf(f_shallowImpacts,'x0,y0,z0,xd0,yd0,zd0,lat,lon,impactAngle_deg,heading_x,heading_y,JC,landingVel_mps,...runTime=%1.0sec\n',toc(ticWhole));  % header

    for kk = 1:size(allImpactData,1)
        fprintf(f_shallowImpacts, '%1.16f,%1.16f,%1.16f,%1.16f,%1.16f,%1.16f,%1.2f,%1.2f,%1.2f,%1.6f,%1.6f,%1.16f,%1.6f\n',...
            allImpactData(kk,1),allImpactData(kk,2),allImpactData(kk,3),allImpactData(kk,4),allImpactData(kk,5),allImpactData(kk,6),allImpactData(kk,7),...
            allImpactData(kk,8),allImpactData(kk,9),allImpactData(kk,10),allImpactData(kk,11),allImpactData(kk,12),allImpactData(kk,13));
    end

    %%% Close file
    fclose(f_shallowImpacts);

end % save_data
% ========================================================================
%%% Closeout
% ========================================================================
tocWhole = toc(ticWhole);

if isequal(computer,'MACI64')
    fprintf('Elapsed time: %1.4f seconds\n',tocWhole)
end














