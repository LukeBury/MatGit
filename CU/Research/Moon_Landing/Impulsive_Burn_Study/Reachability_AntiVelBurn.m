% ========================================================================
%%% Description
% ========================================================================
% In this script, a reference trajectory is divided descretized into some
% number of points. From each of these points, a DV of some specified
% magnitude is applied in the anti-velocity direction, and each resulting
% trajectory is propagated forward in time. Resulting impact conditions are
% analyzed.

% ========================================================================
%%% Initialization
% ========================================================================
clear
clc
close all
ticWhole = tic;

% ========================================================================
%%% Run Switches
% ========================================================================
testCaseOn = 0;

% ========================================================================
%%% Setup and Options
% ========================================================================
% -------------------------------------------------
%%% Setting Paths and Importing Data
% -------------------------------------------------
%%% Set paths based on computer
if isequal(computer,'MACI64')      % Mac
    on_Fortuna         = 0;
    mbinPath = '~/CU_Google_Drive/Documents/MatGit/mbin';
    savepath = '~/CU_Google_Drive/Documents/MatGit/MatlabOutputs';
    computerTag = 'M';
elseif isequal(computer,'GLNXA64') % Fortuna
    on_Fortuna         = 1;
    mbinPath = '/home/lubu8198/MatGit/mbin';
    savepath = '/orc_raid/lubu8198/MatlabOutputs';
    computerTag = 'F';
else 
    warning('This computer will explode in 5 seconds')
    return
end

%%% Add the function paths to matlab
addpath(genpath(mbinPath))

%%% General data on solar system bodies
bodies = getBodyData(mbinPath);

%%% Color options/schemes
colors = get_colors();

% -------------------------------------------------
%%% Initial Conditions
% -------------------------------------------------
X0_n = [1.020461701526617; -0.001741448390222; 0; -0.000666124689955; -0.002025761973389; 0];

% -------------------------------------------------
%%% System
% -------------------------------------------------
%%% 3B system
primary   = bodies.jupiter;
secondary = bodies.europa;

%%% Normalizing constants
rNorm = secondary.a;         % n <-> km
tNorm = 1/secondary.meanMot; % n <-> sec
vNorm = rNorm / tNorm;       % n <-> km/sec

%%% Setting new variables because parfor is picky
R1_n = primary.R/rNorm;
R2_n = secondary.R_n;
MR = secondary.MR;
secondaryName = secondary.name;

%%% Lagrange point of interest
LPoint = 2;

%%% Acquire Collinear Lagrange points
L123 = EquilibriumPoints(secondary.MR,1:3); % [3x3] of L1, L2, L3 normalized BCR coordinates

% -------------------------------------------------
%%% Search settings
% -------------------------------------------------
%%% Delta-V magnitudes
DV_vec_mps = [200]; % m/s

%%% Resolution settings
if testCaseOn == 0
    n_DV_locations     = 1000;
    n_DVs_per_location = 1;
elseif testCaseOn == 1 %%% Lowering resolution for test cases
    n_DV_locations     = 10;
    n_DVs_per_location = 1;
else
    warning('bad testCaseOn setting')
    
end

% -------------------------------------------------
%%% 3-Body-Based System Settings
% -------------------------------------------------
if testCaseOn == 0
    if isequal(secondary.name,'europa')
    %     d
    elseif isequal(secondary.name,'enceladus')
    %     d
    elseif isequal(secondary.name,'titan')
    %     d
    elseif isequal(secondary.name,'triton')
    %     d
    end

elseif testCaseOn == 1 %%% Lowering resolution for test cases
    if isequal(secondary.name,'europa')
%         d
    elseif isequal(secondary.name,'enceladus')
%         d
    elseif isequal(secondary.name,'titan')
%         d
    elseif isequal(secondary.name,'triton')
%         d
    end
end

% -------------------------------------------------
%%% Integration settings
% -------------------------------------------------
%%% Time vector
t0 = 0;
tf_4pi = 4*pi;
tf_6pi = 6*pi;
n_t_ref = 20000;
time0_ref_4pi_n = linspace(t0, tf_4pi, n_t_ref);

%%% Choosing ode tolerance
tol = 1e-13;

%%% Setting integrator options
options_ImpactEscape = odeset('Events',@event_ImpactEscape_CR3Bn,'RelTol',tol,'AbsTol',tol);

% ========================================================================
%%% Calculations prior to looping through positions and burns
% ========================================================================
% -------------------------------------------------
%%% Energy info of X0_n
% -------------------------------------------------
%%% Initial Jacobi Constant
JC_scInitial = JacobiConstantCalculator(secondary.MR, X0_n(1:3)', X0_n(4:6)');

%%% L2 Flyover velocity
[L2FlyoverVelocity_mps] = JC_2_L2FlyoverVelocity(JC_scInitial, secondary.MR, L123(LPoint,:), vNorm);

% -------------------------------------------------
%%% Integrating Reference Trajectory
% -------------------------------------------------
%%% Setting necessary parameters for integration
prms.u    = secondary.MR;
prms.R2_n = secondary.R_n;
prms.L1x  = L123(1,1);
prms.L2x  = L123(2,1);

%%% Integrate reference trajectory until impact
[time_ref_n, X_ref_BCR_n] = ode113(@Int_CR3Bn, time0_ref_4pi_n, X0_n, options_ImpactEscape, prms);

% -------------------------------------------------
%%% Discretize reference trajectory by path distance
% % % % -------------------------------------------------
% % % %%% Use time of impact to create discretized time vector for reference traj
% % % time0_discretize_n = linspace(t0, time_ref_n(end), n_DV_locations);
% % % 
% % % %%% Discretize reference trajectory
% % % [time_discretize_n, X_discretize_BCR_n] = ode113(@Int_CR3Bn, time0_discretize_n, X0_n, options_ImpactEscape, prms);

%%% Discretize reference trajectory by evenly spaced points in path
%%% distance
[X_discretize_BCR_n] = discretize_by_path_distance(X_ref_BCR_n, n_DV_locations+2);

%%% Cutting out first and final state, since we don't care about a burn
%%% from X0_ref or from the surface of the secondary
X_discretize_BCR_n = X_discretize_BCR_n(2:end-1,:);

%%% Find corresponding time vector
time_discretize_n = zeros(n_DV_locations,1);
for kk = 1:n_DV_locations
    closestIndex = rowNorm(X_ref_BCR_n(:,1:3) - X_discretize_BCR_n(kk,1:3)) == min(rowNorm(X_ref_BCR_n(:,1:3) - X_discretize_BCR_n(kk,1:3)));
    time_discretize_n(kk) = time_ref_n(closestIndex);
end

% -------------------------------------------------
%%% Computing reference parameters
% -------------------------------------------------
% --------------------------
%%% Impact Longitude - reference
% --------------------------
%%% Finding lat/lon of impact site
[impactLat_ref_deg, impactLon_ref_deg] = BCR2latlon(X_ref_BCR_n(end,1:3), 'secondary', MR);         

% --------------------------
%%% Impact Angle - reference
% --------------------------
%%% Creating SCR position vector
rImpact_ref_SCR_n = X_ref_BCR_n(end,1:3) - [1-MR,0,0];
rHatImpact_ref_SCR_n = rImpact_ref_SCR_n./norm(rImpact_ref_SCR_n);

%%% Velocity unit vector at impact
vHatImpact_ref_n = X_ref_BCR_n(end,4:6)./norm(X_ref_BCR_n(end,4:6));

%%% Angle between velocity and surface
[impactAngle_ref_deg] = calcImpactAngle(rHatImpact_ref_SCR_n,vHatImpact_ref_n,'degrees');

% -------------------------------------------------
%%% Quick calculations before parallel loop
% -------------------------------------------------
%%% Converting DV to normalized units
DV_vec_n = (DV_vec_mps/1000)/vNorm;

% ========================================================================
%%% Parallel loop through positions and burns
% ========================================================================
%%% Create data structure
positionData = {};

%%% Loop through positions in parallel
parfor position_ii = 1:n_DV_locations
    % -------------------------------------------------
    % Reducing broadcast variables
    % -------------------------------------------------
    DV_vec_n_parfor = DV_vec_n;
    
    % -------------------------------------------------
    % Grabbing current initial position and time from discretized reference
    % trajectory
    % -------------------------------------------------
    X0_ref_ii = X_discretize_BCR_n(position_ii,:);
    t0_ref   = time_discretize_n(position_ii);
    
    %%% Create new time vector for integration
    time0n_ii = linspace(t0_ref,tf_6pi,1000);
    
    % -------------------------------------------------
    % Initialization data structures for current position
    % -------------------------------------------------
    XBurn_n_ii         = zeros(n_DVs_per_location,6);
    R_Burn_SCR_ii      = zeros(n_DVs_per_location,1);
    T_Burn_n_ii        = zeros(n_DVs_per_location,1);
    dJC_ii             = zeros(n_DVs_per_location,1);
    impactLon_deg_ii   = zeros(n_DVs_per_location,1);
    impactAngle_deg_ii = zeros(n_DVs_per_location,1);
    TOF_n_ii           = zeros(n_DVs_per_location,1);
    
    % -------------------------------------------------
    % Applying each burn magnitude to the current reference state and
    % propogating it forward
    % -------------------------------------------------
    for burnMag_kk = 1:n_DVs_per_location
        %%% Preallocating temporary variables
        R_Burn_SCR      = [];
        T_Burn_n        = [];
        dJC             = [];
        impactLon_deg   = [];
        impactAngle_deg = [];
        TOF_n           = [];
            
        %%% Unit vector of reference velocity
        v0_ref_hat = X0_ref_ii(4:6)./norm(X0_ref_ii(4:6));
        
        %%% Burn vector
        burnMagnitude = DV_vec_n_parfor(burnMag_kk);
        burn = -v0_ref_hat .* burnMagnitude;
        
        %%% New state
        XBurn_BCR_ii_kk = [X0_ref_ii(1:3), X0_ref_ii(4:6) + burn]';
        
        %%% Post-burn energy
        JC_ii_kk_sc = JacobiConstantCalculator(MR, XBurn_BCR_ii_kk(1:3)', XBurn_BCR_ii_kk(4:6)');
        
        %%% Integrate new state
        [time_ii_kk_n, X_ii_kk_BCR_n, time_eventImpact, X_eventImpact, index_eventImpact] = ode113(@Int_CR3Bn, time0n_ii, XBurn_BCR_ii_kk, options_ImpactEscape, prms);
        
        % -------------------------------------------------
        % Determining end conditions
        % -------------------------------------------------
        impact1_escape2_orbit3 = []; % 1-impact, 2-escape, 3-still in orbit
        if isempty(X_eventImpact) == 0
            if abs(norm(X_eventImpact(end,1:3)-[1-MR,0,0])-R2_n) < 1e-9 % impact
                impact1_escape2_orbit3 = 1;
            elseif X_eventImpact(end,1) <= prms.L1x || X_eventImpact(end,1) >= prms.L2x % escape
                impact1_escape2_orbit3 = 2;
            end
        elseif isempty(X_eventImpact) == 1 % No event yet, still in orbit
            impact1_escape2_orbit3 = 3;
        end
        
        % -------------------------------------------------
        % Studying end conditions
        % -------------------------------------------------
        if impact1_escape2_orbit3 == 1
            % --------------------------
            % Burn parameters
            % --------------------------
            %%% Distance from secondary at time of burn
            R_Burn_SCR = norm(XBurn_BCR_ii_kk(1:3) - [1-MR; 0; 0]);
            
            %%% Time of burn since beginning of reference trajectory
            T_Burn_n = time0n_ii(1);
            
            %%% Difference in Jacobi constant
            JC_postBurn = JacobiConstantCalculator(MR, XBurn_BCR_ii_kk(1:3)', XBurn_BCR_ii_kk(4:6)');
            dJC = JC_ii_kk_sc - JC_scInitial;
            
            % --------------------------
            % Impact Longitude
            % --------------------------
            %%% Finding lat/lon of impact site
            [impactLat_deg, impactLon_deg] = BCR2latlon(X_eventImpact(end,1:3), 'secondary', MR);     
                        
            % --------------------------
            % Impact Angle
            % --------------------------
            %%% Creating SCR position vector
            rImpact_SCR_n = X_eventImpact(end,1:3) - [1-MR,0,0];
            rHatImpact_SCR_n = rImpact_SCR_n./norm(rImpact_SCR_n);

            %%% Velocity unit vector at impact
            vHatImpact_n = X_eventImpact(end,4:6)./norm(X_eventImpact(end,4:6));

            %%% Angle between velocity and surface
            [impactAngle_deg] = calcImpactAngle(rHatImpact_SCR_n,vHatImpact_n,'degrees');
            
            % --------------------------
            % TOF
            % --------------------------
            %%% Storing endTime
            TOF_n = time_ii_kk_n(end);
               
        elseif impact1_escape2_orbit3 == 2 % escape
            R_Burn_SCR      = NaN;
            T_Burn_n        = NaN;
            dJC             = NaN;
            impactLon_deg   = NaN;
            impactAngle_deg = NaN;
            TOF_n           = NaN;

        elseif impact1_escape2_orbit3 == 3 % orbit
            R_Burn_SCR      = NaN;
            T_Burn_n        = NaN;
            dJC             = NaN;
            impactLon_deg   = NaN;
            impactAngle_deg = NaN;
            TOF_n           = NaN;
            
        end
        
        % -------------------------------------------------
        % Storing results of trajectory
        % -------------------------------------------------
        XBurn_n_ii(burnMag_kk,:)           = XBurn_BCR_ii_kk';
        R_Burn_SCR_ii(burnMag_kk)          = R_Burn_SCR;
        T_Burn_n_ii(burnMag_kk)            = T_Burn_n;
        dJC_ii(burnMag_kk)                 = dJC;
        impactLon_deg_ii(burnMag_kk)       = impactLon_deg;
        impactAngle_deg_ii(burnMag_kk)     = impactAngle_deg;
        TOF_n_ii(burnMag_kk)               = TOF_n;

    end % burnMag_kk = 1:n_DVs_per_location
    
    
    % -------------------------------------------------
    % Storing data from all azimuths/elevations
    % -------------------------------------------------
    positionData{position_ii}.XBurn_n         = XBurn_n_ii;
    positionData{position_ii}.R_Burn_SCR      = R_Burn_SCR_ii;
    positionData{position_ii}.T_Burn_n        = T_Burn_n_ii;
    positionData{position_ii}.dJC             = dJC_ii;
    positionData{position_ii}.impactLon_deg   = impactLon_deg_ii;
    positionData{position_ii}.impactAngle_deg = impactAngle_deg_ii;
    positionData{position_ii}.TOF_n           = TOF_n_ii;
    
end % parfor position_ii = 2:n_DV_locations

% ========================================================================
%%% Post-Processing
% ========================================================================
% -------------------------------------------------
% Combine Data
% -------------------------------------------------
%%% Turning data in one large struct array
positionData_structArray = [positionData{:}];

%%% Seperating fields into vectors & matrices
XBurn_n_mat          = [positionData_structArray(:).XBurn_n];
XBurn_n_mat          = reshape(XBurn_n_mat,6,1000);
XBurn_n_mat          = XBurn_n_mat';

R_Burn_SCR_vec       = [positionData_structArray(:).R_Burn_SCR];
R_Burn_SCR_radii_vec = R_Burn_SCR_vec ./ secondary.R_n;
T_Burn_n_vec         = [positionData_structArray(:).T_Burn_n];
dJC_vec              = [positionData_structArray(:).dJC];
impactLon_deg_vec    = [positionData_structArray(:).impactLon_deg];
impactAngle_deg_vec  = [positionData_structArray(:).impactAngle_deg];
TOF_n_vec            = [positionData_structArray(:).TOF_n];

%%% New data vectors - post burn time of flight
TOF_postBurn_n_vec = TOF_n_vec - T_Burn_n_vec;

%%% New data vectors - Approx Burn efficiency (landing DV saved / burn
%%% DV)... for this analysis, just using position of reference landing
%%% rather than look at every single new landing spot. Technically this
%%% isn't a "fair" comparision fundamentally since the burn changes
%%% where you land
JC_postBurn_vec          = dJC_vec + JC_scInitial;
JC_referenceLandingSite  = JacobiConstantCalculator(secondary.MR, X_ref_BCR_n(end,1:3), [0, 0, 0]);
referenceLandingDV_n     = sqrt(JC_referenceLandingSite - JC_scInitial);
postBurnLandingDV_n_vec  = sqrt(JC_referenceLandingSite - JC_postBurn_vec);
burnEfficiencies_vec     = (referenceLandingDV_n - postBurnLandingDV_n_vec) ./ DV_vec_n(1);

% -------------------------------------------------
% Assign each data point a color based on dJC
% -------------------------------------------------
%%% Finding max and min bounds of post-burn JC change
max_dJC = max(dJC_vec);
min_dJC = min(dJC_vec);

%%% Assigning percentiles of post-burn JC change
dJC_percentiles = (dJC_vec - min_dJC) / (max_dJC - min_dJC);

%%% Making it a column vector
dJC_percentiles = dJC_percentiles(:);

%%% Choosing colors
col1 = colors.std.mag;
col2 = colors.std.cyan;
dcol = col2 - col1;

%%% Assigning color to each data point
dJC_colors_mat = dJC_percentiles.*dcol + col1;

% ========================================================================
%%% Plotting ... (gaps in data represent non-impacts)
% ========================================================================
scatSize = 100;
% -------------------------------------------------
% Reference trajectory
% -------------------------------------------------
figure; hold all
p1 = plot(X_ref_BCR_n(:,1), X_ref_BCR_n(:,2),'r','linewidth',2);
plotCR3BP_Neck(secondary,L123,JC_scInitial,600,200,colors.std.black,1.5)
plotBody2(secondary.R_n, [1-secondary.MR,0,0],colors.std.white,colors.std.blue,2,0)
PlotBoi2('$X_n$','$Y_n$',20,'LaTex')
axis equal
legend([p1],'Reference')

% -------------------------------------------------
% Impact Longitude vs Impact Angle
% -------------------------------------------------
figure; hold all
scatter(fliplr(impactLon_deg_vec), fliplr(impactAngle_deg_vec),scatSize, flipud(dJC_colors_mat),'.');
getColorBarNOW(col1, col2, burnEfficiencies_vec)
p2 = plot(impactLon_ref_deg, impactAngle_ref_deg,'rx','markersize',20,'linewidth',3);
xlim([-180 180])
ylim([0 90])
PlotBoi2('Impact Longitude, $^\circ$','Impact Angle, $^\circ$',18,'LaTex')
legend([p2],'Reference')

% -------------------------------------------------
% Burn Time vs Burn Distance from Secondary
% -------------------------------------------------
figure; hold all
scatter(T_Burn_n_vec./pi, R_Burn_SCR_radii_vec,scatSize, dJC_colors_mat,'.')
getColorBarNOW(col1, col2, burnEfficiencies_vec)
ylim([0 max(R_Burn_SCR_radii_vec)*1.1])
PlotBoi2('Time of Burn, $\pi$','Distance from Secondary at Burn, radii',18,'LaTex')

% -------------------------------------------------
% Change in JC vs Burn Distance from Secondary
% -------------------------------------------------
figure; hold all
scatter(dJC_vec, R_Burn_SCR_radii_vec,scatSize, dJC_colors_mat,'.')
getColorBarNOW(col1, col2, burnEfficiencies_vec)
ylim([0 max(R_Burn_SCR_radii_vec)*1.1])
PlotBoi2('$\Delta$JC','Distance from Secondary at Burn, radii',18,'LaTex')

% -------------------------------------------------
% Impact Longitude vs TOF
% -------------------------------------------------
figure; hold all
scatter(impactLon_deg_vec, TOF_n_vec./pi,scatSize, dJC_colors_mat,'.');
getColorBarNOW(col1, col2, burnEfficiencies_vec)
p2 = plot(impactLon_ref_deg, time_ref_n(end)./pi,'rx','markersize',20,'linewidth',3);
xlim([-180 180])
PlotBoi2('Impact Longitude, $^\circ$','Total TOF, $\pi$',18,'LaTex')
legend([p2],'Reference')

% -------------------------------------------------
% TOF vs Impact Angle
% -------------------------------------------------
figure; hold all
scatter(TOF_n_vec./pi, impactAngle_deg_vec,scatSize, dJC_colors_mat,'.');
getColorBarNOW(col1, col2, burnEfficiencies_vec)
p2 = plot(time_ref_n(end)./pi, impactAngle_ref_deg,'rx','markersize',20,'linewidth',3);
ylim([0 90])
PlotBoi2('Total TOF, $\pi$', 'Impact Angle, $^\circ$', 18, 'LaTex')
legend([p2],'Reference')

% -------------------------------------------------
% Burn Efficiency vs Burn distance from Secondary
% -------------------------------------------------
figure; hold all
scatter(burnEfficiencies_vec, R_Burn_SCR_radii_vec,scatSize,dJC_colors_mat,'.');
getColorBarNOW(col1, col2, burnEfficiencies_vec)
ylim([0 max(R_Burn_SCR_radii_vec)*1.1])
PlotBoi2('Landing $\Delta V$ Saved / Burn $\Delta V$','Distance from Secondary at Burn, radii',18,'LaTex')




% ========================================================================
%%% Various Tests
% ========================================================================
longerTrajIndices   = TOF_n_vec > time_ref_n(end);
XBurns_n_longerTraj = XBurn_n_mat(longerTrajIndices,:);
TBurns_n_longerTraj = T_Burn_n_vec(longerTrajIndices);
TOFs_n_longerTraj   = TOF_n_vec(longerTrajIndices);

longerTrajectories = [];
for kk = 1:length(TBurns_n_longerTraj)
    [time_n_longerTraj, X_BCR_n_longerTraj] = ode113(@Int_CR3Bn, linspace(TBurns_n_longerTraj(kk), TOFs_n_longerTraj(kk),500), XBurns_n_longerTraj(kk,:)', options_ImpactEscape, prms);
    
    longerTrajectories = [longerTrajectories; X_BCR_n_longerTraj(:,1:2); NaN(1,2)];
end

figure; hold all
plot(longerTrajectories(:,1),longerTrajectories(:,2),'m','linewidth',1)
p1 = plot(X_ref_BCR_n(:,1), X_ref_BCR_n(:,2),'k','linewidth',2);
plotCR3BP_Neck(secondary,L123,JC_scInitial,600,200,colors.std.black,1.5)
plotBody2(secondary.R_n, [1-secondary.MR,0,0],colors.std.white,colors.std.blue,2,0)
PlotBoi2('$X_n$','$Y_n$',20,'LaTex')
axis equal
legend([p1],'Reference')

% ========================================================================
%%% Closeout
% ========================================================================
distFig('Screen','Primary')
fh = getLocalFunctions;
fprintf('runtime: %1.1f min\n',toc(ticWhole)/60)


% ========================================================================
% ========================================================================
%%% Functions
% ========================================================================
% ========================================================================
function fh = getLocalFunctions
fh = localfunctions;
end

function getColorBarNOW(col1, col2, burnEfficiencies_vec)
cbar1 = colorbar;
colormap(colorScale([col1; col2],100))
cbar1.FontName     = 'Arial';
cbar1.FontSize     = 10;
cbar1.Ticks        = [0, 0.25, 0.5, 0.75, 1];
if min(burnEfficiencies_vec) < 0
    cbar1.TickLabels = num2cell([min(burnEfficiencies_vec), 0.25, 0.5, 0.75, 1]);
else
    cbar1.TickLabels = num2cell([0, 0.25, 0.5, 0.75, 1]);
end
cbar1.Label.String = {'Burn Efficiency'};
cbar1.Label.Rotation = 0; % 0 = horizontal, 90 = vertical
cbar1.Label.Position = [0.7, 1.05, 0];

end












