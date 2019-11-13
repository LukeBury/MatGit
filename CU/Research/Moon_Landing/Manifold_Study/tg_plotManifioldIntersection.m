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
%%% General data on solar system bodies
bodies = getBodyData(mbinPath);

%%% Color options/schemes
colors = get_colors();

%%% Periodic orbit ICs
PO_ICs = get_PO_ICs();

% ========================================================================
%%% Run Switches
% ========================================================================
run_p_unstable = 1;
run_p_stable   = 1;
run_m_unstable = 1;
run_m_stable   = 1;

% ========================================================================
%%% Setup
% ========================================================================
% -------------------------------------------------
%%% Choose data
% -------------------------------------------------
% family = 'Jupiter_Europa.CR3BP.L2_EasternAxial.txt';
family = 'Jupiter_Europa.CR3BP.L2_Lyapunov.txt';
% family = 'Jupiter_Europa.CR3BP.L2_NHalo.txt';
% family = 'Jupiter_Europa.CR3BP.L2_SHalo.txt';
% family = 'Jupiter_Europa.CR3BP.L2_Vertical.txt';
% family = 'Jupiter_Europa.CR3BP.L2_WesternAxial.txt';
% family = 'Saturn_Enceladus.CR3BP.L1_Lyapunov.txt';
% family = 'Saturn_Enceladus.CR3BP.L1_NHalo.txt';
% family = 'Saturn_Enceladus.CR3BP.L1_SHalo.txt';
% family = 'Saturn_Enceladus.CR3BP.L2_EasternAxial.txt';
% family = 'Saturn_Enceladus.CR3BP.L2_Lyapunov.txt';
% family = 'Saturn_Enceladus.CR3BP.L2_NHalo.txt';
% family = 'Saturn_Enceladus.CR3BP.L2_SHalo.txt';
% family = 'Saturn_Enceladus.CR3BP.L2_Vertical.txt';
% family = 'Saturn_Enceladus.CR3BP.L2_WesternAxial.txt';
% family = 'Saturn_Enceladus.CR3BP.SLeadingSneakerHeel.txt';
% family = 'Saturn_Enceladus.CR3BP.SLeadingSneakerToe.txt';
% family = 'Saturn_Enceladus.CR3BP.unlabeled1LeadingS.txt';
% family = 'Saturn_Titan.CR3BP.L2_Lyapunov.txt';
% family = 'Saturn_Titan.CR3BP.L2_Vertical.txt';

%%% Path from mbin to data
dataPathFromMBin = '/Data/InitialConditions/PO_Families/';

%%% PO data file
PO_datafile = [mbinPath, dataPathFromMBin, family];

% -------------------------------------------------
%%% Set up the system
% -------------------------------------------------
%%% Dynamically assign the primary/secondary
[primary, secondary] = getPlanetAndMoon(family,bodies);

%%% Normalizing constants
rNorm = secondary.a;         % n <-> km
tNorm = 1/secondary.meanMot; % n <-> sec
vNorm = rNorm / tNorm;       % n <-> km/sec

%%% prms for integration
prms.u    = secondary.MR;
prms.R2 = secondary.R_n;

%%% Equillibrium Points
rLPs_n = EquilibriumPoints(prms.u);
prms.L1x  = rLPs_n(1,1);

% -------------------------------------------------
%%% Integration options
% -------------------------------------------------
tol            = 1e-13;
options        = odeset('RelTol',tol,'AbsTol',tol);
% options_impact = odeset('Events',@event_ImpactorL1Escape_CR3Bn,'RelTol',tol,'AbsTol',tol);
options_impact = odeset('Events',@event_Impact_CR3Bn,'RelTol',tol,'AbsTol',tol);

% -------------------------------------------------
%%% Load data
% -------------------------------------------------
%%% Load the data file
PO_data = dlmread(PO_datafile,',',1,0);

%%% Number of ICs
n_POs = size(PO_data,1);

%%% Set column specifiers
c_X0   = 1:6;
c_Tp   = 7;
c_JC   = 8;
c_L2EV = 9;

% --------------------------
% PO Options
% --------------------------
%%% Index of PO from data to study
% PO_index = 40; % 'JupiterEuropa.CR3BP.L2_NHalo.txt'
PO_index = 10;

%%% Number of manifolds to integrate from PO
n_manifolds = 100;

%%% Step from PO to manifold
manifoldStep = [1, 0, 0, 0, 0, 0] .* 1e-9;
% manifoldStep = [1, 0, 0, 0, 0, 0] .* 1e-11;
warning('Test stepping in different directions')

%%% Manifold integration time
tfn_manifold = 2.3*pi;

% --------------------------
% Set color scheme for impact angles
% --------------------------
%%% Colors for impact angles of 0 and 90 degrees
color_0  = colors.std.cyan;
color_90 = colors.std.mag;

%%% Number of colors in scale
n_impactAngleColors = 4;

%%% Get nx3 matrix of colors for impact angles
impactAngleColors = colorScale([color_0; color_90], n_impactAngleColors);

% ========================================================================
%%% Integrate PO to find manifold ICs
% ========================================================================
%%% Get X0 and T of current PO
X0n_PO = PO_data(PO_index,c_X0)';
Tpn_PO = PO_data(PO_index,c_Tp);

%%% Get Jacobi constant of this IC - find flythrough and landing velocities
JC_scInitial = getJacobiConstant_ZH(X0n_PO(1:6)', prms);
L2FlythroughVelocity_mps = JC_2_L2FlyoverVelocity(JC_scInitial, prms, rLPs_n(2,:), vNorm);
JC_L2 = getJacobiConstant_ZH([rLPs_n(2,:), 0, 0, 0], prms);
minLandingVelocity_mps = JC_2_approxLandingVelocity(JC_L2, prms, prms.R2, vNorm);
approxLandingVelocity_mps = JC_2_approxLandingVelocity(JC_scInitial, prms, prms.R2, vNorm);
landingVelocityPenalty_mps = approxLandingVelocity_mps - minLandingVelocity_mps;
fprintf('L2 flythrough velocity:   %1.0f m/s\n', L2FlythroughVelocity_mps)
fprintf('Landing velocity penalty: %1.0f m/s\n', landingVelocityPenalty_mps)

%%% Integrate
[test, Xn_PO] = ode113(@Int_CR3Bn, linspace(0,Tpn_PO,n_manifolds+1), X0n_PO, options, prms);

%%% Preallocate structure for holding manifold ICs
X0n_manifolds_p = cell(n_manifolds,1);
X0n_manifolds_m = cell(n_manifolds,1);

%%% Create initial conditions of manifolds
for kk = 1:n_manifolds
    %%% Grab IC of current position on PO
    X0n_POi = Xn_PO(kk,:);
    
    %%% Store new IC for each direction of step from the PO
    X0n_manifolds_p{kk} = X0n_POi + manifoldStep;
    X0n_manifolds_m{kk} = X0n_POi - manifoldStep;
end

% ========================================================================
%%% Integrate each manifold forward and backward in time, storing
%%% trajectories
% ========================================================================

%%% Preallocate structure for holding manifold trajectories
traj_manifolds_p_unstable = cell(n_manifolds,1);
traj_manifolds_p_stable   = cell(n_manifolds,1);
traj_manifolds_m_unstable = cell(n_manifolds,1);
traj_manifolds_m_stable   = cell(n_manifolds,1);

% -------------------------------------------------
%%% Loop through initial manifold conditions and integrate
% -------------------------------------------------
for kk = 1:n_manifolds
    % --------------------------
    %%% Integrate and store trajectories
    % --------------------------
    %%% Plus-step,  unstable manifold
    if run_p_unstable
        [t_manifold, Xn_manifold] = ode113(@Int_CR3Bn, [0, tfn_manifold], X0n_manifolds_p{kk}, options_impact, prms);
        traj_manifolds_p_unstable{kk}.Xn = Xn_manifold;
        traj_manifolds_p_unstable{kk}.tn = t_manifold;
    end
    
    %%% Plus-step,  stable manifold
    if run_p_stable 
        [t_manifold, Xn_manifold] = ode113(@Int_CR3Bn, [tfn_manifold, 0], X0n_manifolds_p{kk}, options_impact, prms);
        traj_manifolds_p_stable{kk}.Xn = flipud(Xn_manifold);
        traj_manifolds_p_stable{kk}.tn = flipud(t_manifold);
    end
    
    %%% Minus-step, unstable manifold
    if run_m_unstable 
        [t_manifold, Xn_manifold] = ode113(@Int_CR3Bn, [0, tfn_manifold], X0n_manifolds_m{kk}, options_impact, prms);
        traj_manifolds_m_unstable{kk}.Xn = Xn_manifold;
        traj_manifolds_m_unstable{kk}.tn = t_manifold;
    end
    
    %%% Minus-step, unstable manifold
    if run_m_stable 
        [t_manifold, Xn_manifold] = ode113(@Int_CR3Bn, [tfn_manifold, 0], X0n_manifolds_m{kk}, options_impact, prms);
        traj_manifolds_m_stable{kk}.Xn = flipud(Xn_manifold);
        traj_manifolds_m_stable{kk}.tn = flipud(t_manifold);
    end
    % --------------------------
    %%% If unstable manifold impacted, store impact conditions
    % --------------------------
    if run_p_unstable
        [latLon_deg, impactAngle_deg, latLonHeading_hat, impactColor] = get_impactConditions(traj_manifolds_p_unstable{kk}, prms, rNorm, impactAngleColors);
        traj_manifolds_p_unstable{kk}.impactLatLon_deg   = latLon_deg;
        traj_manifolds_p_unstable{kk}.impactAngle_deg    = impactAngle_deg;
        traj_manifolds_p_unstable{kk}.latLonHeading_hat  = latLonHeading_hat;
        traj_manifolds_p_unstable{kk}.impactColor        = impactColor;
    end
    
    if run_m_unstable 
        [latLon_deg, impactAngle_deg, latLonHeading_hat, impactColor] = get_impactConditions(traj_manifolds_m_unstable{kk}, prms, rNorm, impactAngleColors);
        traj_manifolds_m_unstable{kk}.impactLatLon_deg   = latLon_deg;
        traj_manifolds_m_unstable{kk}.impactAngle_deg    = impactAngle_deg;
        traj_manifolds_m_unstable{kk}.latLonHeading_hat  = latLonHeading_hat;
        traj_manifolds_m_unstable{kk}.impactColor        = impactColor;
    end
end

% % % % % % % ========================================================================
% % % % % % %%% Post-processing the trajectories
% % % % % % % ========================================================================
% % % % % % % -------------------------------------------------
% % % % % % %%% Determine impact conditions
% % % % % % % -------------------------------------------------
% % % % % % for kk = 1:n_manifolds
% % % % % %     d
% % % % % % end
% % % % % % 
% % % % % % 
% % % % % % if abs(norm(X_eventImpact(1:3)-[1-MR,0,0])-R2_n) < 1e-9 % impact

% ========================================================================
%%% Plotting
% ========================================================================
% -------------------------------------------------
%%% Plot system with manifolds
% -------------------------------------------------
figure(1); hold all
plotSecondary(secondary)
% plotCR3BP_Neck(secondary, rLPs_n(1:3,:), prms, JC_scInitial, 600, 200, colors.std.black, 1.5)
if n_manifolds == 1
    plot3(Xn_PO(1,1),Xn_PO(1,2),Xn_PO(1,3),'ko','linewidth',2,'markersize',6,'markerfacecolor',colors.std.grey,'markeredgecolor',colors.std.black)
else
    plot3(Xn_PO(:,1),Xn_PO(:,2),Xn_PO(:,3),'ko','linewidth',2,'markersize',6,'markerfacecolor',colors.std.grey,'markeredgecolor',colors.std.black)
end

for kk = 1:n_manifolds
    if run_p_unstable
        plot3(traj_manifolds_p_unstable{kk}.Xn(:,1), traj_manifolds_p_unstable{kk}.Xn(:,2), traj_manifolds_p_unstable{kk}.Xn(:,3), 'r');
    end
    
    if run_m_unstable
        plot3(traj_manifolds_m_unstable{kk}.Xn(:,1), traj_manifolds_m_unstable{kk}.Xn(:,2), traj_manifolds_m_unstable{kk}.Xn(:,3), 'r', 'linewidth', 1);
    end
    
    if run_p_stable
        plot3(traj_manifolds_p_stable{kk}.Xn(:,1),   traj_manifolds_p_stable{kk}.Xn(:,2),   traj_manifolds_p_stable{kk}.Xn(:,3), 'g');
    end
    
    if run_m_stable
    plot3(traj_manifolds_m_stable{kk}.Xn(:,1), traj_manifolds_m_stable{kk}.Xn(:,2), traj_manifolds_m_stable{kk}.Xn(:,3), 'g');
    end

end

axis equal
view(0, 90)
xlim([0.975, 1.035])
ylim([-1, 1].*0.02)
% PlotBoi3('$X_n$', '$Y_n$', '$Z_n$', 18, 'LaTex')
PlotBoi3_CR3Bn(20)

% -------------------------------------------------
%%% Plot impact map
% -------------------------------------------------
figure(2); hold all
for kk = 1:n_manifolds
% %     if isempty(traj_manifolds_p_unstable{kk}.impactLatLon_deg) == 0
% %         plot(traj_manifolds_p_unstable{kk}.impactLatLon_deg(2), traj_manifolds_p_unstable{kk}.impactLatLon_deg(1),...
% %             'o','markersize',8,'markeredgecolor',colors.std.red,'markerfacecolor',colors.std.ltred)
% %     end
    
    if run_m_unstable
        if isempty(traj_manifolds_m_unstable{kk}.impactLatLon_deg) == 0
            plot(traj_manifolds_m_unstable{kk}.impactLatLon_deg(2), traj_manifolds_m_unstable{kk}.impactLatLon_deg(1),...
                'o','markersize',8,'markeredgecolor',colors.std.red,'markerfacecolor',colors.std.ltred)

            quiver(traj_manifolds_m_unstable{kk}.impactLatLon_deg(2), traj_manifolds_m_unstable{kk}.impactLatLon_deg(1), ...
                traj_manifolds_m_unstable{kk}.latLonHeading_hat(2), traj_manifolds_m_unstable{kk}.latLonHeading_hat(1),...
                10, 'color', traj_manifolds_m_unstable{kk}.impactColor,'linewidth',3)
        end
        
        % --------------------------
        %%% Colorbar
        % --------------------------
        cbar1 = colorbar;
        colormap(impactAngleColors)
        cbar1.FontName     = 'Arial';
        cbar1.FontSize     = 10;
        cbar1.Ticks        = [0, 0.33, 0.66, 1];
        cbar1.TickLabels   = num2cell([0, 30, 60, 90]);
        cbar1.Label.String = {'Impact Angle, deg'};
        cbar1.Label.Rotation = 0; % 0 = horizontal, 90 = vertical
        cbar1.Label.Position = [0.7, 1.05, 0];

    end
    
end

xlim([-1, 1] .* 180)
ylim([-1, 1] .* 90)
PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 18, 'LaTex')


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
fprintf('Elapsed time: %1.4f seconds\n',tocWhole)



% ========================================================================
%%% Functions
% ========================================================================
function [primary, secondary] = getPlanetAndMoon(family,bodies)
%%% Assign primary
if contains(lower(family),'earth')
    primary = bodies.earth;
elseif contains(lower(family),'jupiter')
    primary = bodies.jupiter;
elseif contains(lower(family),'saturn')
    primary = bodies.saturn;
elseif contains(lower(family),'uranus')
    primary = bodies.uranus;
elseif contains(lower(family),'neptune')
    primary = bodies.neptune;
else
    error('Primary body data not found')
end

%%% Assign secondary
if contains(lower(family),'moon')
    secondary = bodies.moon;
elseif contains(lower(family),'io')
    secondary = bodies.io;
elseif contains(lower(family),'europa')
    secondary = bodies.europa;
elseif contains(lower(family),'ganymede')
    secondary = bodies.ganymede;
elseif contains(lower(family),'callisto')
    secondary = bodies.callisto;
elseif contains(lower(family),'enceladus')
    secondary = bodies.enceladus;
elseif contains(lower(family),'titan')
    secondary = bodies.titan;
elseif contains(lower(family),'triton')
    secondary = bodies.triton;
else
    error('Secondary body data not found')
end
end


function [latLon_deg, impactAngle_deg, latLonHeading_hat, impactColor] = get_impactConditions(traj, prms, rNorm, impactAngleColors)
if abs(norm(traj.Xn(end,1:3)-[1-prms.u,0,0])-prms.R2) < ((1e-3)/rNorm) % If trajectory impacted 
    %%% Determine and store impact lat/lon
    [lat_deg, lon_deg] = BCR2latlon(traj.Xn(end,1:3), 'secondary', prms.u);
    latLon_deg = [lat_deg, lon_deg];

    %%% Determine and store impact angle
    % Creating SCR position vector
    rImpact_SCR_n = traj.Xn(end,1:3) - [1-prms.u,0,0];
    rImpact_SCR_n = rImpact_SCR_n./norm(rImpact_SCR_n);

    % Velocity unit vector at impact
    vHatImpact_n = traj.Xn(end,4:6)./norm(traj.Xn(end,4:6));

    % Angle between velocity and surface at impact
    [impactAngle_deg] = calcImpactAngle(rImpact_SCR_n,vHatImpact_n,'degrees');
    
    %%% Heading at impact
    r_SCR = traj.Xn(end,1:3) - [1-prms.u, 0, 0];
    rHat_SCR = r_SCR ./ norm(r_SCR);
    headingVec = traj.Xn(end,4:6) - rHat_SCR * dot(traj.Xn(end,4:6),rHat_SCR);
    headingHat = headingVec ./ norm(headingVec);
    refPoint_SCR = r_SCR + headingHat .* ((1e-3)/rNorm);
    refPoint_BCR = refPoint_SCR + [1-prms.u, 0, 0];
    [latHeadingRef_deg, lonHeadingRef_deg] = BCR2latlon(refPoint_BCR, 'secondary', prms.u);
    latLonHeading = [latHeadingRef_deg, lonHeadingRef_deg] - latLon_deg;
    latLonHeading_hat = latLonHeading ./ norm(latLonHeading);
    
    %%% Color for impact angle
    n_impactAngleColorScale = size(impactAngleColors, 1);
    colorIndex_unrounded = interp1(linspace(0,90,n_impactAngleColorScale),linspace(1,n_impactAngleColorScale,n_impactAngleColorScale), impactAngle_deg);
    impactColor = impactAngleColors(round(colorIndex_unrounded), :);
else
    latLon_deg        = [];
    impactAngle_deg   = [];
    latLonHeading_hat = [];
    impactColor       = [];
end
end









