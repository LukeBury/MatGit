% ========================================================================
%%% Description
% ========================================================================
% This script is designed to plot the results generated in in
% getShallowImpactAngles.m. This script will generate a lat/lon plot of all
% the shallow impact-angle locations achieved by manifolds of a certain
% family of periodic orbits. Results can be color coded by energy

% Created: 11/13/19
% Author : Luke Bury, luke.bury@colorado.edu
% ========================================================================
%%% Initialization
% ========================================================================
clear
clc
close all
mbinPath    = '~/CU_Google_Drive/Documents/MatGit/mbin';
% dataBinpath = '~/CU_Google_Drive/Documents/MatGit/MatlabOutputs/ShallowImpactDataSets';
dataBinpath = '~/CU_Google_Drive/Documents/MatGit/MatlabOutputs/ManifoldShallowImpactDataSets';
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
useImpactSpeedColor = 1;

plot_shallowImpactTrajectories  = 0;

plot_latitudeReflection = 1;

% ========================================================================
%%% Setup
% ========================================================================
% -------------------------------------------------
%%% Preferences
% -------------------------------------------------
headingLength = 8;

% -------------------------------------------------
%%% Choose data file
% -------------------------------------------------
% dataFile = 'shallowImpacts.F.Earth_Moon.CR3BP.L1_Lyapunov.nodes2000.txt';
% dataFile = 'shallowImpacts.F.Earth_Moon.CR3BP.L1_Vertical.nodes2000.txt';
% dataFile = 'shallowImpacts.F.Earth_Moon.CR3BP.L1_SHalo.nodes2000.txt';
% dataFile = 'shallowImpacts.F.Earth_Moon.CR3BP.L1_NHalo.nodes2000.txt';
% dataFile = 'shallowImpacts.F.Earth_Moon.CR3BP.L1_All.nodes2000.txt'; 

% dataFile = 'shallowImpacts.F.Jupiter_Europa.CR3BP.L2_Lyapunov.nodes2000.txt';
% dataFile = 'shallowImpacts.F.Jupiter_Europa.CR3BP.L2_Vertical.nodes2000.txt';
% dataFile = 'shallowImpacts.F.Jupiter_Europa.CR3BP.L2_SHalo.nodes2000.txt';
% dataFile = 'shallowImpacts.F.Jupiter_Europa.CR3BP.L2_NHalo.nodes2000.txt';
% dataFile = 'shallowImpacts.F.Jupiter_Europa.CR3BP.L2_All.nodes2000.txt';

% dataFile = 'shallowImpacts.F.Jupiter_Ganymede.CR3BP.L2_Lyapunov.nodes2000.txt';
% dataFile = 'shallowImpacts.F.Jupiter_Ganymede.CR3BP.L2_Vertical.nodes2000.txt';
% dataFile = 'shallowImpacts.F.Jupiter_Ganymede.CR3BP.L2_SHalo.nodes2000.txt';
% dataFile = 'shallowImpacts.F.Jupiter_Ganymede.CR3BP.L2_NHalo.nodes2000.txt';
% dataFile = 'shallowImpacts.F.Jupiter_Ganymede.CR3BP.L2_All.nodes2000.txt';

% dataFile = 'shallowImpacts.F.Saturn_Enceladus.CR3BP.L2_Lyapunov.nodes2000.txt';
% dataFile = 'shallowImpacts.F.Saturn_Enceladus.CR3BP.L2_Vertical.nodes2000.txt'; % n_shallowImpacts = 197;
% dataFile = 'shallowImpacts.F.Saturn_Enceladus.CR3BP.L2_SHalo.nodes2000.txt';
% dataFile = 'shallowImpacts.F.Saturn_Enceladus.CR3BP.L2_NHalo.nodes2000.txt';
% dataFile = 'shallowImpacts.F.Saturn_Enceladus.CR3BP.L2_All.nodes2000.txt';

% dataFile = 'shallowImpacts.F.Neptune_Triton.CR3BP.L2_Lyapunov.nodes2000.txt';
% dataFile = 'shallowImpacts.F.Neptune_Triton.CR3BP.L2_Vertical.nodes2000.txt';
% dataFile = 'shallowImpacts.F.Neptune_Triton.CR3BP.L2_SHalo.nodes2000.txt';
% dataFile = 'shallowImpacts.F.Neptune_Triton.CR3BP.L2_NHalo.nodes2000.txt';
% dataFile = 'shallowImpacts.F.Neptune_Triton.CR3BP.L2_All.nodes2000.txt';



dataFile = 'shallowImpacts.M.Jupiter_Europa.CR3BP.L2_SHalo.nodes2000.txt';
% dataFile = 'shallowImpacts.M.Jupiter_Europa.CR3BP.L2_SHalo.nodes2000_EnceladusImposter.txt';
% dataFile = 'shallowImpacts.M.Jupiter_Europa.CR3BP.L2_SHalo.nodes2000_MoonImposter.txt';



%%% Create path to dataFile
dataFilePath = [dataBinpath, '/',dataFile];

% -------------------------------------------------
%%% Set up system
% -------------------------------------------------
%%% Set primary & secondary
[primary, secondary] = assignPrimaryAndSecondary_CR3BP(dataFile, bodies);

%%% Normalizing constants
rNorm = secondary.a;         % n <-> km
tNorm = 1/secondary.meanMot; % n <-> sec
vNorm = rNorm / tNorm;       % n <-> km/sec

%%% Setting parameter structure
prms.u     = secondary.MR;
prms.rNorm = rNorm;
prms.R1    = primary.R / rNorm;
prms.R2    = secondary.R_n;
prms.n = 1;

if contains(dataFile,'.CR3BP_J2pJ4pJ6pJ2s.')
    prms.J2p = primary.J2; prms.J4p = primary.J4; prms.J6p = primary.J6; prms.J2s = secondary.J2;
    warning('Need to calculate prms.n')
end

% -------------------------------------------------
%%% Colors for impact speed
% -------------------------------------------------
color1 = colors.grn2;
color2 = colors.blue2;
% color1 = colors.blue2;
% color2 = colors.red3;

colorDiff = color2 - color1;
colorBarColors = colorScale([color1; color2],100);
% -------------------------------------------------
%%% Load PO family data
% -------------------------------------------------
%%% Load shallow angle impact data
shallowAngleImpactData = dlmread(dataFilePath,',',1,0);

% 989
% shallowAngleImpactData = shallowAngleImpactData(1:300, :);

%%% Column specifiers
c_x0_n              = 1;
c_y0_n              = 2;
c_z0_n              = 3;
c_xd0_n             = 4;
c_yd0_n             = 5;
c_zd0_n             = 6;
c_Tf_n              = 7;
c_lat               = 8;
c_lon               = 9;
c_impactAngle_deg   = 10;
c_impactHeading_lat = 11;
c_impactHeading_lon = 12;
c_JC                = 13;
c_landingVel_mps    = 14;

% -------------------------------------------------
%%% Integrator options
% -------------------------------------------------
tol                      = 1e-13;
options                  = odeset('RelTol',tol,'AbsTol',tol);

% ========================================================================
%%% Make plot of shallow-angle impacts
% ========================================================================
%%% Get number of impacts
n_shallowImpacts = size(shallowAngleImpactData,1);
% 989
% n_shallowImpacts = 197;

%%% Get min and max landing speeds
maxLandingSpeed_mps = max(shallowAngleImpactData(1:n_shallowImpacts, c_landingVel_mps));
minLandingSpeed_mps = min(shallowAngleImpactData(1:n_shallowImpacts, c_landingVel_mps));


%%% Generate figure
% figure('position', [156 385 560 420]); hold all
figure; hold all
xlim([-1, 1] .* 185)
ylim([-1, 1] .* 95)
PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')

if useImpactSpeedColor
%     cbar1 = colorbar;
%     cbar1.FontName     = 'Arial';
%     cbar1.FontSize     = 10;
%     cbar1.Ticks        = linspace(1e-5, 1, 5);
%     cbar1.TickLabels   = num2cell(round(linspace(minLandingSpeed_mps, maxLandingSpeed_mps, 5)));
%     cbar1.Label.String = {'Impact Speed, m/s'};
%     cbar1.Label.Rotation = 0; % 0 = horizontal, 90 = vertical
%     cbar1.Label.Position = [0.7, 1.05, 0];
%     colormap(colorBarColors)

    cbar1 = colorbar;
    cbar1.FontName     = 'Arial';
    cbar1.FontSize     = 14;
    cbar1.Ticks        = linspace(1e-5, 1, 5);
    cbar1.TickLabels   = num2cell(round(linspace(minLandingSpeed_mps, maxLandingSpeed_mps, 5)));
    cbar1.Label.String = {'Impact Speed, m/s'};
    cbar1.Label.Rotation = 0; % 0 = horizontal, 90 = vertical
    cbar1.Label.Position = [0.7, 1.07, 0];
    colormap(colorBarColors)
end
% 
% for impact_i = 1:n_shallowImpacts
%     plot(shallowAngleImpactData(impact_i, c_lon), shallowAngleImpactData(impact_i, c_lat), '.','markersize',17,...
%         'markeredgecolor',colors.std.black,'markerfacecolor',colors.std.grey)
%     
%     quiver(shallowAngleImpactData(impact_i, c_lon), shallowAngleImpactData(impact_i, c_lat),...
%         shallowAngleImpactData(impact_i, c_impactHeading_x), shallowAngleImpactData(impact_i, c_impactHeading_y),...
%         10,'color',colors.std.cyan,'linewidth',3)
% 
% end

for impact_i = 1:n_shallowImpacts
    % --------------------------
    % For plotting the heading directions, since the lat/lon plot isn't
    % equal, we need to artificially increase the magnitude as headings
    % become closer to purely horizontal (since '1' in the longitude-axis
    % is half as long as '1' in the latitude axis). To do this we calculate
    % the angle between the heading and [1,0], and use the cosine of that
    % angle to increase the heading's magnitude
    % --------------------------
    dx = shallowAngleImpactData(impact_i,c_impactHeading_lon);
    dy = shallowAngleImpactData(impact_i,c_impactHeading_lat);
    
    headingVec = [dx, dy]./norm([dx, dy]) .*headingLength;
    angle = acos((dot(headingVec, [1,0])./norm(headingVec)));
    scaling = 1 + abs(cos(angle));
    
    headingVec = headingVec.*scaling;
    
    
    % --------------------------
    % Plot the heading line
    % --------------------------
    xs = [shallowAngleImpactData(impact_i,c_lon), shallowAngleImpactData(impact_i,c_lon) + headingVec(1)];
    ys = [shallowAngleImpactData(impact_i,c_lat), shallowAngleImpactData(impact_i,c_lat) + headingVec(2)];
    if useImpactSpeedColor
        landingSpeedFrac = (shallowAngleImpactData(impact_i, c_landingVel_mps) - minLandingSpeed_mps) / (maxLandingSpeed_mps - minLandingSpeed_mps);
        impactSpeedColor_i = color1 + colorDiff.*landingSpeedFrac;
        plot(xs, ys, 'color', impactSpeedColor_i, 'linewidth', 1.5)
        
        if plot_latitudeReflection
            plot(xs, -ys, 'color', impactSpeedColor_i, 'linewidth', 1.5)
        end
    else
        plot(xs, ys, 'color', colors.ltred, 'linewidth', 1.5)
        
        if plot_latitudeReflection
            plot(xs, -ys, 'color', colors.ltred, 'linewidth', 1.5)
        end
    end
    
    
    
    
    % --------------------------
    % Plot the lat/lon impact point (iteratively)
    % --------------------------
    if (n_shallowImpacts ~= size(shallowAngleImpactData,1)) || (useImpactSpeedColor == 1)
        if useImpactSpeedColor
            plot3(shallowAngleImpactData(impact_i, c_lon), shallowAngleImpactData(impact_i, c_lat), impact_i, '.','markersize',8,...
                'markeredgecolor',colors.black,'markerfacecolor',impactSpeedColor_i)
            if plot_latitudeReflection
                plot3(shallowAngleImpactData(impact_i, c_lon), -shallowAngleImpactData(impact_i, c_lat), impact_i, '.','markersize',8,...
                'markeredgecolor',colors.black,'markerfacecolor',impactSpeedColor_i)
            end
        else
            plot3(shallowAngleImpactData(impact_i, c_lon), shallowAngleImpactData(impact_i, c_lat), impact_i, '.','markersize',8,...
                'markeredgecolor',colors.red,'markerfacecolor',colors.red)
            if plot_latitudeReflection
                plot3(shallowAngleImpactData(impact_i, c_lon), -shallowAngleImpactData(impact_i, c_lat), impact_i, '.','markersize',8,...
                'markeredgecolor',colors.red,'markerfacecolor',colors.red)
            end
        end
        
    end
    
end

% --------------------------
% Plot the lat/lon impact point (all at once)
% --------------------------
if n_shallowImpacts == size(shallowAngleImpactData,1) && (useImpactSpeedColor == 0)
    plot(shallowAngleImpactData(:, c_lon), shallowAngleImpactData(:, c_lat), '.','markersize',8,...
        'markeredgecolor',colors.red,'markerfacecolor',colors.red)
end

% ========================================================================
%%% Plot trajectories
% ========================================================================
if plot_shallowImpactTrajectories
    trajectories = cell(n_shallowImpacts,1);
    
    X0s_n = shallowAngleImpactData(:,c_x0_n:c_zd0_n);
    Tfs_n = shallowAngleImpactData(:,c_Tf_n);
    parfor impact_i = 1:n_shallowImpacts
        X0_n = X0s_n(impact_i,:)';
        Tf_n = Tfs_n(impact_i);

        [~, X_n] = ode113(@Int_CR3Bn, [0, Tf_n], X0_n, options, prms);

        trajectories{impact_i} = X_n;
    end

    figure('position', [717 385 560 420]); hold all
    PlotBoi3_CR3Bn(20)
    axis equal
    plotSecondary(secondary)
    for impact_i = 1:n_shallowImpacts
        plot3(trajectories{impact_i}(:,1), trajectories{impact_i}(:,2), trajectories{impact_i}(:,3), 'r', 'linewidth', 0.5);
    end

%     for impact_i = 1:n_shallowImpacts
%         if impact_i <= 197
%             p1 = plot3(trajectories{impact_i}(:,1), trajectories{impact_i}(:,2), trajectories{impact_i}(:,3), 'r', 'linewidth', 0.5);
%         else
%             p2 = plot3(trajectories{impact_i}(:,1), trajectories{impact_i}(:,2), trajectories{impact_i}(:,3), 'color',colors.drkgrn, 'linewidth', 0.5);
%         end
%     end
%     [hleg, hobj, hout, mout] = legend([p1, p2], 'Low Energy Impacts', 'High Energy Impacts');
%     set(hobj,'linewidth',2);
%     set(hobj,'FontSize',12);
end



% % 385, 994
% % 
% % if plot_shallowImpactTrajectories
% %     trajectories = cell(n_shallowImpacts,1);
% %     
% %     X0s_n = shallowAngleImpactData(:,c_x0_n:c_zd0_n);
% %     Tfs_n = shallowAngleImpactData(:,c_Tf_n);
% %     for impact_i = [385]
% %         X0_n = X0s_n(impact_i,:)';
% %         Tf_n = Tfs_n(impact_i);
% % 
% %         [~, X_n] = ode113(@Int_CR3Bn, [0, Tf_n], X0_n, options, prms);
% % 
% %         trajectories{impact_i} = X_n;
% %     end
% % 
% %     figure('position', [717 385 560 420]); hold all
% %     PlotBoi3_CR3Bn(20)
% %     axis equal
% %     plotSecondary(secondary)
% %     for impact_i = [385]
% %         plot3(trajectories{impact_i}(:,1), trajectories{impact_i}(:,2), trajectories{impact_i}(:,3), 'r', 'linewidth', 1.5)
% %     end
% % 
% % end



if contains(dataFile, 'Lyapunov')
    %%% Generate figure
    figure; hold all
    xlim([-1, 1] .* 185)
    % ylim([-1, 1] .* 95)
    PlotBoi2('Longitude, $^\circ$', 'Landing Speed, $m/s$', 18, 'LaTex')

    for impact_i = 1:n_shallowImpacts
        % --------------------------
        % For plotting the heading directions, since the lat/lon plot isn't
        % equal, we need to artificially increase the magnitude as headings
        % become closer to purely horizontal (since '1' in the longitude-axis
        % is half as long as '1' in the latitude axis). To do this we calculate
        % the angle between the heading and [1,0], and use the cosine of that
        % angle to increase the heading's magnitude
        % --------------------------
        dx = shallowAngleImpactData(impact_i,c_impactHeading_lon);
        dy = shallowAngleImpactData(impact_i,c_impactHeading_lat);

        headingVec = [dx, dy]./norm([dx, dy]) .*headingLength;
        angle = acos((dot(headingVec, [1,0])./norm(headingVec)));
        scaling = 1 + abs(cos(angle));

        headingVec = headingVec.*scaling;
        
        
        
        
        
        
        
        % --------------------------
        % Plot the heading line
        % --------------------------
    %     if (shallowAngleImpactData(impact_i,c_lon) < -160.139) && (shallowAngleImpactData(impact_i,c_lon) > -160.141)
    %         impact_i
    %     end
        xs = [shallowAngleImpactData(impact_i,c_lon), shallowAngleImpactData(impact_i,c_lon) + headingVec(1)];
        ys = [shallowAngleImpactData(impact_i,c_landingVel_mps), shallowAngleImpactData(impact_i,c_landingVel_mps)];
        
        if headingVec(1) > 0
            p1 = plot(xs, ys, 'color', colors.mag, 'linewidth', 1.5);
        elseif headingVec(1) < 0
            p2 = plot(xs, ys, 'color', colors.drkgrn2, 'linewidth', 1.5);
        end
        
        

        % --------------------------
        % Plot the lat/lon impact point (iteratively)
        % --------------------------

        plot(shallowAngleImpactData(impact_i, c_lon), shallowAngleImpactData(impact_i, c_landingVel_mps), '.','markersize',8,...
            'markeredgecolor',colors.black,'markerfacecolor',colors.black)

    end
    legend([p1, p2],'East','West')
end



% ========================================================================
%%% Closeout
% ========================================================================
tocWhole = toc(ticWhole);
fprintf('Elapsed time: %1.4f seconds\n',tocWhole)
















