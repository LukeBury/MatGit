clear
clc
close all
addpath(genpath('/Users/lukebury/Documents/MATLAB/mbin'))

% ========================================================================
%%% Run Switches
% ========================================================================
% ----------------------
% Functions - Single Runs
% ----------------------
%%% Primary Functions
singleSimulation = 0;
        %%% Secondary Functions
        firstPeriapsisAnalysis_sing = 0;
        
        %%% Plots and Extras
        latlonPlot     = 1;
        plotBCR_n      = 1;
        plotSCI_n      = 0;
        jacobiAnalysis = 0;

% ----------------------
% Functions - Multiple Runs
% ----------------------
%%% Primary Functions
multSimulation = 1;
        %%% Secondary Functions
        firstPeriapsisAnalysis_mult = 1;
        
        %%% Plots and Extras
        plotBCR_n_mult  = 1;
        fpaVSspeed_mult = 1;
        latlonPlot_mult = 1;

% ========================================================================
%%% Importing Data
% ========================================================================
%%% General data on solar system bodies
bodies = getBodyData();

%%% Color options/schemes
colors = get_colors();

% ========================================================================
%%% Constant / Conditions
% ========================================================================
%%% Selecting moon
secondary = bodies.europa;
%%% for R_n:
% enceladus ~ mimas
% oberon ~ callisto

%%% Normalizing constants
rNorm = secondary.a; % km
tNorm = 1/secondary.meanMot; % sec
vNorm = rNorm / tNorm; % km/sec

%%% Acquire Lagrange points
L_points_n = EquilibriumPoints(secondary.MR);
L1_n = L_points_n(1,:); % [1x3] normalized barycentric coordinates
L2_n = L_points_n(2,:); % [1x3] normalized barycentric coordinates


% ========================================================================
%%% Single Run Simulation
% ========================================================================
if singleSimulation == 1
    %%% Nudges from equillibrium point
    % Enceladus North Pole
%     positionNudge = [0, 0, 0]./rNorm; % enter in km
%     velocityNudge = R3([-80, 0, 55],40*pi/180)./1000./vNorm; % enter in m/s
    positionNudge = [0, 10, 0]./rNorm; % enter in km
    velocityNudge = [0, -10, 50]./1000./vNorm; % enter in m/s

    %%% Initial State
    x0_n = L2_n + positionNudge;
    xdot0_n = [0, 0, 0] + velocityNudge;
    X0_n = [x0_n xdot0_n]';

    %%% Setting time vector and normalizing 
    ti = 0; % sec
    dt = 60; % sec
    tf = 3600*24*5; % sec
    time0_n = [ti:dt:tf] ./ tNorm;

    %%% Choosing ode45 tolerance
    tol = 1e-10;
    
    %%% Setting integrator options
    if firstPeriapsisAnalysis_sing == 0 % Stop on impact
        options = odeset('Events',@impactEvent_CR3Bn,'RelTol',tol,'AbsTol',tol);
    elseif firstPeriapsisAnalysis_sing == 1 % Don't Stop on Impact
        options = odeset('RelTol',tol,'AbsTol',tol);
    end

    %%% Propagating the State
    [time_n, X_BCR_n] = ode45(@Int_CR3Bn, time0_n, X0_n, options, secondary.MR, secondary.R_n);
    
    % -----------------------------
    % Jacobi Constant Analysis
    % -----------------------------
    if jacobiAnalysis == 1
        [JCs] = JacobiConstantCalculator(secondary.MR,X_BCR_n(:,1:3),X_BCR_n(:,4:6));

        figure
        plot(percentchange(JCs))
        PlotBoi2('Time Indices','JC % Change',14)
    end

    % -----------------------------
    % Latitude & Longitude
    % -----------------------------
    %%% Acquire normalized secondary-centered-rotating coordinate (SCR_n)
    r_SCR_n = X_BCR_n(:,1:3) - [1-secondary.MR, 0, 0];

    %%% Acquire Lat/Lon of landing site on secondary
    angleUnit = 'degrees';
    [latitudes, longitudes] = ECEF2latlon(r_SCR_n,angleUnit);

    % -----------------------------
    % Closest Approach and Altitude
    % -----------------------------
    CAindex = find(rownorm(r_SCR_n) == min(rownorm(r_SCR_n)));
    altitudes = rownorm(r_SCR_n).*rNorm - secondary.R;

    % -----------------------------
    % Acquiring Normalized Inertial Results
    % -----------------------------
    [ r_SCI_n ] = r_BCR2BCI( r_SCR_n, time_n.*tNorm, secondary.meanMot );
    
    % -----------------------------
    % First Periapse Analysis
    % -----------------------------
    if firstPeriapsisAnalysis_sing == 1
        periapsisFound = 0;
        firstPeriapsis_sing = [];
        
        %%% Creating SCR positions
        r_SCR_n_sing = X_BCR_n(:,1:3) - [1-secondary.MR, 0, 0];
        
        secondaryDistances = rownorm(r_SCR_n_sing);
        for kk = 2:size(secondaryDistances,1)
            if secondaryDistances(kk-1) < secondaryDistances(kk) && norm(secondaryDistances(kk) - secondaryDistances(kk-1)) > tol
                periapsisFound = 1;
                firstPeriapsis_sing = X_BCR_n(kk,1:3);
                break
            end
        end
        if periapsisFound == 0 && secondaryDistances(1) > secondaryDistances(end)
            firstPeriapsis_sing = [firstPeriapsis_sing; X_BCR_n(end,1:3)];
        elseif periapsisFound == 0 && secondaryDistances(1) < secondaryDistances(end)
            firstPeriapsis_sing = [firstPeriapsis_sing; X_BCR_n(1,1:3)];
        end
        
        %%% Plotting
        figure; hold all
        plot3(firstPeriapsis_sing(:,1),firstPeriapsis_sing(:,2),firstPeriapsis_sing(:,3),'rx','markersize',5,'linewidth',2)
        plot3(X_BCR_n(:,1),X_BCR_n(:,2),X_BCR_n(:,3),'g','linewidth',1.0)
%         plot3(L1_n(1),L1_n(2),L1_n(3),'k>','linewidth',1,'markersize',10)
%         plot3(L2_n(1),L2_n(2),L2_n(3),'k<','linewidth',1,'markersize',10)
        plotBody3(secondary.R_n, [1-secondary.MR, 0, 0], secondary.color)
        alpha(.3)
        PlotBoi3('X_n','Y_n','Z_n',14)
        axis equal
        xlim([1-secondary.MR-5*secondary.R_n, 1-secondary.MR+5*secondary.R_n])
        ylim([-5*secondary.R_n, 5*secondary.R_n])
        zlim([-5*secondary.R_n, 5*secondary.R_n])
        view(0,90)
    end

    % ----------------------------------------------------------
    %%% Plots and Extras
    % ----------------------------------------------------------
    % -----------------------------
    % Normalized Secondary-Centric Rotating Plot
    % -----------------------------
    if plotBCR_n == 1
        runTime = time_n(end) * tNorm / 3600;
        figure; hold all
        plot3(X_BCR_n(:,1),X_BCR_n(:,2),X_BCR_n(:,3),'linewidth',2)
        plotBody3(secondary.R_n, [1-secondary.MR, 0, 0], secondary.color)
%         colors.std.turq
        PlotBoi3('X_n','Y_n','Z_n',14)
        title(sprintf('%s - BCR, Runtime: %1.1f hours', secondary.name, runTime))
        axis equal
        view(0,90)

        %%% Ending location
        if time_n(end) ~= time0_n(end) % Crash
            p1 = plot3(X_BCR_n(end,1),X_BCR_n(end,2),X_BCR_n(end,3),'rx','markersize',10,'linewidth',2);
        else
            p1 = plot3(X_BCR_n(end,1),X_BCR_n(end,2),X_BCR_n(end,3),'go','markersize',10,'linewidth',2,'markerfacecolor','g');
        end

        %%% Closest approach location
        p2 = plot3(X_BCR_n(CAindex,1),X_BCR_n(CAindex,2),X_BCR_n(CAindex,3),'^','markeredgecolor',colors.std.black,'markerfacecolor',colors.std.ltred,'markersize',10);

        %%% legend
        legend([p1 p2],'End','C/A')
    end

    % -----------------------------
    % Normalized Secondary-Centric Inertial Plot
    % -----------------------------
    if plotSCI_n == 1
        runTime = time_n(end) * tNorm / 3600;
        figure; hold all
        plot3(r_SCI_n(:,1),r_SCI_n(:,2),r_SCI_n(:,3),'linewidth',2)
        plotBody3(secondary.R_n, [0, 0, 0], secondary.color)
        PlotBoi3('X_n','Y_n','Z_n',14)
        title(sprintf('%s - SCI, Runtime: %1.1f hours',secondary.name, runTime))
        axis equal
        view(0,90)

        %%% Ending location
        if time_n(end) ~= time0_n(end) % Crash
            p1 = plot3(r_SCI_n(end,1),r_SCI_n(end,2),r_SCI_n(end,3),'rx','markersize',10,'linewidth',2);
        else
            p1 = plot3(r_SCI_n(end,1),r_SCI_n(end,2),r_SCI_n(end,3),'go','markersize',10,'linewidth',2,'markerfacecolor','g');
        end

        %%% Closest approach location
        p2 = plot3(r_SCI_n(CAindex,1),r_SCI_n(CAindex,2),r_SCI_n(CAindex,3),'^','markeredgecolor',colors.std.black,'markerfacecolor',colors.std.ltred,'markersize',10);

        %%% legend
        legend([p1 p2],'End','C/A')
    end

    % -----------------------------
    % Plot of Lat/Lon Groundtrack on Secondary
    % -----------------------------
    if latlonPlot == 1
        %%% Plotting
        figure
        % -----------
        % Altitude
        % -----------
        subplot(2,1,1); hold all
        plot(time_n.*tNorm, altitudes, 'linewidth',1.5)
        PlotBoi2('Time, sec','Altitude, km',14)

        %%% Closest approach location
        plot(time_n(CAindex).*tNorm, altitudes(CAindex),'^','markeredgecolor',colors.std.black,'markerfacecolor',colors.std.ltred,'markersize',5)

        %%% Ending location
        if time_n(end) ~= time0_n(end) % Crash
            plot(time_n(end).*tNorm,altitudes(end),'rx','markersize',6,'linewidth',2)
        else
            plot(time_n(end).*tNorm,altitudes(end),'go','markersize',6,'markerfacecolor','g','markeredgecolor','k')
        end

        % -----------
        % Lat/Lon
        % -----------
        subplot(2,1,2); hold all
        plot(longitudes,latitudes,'linewidth',1.5)
        grid on
        if strcmp(angleUnit,'degrees') == 1
            PlotBoi2('Longitude, °', 'Latitude, °', 14)
            xlim([-180 180])
            ylim([-90 90])
        elseif strcmp(angleUnit,'radians') == 1
            PlotBoi2('Longitude, rad', 'Latitude, rad', 14)
            xlim([-pi pi])
            ylim([-pi/2 pi/2])
        end

        %%% Ending location
        if time_n(end) ~= time0_n(end) % Crash
            plot(longitudes(end),latitudes(end),'rx','markersize',10,'linewidth',2)
        else
            plot(longitudes(end),latitudes(end),'go','markersize',10,'markerfacecolor','g','markeredgecolor','k')
        end

        %%% Closest approach location
        plot(longitudes(CAindex),latitudes(CAindex),'^','markeredgecolor',colors.std.black,'markerfacecolor',colors.std.ltred,'markersize',10)

    end
 
    
    
end




% ========================================================================
%%% Multiple Run Simulation
% ========================================================================
if multSimulation == 1
    %%% Setting time vector and normalizing 
%     ti_mult = 0; % sec
%     dt_mult = 60*1; % sec 
%     tf_mult = 3600*24; % sec 
%     time0_n_mult = [ti_mult:dt_mult:tf_mult] ./ tNorm;
    
    time0_n_mult = [0:0.004629904709320:6.667062781420788];
    
%     %% Set ranges for initial positions
%     xrange_mult = linspace(0, 0, 1)./rNorm; % enter in km
%     yrange_mult = linspace(0, 0, 1)./rNorm; % enter in km
%     zrange_mult = linspace(0, 0, 1)./rNorm; % enter in km
%     
%     %%% Set ranges for initial velocities
%     dxrange_mult = linspace(-10, 0, 5)./1000./vNorm; % enter in m/s
%     dyrange_mult = linspace(0, 0, 1)./1000./vNorm; % enter in m/s
%     dzrange_mult = linspace(-40, 40, 10)./1000./vNorm; % enter in m/s
    
    %%% Set ranges for initial positions
    xrange_mult = linspace(0, 0, 1); 
    yrange_mult = linspace(0, 0, 1); 
    zrange_mult = linspace(0, 0, 1); 
    
    %%% Set ranges for initial velocities
    dxrange_mult = linspace(-7.918473792209050e-04, 0, 5);
    dyrange_mult = linspace(0, 0, 1);
    dzrange_mult = linspace(-0.003167389516884, 0.003167389516884, 5);
    
    %%% Create ordered matrices of intial conditions
    positionNudgeMatrix = [];
    for xx = 1:length(xrange_mult)
        for yy = 1:length(yrange_mult)
            for zz = 1:length(zrange_mult)
                positionNudgeMatrix = [positionNudgeMatrix; xrange_mult(xx), yrange_mult(yy), zrange_mult(zz)];
            end
        end
    end
    
    velocityNudgeMatrix = [];
    for xx = 1:length(dxrange_mult)
        for yy = 1:length(dyrange_mult)
            for zz = 1:length(dzrange_mult)
                velocityNudgeMatrix = [velocityNudgeMatrix; dxrange_mult(xx), dyrange_mult(yy), dzrange_mult(zz)];
            end
        end
    end
    
    %%% Preallocating
    r_BCR_mult = [];
    endSpeeds_n_mult = [];
    endFPAs_mult = [];
    latlons_mult = [];
    firstPeriapsis_mult = [];
    anyLandings = 0;
    periapsisFound = 0;
    
    for pp = 1:size(positionNudgeMatrix,1)
        for vv = 1:size(velocityNudgeMatrix,1)
            
            %%% Adding nudges to the initial state
            x0_n_mult = L2_n + positionNudgeMatrix(pp,:);
            xdot0_n_mult = [0, 0, 0] + velocityNudgeMatrix(vv,:);
            X0_n_mult = [x0_n_mult xdot0_n_mult]';

            %%% Choosing ode45 tolerance
            tol = 1e-10;

            %%% Setting integrator options
            if firstPeriapsisAnalysis_mult == 0 % Stop on impact
                options = odeset('Events',@impactEvent_CR3Bn,'RelTol',tol,'AbsTol',tol);
            elseif firstPeriapsisAnalysis_mult == 1 % Don't Stop on Impact
                options = odeset('RelTol',tol,'AbsTol',tol);
            end

            %%% Propagating the State
            [time_n_mult, X_BCR_n_mult] = ode45(@Int_CR3Bn, time0_n_mult, X0_n_mult, options, secondary.MR, secondary.R_n);
            
            %%% Storing new state positions
            r_BCR_mult = [r_BCR_mult; X_BCR_n_mult(:,1:3); nan(1,3)];
            
            %%% Creating SCR positions
            r_SCR_n_mult = X_BCR_n_mult(:,1:3) - [1-secondary.MR, 0, 0];

            %%% Analyzing Impacts
            if time_n_mult(end) ~= time0_n_mult(end) % Crash
                anyLandings = 1;
                
                %%% Calculate end speed and flight path angle
                [ endFPA ] = rv2fpa( X_BCR_n_mult(end,1:3), X_BCR_n_mult(end,4:6)); % rad
                endSpeed_n = norm(X_BCR_n_mult(end,4:6));
                
                %%% Storing new end speed and flight path angles
                endSpeeds_n_mult = [endSpeeds_n_mult; endSpeed_n];
                endFPAs_mult = [endFPAs_mult; endFPA];

                %%% Acquire Lat/Lon of landing site
                angleUnit = 'degrees';
                [lat, lon] = ECEF2latlon(r_SCR_n_mult(end,1:3),angleUnit);
                
                %%% Storing Lat/Lon of landing site
                latlons_mult = [latlons_mult; lat, lon];
            end
            
            %%% Periapsis Analysis
            if firstPeriapsisAnalysis_mult == 1
                secondaryDistances = rownorm(r_SCR_n_mult);
                for kk = 2:size(secondaryDistances,1)
                    if secondaryDistances(kk-1) < secondaryDistances(kk) && norm(secondaryDistances(kk) - secondaryDistances(kk-1)) > tol
                        periapsisFound = 1;
                        firstPeriapsis_mult = [firstPeriapsis_mult; X_BCR_n_mult(kk,1:3)];
                        break
                    end
                end
                if periapsisFound == 0 && secondaryDistances(1) > secondaryDistances(end)
                    firstPeriapsis_mult = [firstPeriapsis_mult; X_BCR_n_mult(end,1:3)];
                elseif periapsisFound == 0 && secondaryDistances(1) < secondaryDistances(end)
                    firstPeriapsis_mult = [firstPeriapsis_mult; X_BCR_n_mult(1,1:3)];
                end
            end
        end
    end
    
    % Clearing iterative variables
    clear lat lon endFPA endSpeed_n X_BCR_n_mult
    
    %%% Plotting Results of Multiple Runs
    if plotBCR_n_mult == 1
        % BCR Trajectories
        figure; hold all
        plot3(r_BCR_mult(:,1),r_BCR_mult(:,2),r_BCR_mult(:,3),'g','linewidth',1.0)
        plot3(L2_n(1),L2_n(2),L2_n(3),'k^','markersize',5,'linewidth',1.5)
        plot3(L1_n(1),L1_n(2),L1_n(3),'k^','markersize',5,'linewidth',1.5)
        plotBody3(secondary.R_n, [1-secondary.MR, 0, 0], secondary.color)
        PlotBoi3('X_n','Y_n','Z_n',14)
        axis equal
%         xlim([1-secondary.MR-5*secondary.R_n, 1-secondary.MR+5*secondary.R_n])
%         ylim([-5*secondary.R_n, 5*secondary.R_n])
        xlim([L1_n(1)-2*secondary.R_n, L2_n(1)+2*secondary.R_n])
        ylim([-12*secondary.R_n, 12*secondary.R_n])
        zlim([-5*secondary.R_n, 5*secondary.R_n])
        title(secondary.name)
        view(0,90)
        if firstPeriapsisAnalysis_mult == 1
            alpha(.6)
        end
    end
    
    if anyLandings == 1
        if fpaVSspeed_mult == 1
            % Terminal FPA vs Speed
            figure; hold all
            plot(endFPAs_mult*180/pi,endSpeeds_n_mult, 'bo','markersize',10,'linewidth',1.5)
            PlotBoi2('Flight Path Angle at Crash, °','Normalized Speed at Crash',12)
        end
        
        if latlonPlot_mult == 1
            % Terminal Lat/Lon
            figure; hold all
            plot(latlons_mult(:,2), latlons_mult(:,1), 'rx','markersize',5,'linewidth',1.5);
            PlotBoi2('Longitude, °', 'Latitude, °',14)
            xlim([-180 180])
            ylim([-90 90])
            grid on
        end
    end
    
    if firstPeriapsisAnalysis_mult == 1
        figure
        subplot(1,2,1); hold all
        plot3(firstPeriapsis_mult(:,1),firstPeriapsis_mult(:,2),firstPeriapsis_mult(:,3),'rx','markersize',5,'linewidth',1)
        PlotBoi3('X_n','Y_n','Z_n',14)
        axis equal
        xlim([1-secondary.MR-5*secondary.R_n, 1-secondary.MR+5*secondary.R_n])
        ylim([-5*secondary.R_n, 5*secondary.R_n])
        zlim([-5*secondary.R_n, 5*secondary.R_n])
        
        subplot(1,2,2); hold all
        plot3(firstPeriapsis_mult(:,1),firstPeriapsis_mult(:,2),firstPeriapsis_mult(:,3),'rx','markersize',5,'linewidth',1)
        plot3(r_BCR_mult(:,1),r_BCR_mult(:,2),r_BCR_mult(:,3),'g','linewidth',1.0)
        plotBody3(secondary.R_n, [1-secondary.MR, 0, 0], secondary.color)
        alpha(.3)
        PlotBoi3('X_n','Y_n','Z_n',14)
        axis equal
        xlim([1-secondary.MR-5*secondary.R_n, 1-secondary.MR+5*secondary.R_n])
        ylim([-5*secondary.R_n, 5*secondary.R_n])
        zlim([-5*secondary.R_n, 5*secondary.R_n])
    end
    
end


















