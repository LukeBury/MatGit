clear
clc
close all
addpath(genpath('/Users/CU_Google_Drive/lukebury/Documents/MATLAB/mbin'))
ticWhole = tic;
% ========================================================================
%%% Run-Switches
% ========================================================================
plot_postBurnTrajectories = 0;
plot_burnQuivers          = 0;

saveFigures               = 0;

% ========================================================================
%%% Setup
% ========================================================================
% -------------------------------------------------
% Choose system
% -------------------------------------------------
%%% General data on solar system bodies
bodies = getBodyData();
%%% Color options/schemes
colors = get_colors();

%%% 3B system
primary   = bodies.jupiter;
secondary = bodies.europa;

%%% Normalizing constants
rNorm = secondary.a;         % n <-> km
tNorm = 1/secondary.meanMot; % n <-> sec
vNorm = rNorm / tNorm;       % n <-> km/sec

%%% Setting new variables because parfor is picky
R2_n = secondary.R_n;
MR = secondary.MR;
secondaryName = secondary.name;

% -------------------------------------------------
% Selecting initial conditions
% -------------------------------------------------
%%% Collinear lagrange point of interest
Lpoint = 1;  % 1 or 2

%%% Acquire Collinear Lagrange points
L123 = EquilibriumPoints(secondary.MR,1:3); % [3x3] of L1, L2, L3 normalized BCR coordinates
Lpoint_x = L123(Lpoint,1);

%%% Choosing number of burn locations along reference trajectory (nDVs) and
%%% number of burns per location (nAngles)
nAngles = 3; % 361
nDVs    = 3; % 200

%%% How fast the SC would be traveling over the Lagrange point
dvLp_mps = 1; % Meters per second

%%% Chooseing delta-V magnitude and converting to normalized units
dVMag_mps = 50; % meters per second
dVMag_n = dVMag_mps/1000/vNorm;

%%% Choosing number of post-burn JC bins for results and colors for those
%%% bins
binCount_JC = 4;
binColors = [colors.sch.d4_1(2,:);colors.sch.d4_1(1,:);colors.sch.d4_1(4,:);colors.sch.d4_1(3,:)];

%%% Neck positions
nPositions = 5;
ys = linspace(-0.9*R2_n, -0.7*R2_n, nPositions);

%%% Contour plotting options
nContoursPerBin = 4;

%%% Selecting time vector
t_i = 0; % sec
t_f = pi;
dt = t_f/10000;

%%% Store initial conditions
nICs = 5;
ICs{nICs} = [];
for kk = 1:nICs
    ICs{kk}.primary = primary;
    ICs{kk}.secondary = secondary;
    ICs{kk}.Lpoint = Lpoint;
    ICs{kk}.t_f = t_f;
    ICs{kk}.y0 = ys(kk);
    ICs{kk}.dvLp_mps = dvLp_mps;
    ICs{kk}.dVMag_mps   = dVMag_mps;
end


% % -------------------------------------------------
% % Finding upper-y-value of neck at L-point
% % -------------------------------------------------
% %%% Set values
% c = JC_scInitial;
% u = secondary.MR;
% x = L123(Lpoint,1);
% z = 0;
% 
% %%% JC function equal to zero
% f = @(y) x^2 + y.^2 - c + 2*(1-u)./sqrt((x+u)^2+y.^2+z^2) + 2*u./sqrt((x-1+u)^2+y.^2+z^2);
% 
% %%% Find the root of this function in the appropriate range
% y_neck = fzero(f,[0 6*secondary.R_n]);
% 
% %%% Clear variables
% clear c u x f z
% 
% 
% % -------------------------------------------------
% % Finding upper-z-value of neck at L-point
% % -------------------------------------------------
% %%% Set values
% c = JC_scInitial;
% u = secondary.MR;
% x = L123(Lpoint,1);
% y = 0;
% 
% %%% JC function equal to zero
% f = @(z) x^2 + y.^2 - c + 2*(1-u)./sqrt((x+u)^2+y^2+z.^2) + 2*u./sqrt((x-1+u)^2+y^2+z.^2);
% 
% %%% Find the root of this function in the appropriate range
% z_neck = fzero(f,[0 6*secondary.R_n]);
% 
% %%% Clear variables
% clear c u x f y

% -------------------------------------------------
% Define enumerated cell-column indices
% -------------------------------------------------
ci_X_impact        = 1;
ci_t_burn          = 2;
ci_burnQuivers     = 3;
ci_burnHat         = 4;
ci_postBurnJC      = 5;
ci_JCbin           = 6;
ci_t_impact        = 7;
ci_impactLatLons   = 8;
ci_impactSpeeds    = 9;
ci_impactAngles    = 10;
ci_X_BCRn_postBurn = 11;

% ========================================================================
%%% Loop through conditions, store results in bins, and plot
% ========================================================================
parfor ii = 1:length(ICs)
ticLoop = tic;
    % -------------------------------------------------
    % Unpack initial conditions
    % -------------------------------------------------
    primary = ICs{ii}.primary;
    secondary = ICs{ii}.secondary;
    Lpoint = ICs{ii}.Lpoint;
    t_f = ICs{ii}.t_f;
    y0 = ICs{ii}.y0;
    dvLp_mps = ICs{ii}.dvLp_mps;
    dVMag_mps = ICs{ii}.dVMag_mps;
    
    dVMag_n = dVMag_mps/1000/vNorm;
    
    % -------------------------------------------------
    % Reducing broadcast variables
    % -------------------------------------------------
    %%% Color options/schemes
    colors = get_colors();
    colorMatrix = binColors;

    %%% Reacquire Collinear Lagrange points
    L123 = EquilibriumPoints(MR,1:3); % [3x3] of L1, L2, L3 normalized BCR coordinates
    Lpoint_x = L123(Lpoint,1);
    
    % -------------------------------------------------
    % Setting Jacobi value of spacecraft relative to L-point
    % -------------------------------------------------
    %%% Converting this velocity to a jacobi constant value to be
    %%% differenced from the value of the chosen L-point
    dJC_vel_kps = dvLp_mps/1000;
    dJC_Lp = (dJC_vel_kps/vNorm)^2;

    %%% Jacobi constant of Lagrange point
    [JC_Lp] = JacobiConstantCalculator(secondary.MR,L123(Lpoint,:),[0,0,0]);

    %%% s/c starting JC (JC_scDesired) is lower than JC_Lp because the dJC_Lp is accounted for
    %%% in velocity
    JC_scInitial = JC_Lp-dJC_Lp;
    
    % -------------------------------------------------
    % Integrate Integrate Integrate Integrate Woah
    % -------------------------------------------------
    %%% Setting initial state
    r0_n = [Lpoint_x, ys(ii), 0];
    
    %%% With JC0 defined, starting velocity is a function of position. So first
    %%% we must calculate the JC of the stationary starting position
    JC_initialPos = JacobiConstantCalculator(MR,r0_n,[0,0,0]);

    %%% Starting velocity is found from difference between s/c JC (JC_scDesired) and the
    %%% JC of the stationary starting position (JC_initialPos)
    dJC_forInitialVelocity = JC_initialPos - JC_scInitial;

    %%% Setting variables so they're not considered 'temporary variables'
    v0Hat = [];
    v0_n  = [];
    
    %%% Find necessary velocity magnitude
    if dJC_forInitialVelocity < 0
        warning('Spacecraft starting in a forbidden zone')
    elseif dJC_forInitialVelocity >= 0
        v0i = sqrt(abs(dJC_forInitialVelocity));

        % Find direction for velocity
        if Lpoint == 1
            v0Hat = R3([1, 0, 0],0);
        elseif Lpoint == 2
            v0Hat = R3([-1, 0, 0],0);
        end

        % Create initial velocity vector
        v0_n = v0Hat .* v0i;
    end

    X0_n = [r0_n, v0_n];
    
    %%% Setting time vector
    time0_n = t_i:dt:t_f;

    %%% Choosing ode45 tolerance
    tol = 1e-10;

    %%% Setting integrator options
    options = odeset('RelTol',tol,'AbsTol',tol);
    options_Impact = odeset('Events',@event_Impact_CR3Bn,'RelTol',tol,'AbsTol',tol);

    %%% Propagating reference trajectory
    [time_n_ref, X_BCR_n_ref] = ode45(@Int_CR3Bn, time0_n, X0_n, options_Impact, MR, R2_n);
    
    % -------------------------------------------------
    % Calculate/Store reference data
    % -------------------------------------------------
    %%% Structure for storing reference info
    refTraj = struct();
    
    %%% Storing data we already have
    refTraj.X_BCR_n = X_BCR_n_ref;
    refTraj.t = time_n_ref;
    refTraj.JC = JC_scInitial;

    %%% Creating SCR position
    refEnd_SCR = refTraj.X_BCR_n(end,1:3) - [1-MR, 0, 0];

    %%% Finding lat/lon of reference end site
    [ref_lat, ref_lon] = ECEF2latlon(refEnd_SCR,'degrees','stupidMoon');

    %%% Storing lat/lon of reference end site
    refTraj.endLatLons = [ref_lat, ref_lon];
    
    %%% End speed of reference
    refTraj.endSpeed = norm(refTraj.X_BCR_n(end,4:6));
    refTraj.endSpeed_mps = refTraj.endSpeed * (vNorm*1000);
    
    %%% End "impact" angle of reference
    vHatRefEnd = refTraj.X_BCR_n(end,4:6)./norm(refTraj.X_BCR_n(end,4:6));
    % Angle between velocity vector and "surface"
    A = R3(refEnd_SCR,pi/2);
    B = vHatRefEnd;
    refTraj.impactAngle = acos(dot(A,B)/(norm(A)*norm(B)));
    if refTraj.impactAngle > pi/2
        refTraj.impactAngle = pi - refTraj.impactAngle;
    end
    
    % -------------------------------------------------
    % Find indices of reference trajectory to burn at
    % -------------------------------------------------
    %%% Calculate distance traveled at each state
    distTraveled = 0;
    sum = 0;
    for kk = 2:size(X_BCR_n_ref,1)
        distTraveled(kk) = distTraveled(kk-1) + norm(X_BCR_n_ref(kk,1:3)-X_BCR_n_ref(kk-1,1:3));
    end

    %%% Set target distances to burn at
    burnDistances = linspace(0,distTraveled(end-1),nDVs);

    %%% For each of these distances, find closest index of reference
    %%% trajectory (since they won't match exactly)
    burnIndices = [];
    for kk = 1:length(burnDistances)
        ti = find(abs(distTraveled - burnDistances(kk)) == min(abs(distTraveled - burnDistances(kk))));
        burnIndices = [burnIndices, ti];
    end

    % -------------------------------------------------
    % Set various heading angles for burns
    % -------------------------------------------------
    headingAngles = 0:(2*pi/nAngles):(2*pi)-(2*pi/nAngles);
    
    traj_i = 0;
    impactTraj = {};
    for bi = burnIndices
        for dV_headingAngle = headingAngles
            %%% Advance the trajecotry index
            traj_i = traj_i + 1;
            
            %%% Assigning time of burn
            t_dV = time_n_ref(bi);

            %%% Assigning pre-dV state
            X_predV_n = X_BCR_n_ref(bi,:);

            %%% Adding tangent dV to states
            % Finding vHat
            vHat = X_predV_n(4:6)./norm(X_predV_n(4:6));
            vHat = R3(vHat,dV_headingAngle);
            % Determining dV
            dV_n = vHat .* dVMag_n;
            % Creating post-dV state
            X0_postdV_n = [X_predV_n(1:3), X_predV_n(4:6)+dV_n];
            
            %%% Integrate
            [time_n_postBurn, X_BCR_n_postBurn, time_eventImpact, X_eventImpact, index_eventImpact] = ode45(@Int_CR3Bn,...
                [t_dV:dt:t_f], X0_postdV_n, options_Impact, MR, R2_n);
            
            %%% Post-burn JC
            JC_postBurn = JacobiConstantCalculator(MR, [X_BCR_n_postBurn(1,1), X_BCR_n_postBurn(1,2),...
                X_BCR_n_postBurn(1,3)], [X_BCR_n_postBurn(1,4), X_BCR_n_postBurn(1,5), X_BCR_n_postBurn(1,6)]);
            
            % -------------------------------------------------
            % Calculating Lat/Lon impact locations
            % -------------------------------------------------
            impactLatLons = [];
            for kk = 1:size(X_eventImpact,1)
                %%% Creating SCR position
                rImpact_SCR = X_eventImpact(kk,1:3) - [1-MR, 0, 0];

                %%% Finding lat/lon of impact site
                [lat, lon] = ECEF2latlon(rImpact_SCR,'degrees','stupidMoon');

                %%% Storing lat/lon of impact site
                impactLatLons = [impactLatLons; lat, lon];
            end
            
            % -------------------------------------------------
            % Calculating impact velocities and angles
            % -------------------------------------------------
            impactSpeeds = zeros(size(X_eventImpact,1));
            impactAngles = zeros(size(X_eventImpact,1));

            for kk = 1:size(X_eventImpact,1)
                %%% Norm of impact velocity
                impactSpeeds(kk) = norm(X_eventImpact(kk,4:6));

                %%% Creating SCR position vector
                rImpact_SCR_n = X_eventImpact(kk,1:3) - [1-MR,0,0];
                rImpact_SCR_n = rImpact_SCR_n./norm(rImpact_SCR_n);

                %%% Velocity unit vector at impact
                vHatImpact_n = X_eventImpact(kk,4:6)./norm(X_eventImpact(kk,4:6));

                %%% Angle between velocity and surface
                A = R3(rImpact_SCR_n,pi/2);
                B = vHatImpact_n;
                impactAngles(kk) = acos(dot(A,B)/(norm(A)*norm(B)));
                if impactAngles(kk) > pi/2
                    impactAngles(kk) = pi - impactAngles(kk);
                end
            end

            % -------------------------------------------------
            % Storing information from this trajectory
            % -------------------------------------------------
            impactTraj{traj_i}(ci_X_impact)        = {X_eventImpact};
            impactTraj{traj_i}(ci_t_burn)          = {t_dV};
            impactTraj{traj_i}(ci_burnQuivers)     = {[X0_postdV_n(1:3), dV_n]};
            impactTraj{traj_i}(ci_burnHat)         = {vHat};
            impactTraj{traj_i}(ci_postBurnJC)      = {JC_postBurn};
            impactTraj{traj_i}(ci_t_impact)        = {time_eventImpact};
            impactTraj{traj_i}(ci_impactLatLons)   = {impactLatLons};
            impactTraj{traj_i}(ci_impactSpeeds)    = {impactSpeeds};
            impactTraj{traj_i}(ci_impactAngles)    = {impactAngles};
            if plot_postBurnTrajectories == 1
                impactTraj{traj_i}(ci_X_BCRn_postBurn) = {X_BCR_n_postBurn};
            end
        end
    end
    
    % -------------------------------------------------
    % Creating JC bins for this set of ICs and store results in bins
    % -------------------------------------------------
    %%% Organize impactTraj into single cell array
    impactTrajs = vertcat(impactTraj{:});
    
    %%% Determine minimum and maximum JC after burn
    JC_min = min([impactTrajs{:,ci_postBurnJC}]);
    JC_max = max([impactTrajs{:,ci_postBurnJC}]);
    
    %%% Create cutoff values for bins
    bins_JC = linspace(JC_min, JC_max, binCount_JC+1);
    bins_JC(end) = bins_JC(end)+1; % make last bin big so everything fits in first set
    
    %%% Preallocate structure array
    trajBin = struct();
    trajBin(binCount_JC).X_impact            = [];
    trajBin(binCount_JC).impactLatLons       = [];
    trajBin(binCount_JC).impactSpeeds        = [];
    trajBin(binCount_JC).impactSpeeds_mps    = [];
    trajBin(binCount_JC).t_impact            = [];
    trajBin(binCount_JC).impactAngles        = [];
    trajBin(binCount_JC).postBurnJC          = [];
    trajBin(binCount_JC).t_burn              = [];
    trajBin(binCount_JC).burnQuivers         = [];
    trajBin(binCount_JC).burnHat             = [];
    if plot_postBurnTrajectories == 1
        trajBin(binCount_JC).X_BCRn_postBurn = [];
    end
    
    %%% Assign bins
    for kk = 1:size(impactTrajs,1)
        if isempty(impactTrajs{kk,ci_X_impact}) == 0
            %%% Assign current bin number
            impactTrajs{kk,ci_JCbin} = discretize(impactTrajs{kk,ci_postBurnJC},bins_JC);
            bi = impactTrajs{kk,ci_JCbin};

            %%% Store data in structure
            trajBin(bi).X_impact             = [trajBin(bi).X_impact;         impactTrajs{kk,ci_X_impact}];
            trajBin(bi).impactLatLons        = [trajBin(bi).impactLatLons;    impactTrajs{kk,ci_impactLatLons}];
            trajBin(bi).impactSpeeds         = [trajBin(bi).impactSpeeds;     impactTrajs{kk,ci_impactSpeeds}];
            trajBin(bi).t_impact             = [trajBin(bi).t_impact;         impactTrajs{kk,ci_t_impact}];
            trajBin(bi).impactAngles         = [trajBin(bi).impactAngles;     impactTrajs{kk,ci_impactAngles}];
            trajBin(bi).postBurnJC           = [trajBin(bi).postBurnJC;       impactTrajs{kk,ci_postBurnJC}];
            trajBin(bi).t_burn               = [trajBin(bi).t_burn;           impactTrajs{kk,ci_t_burn}];
            trajBin(bi).burnQuivers          = [trajBin(bi).burnQuivers;      impactTrajs{kk,ci_burnQuivers}];
            trajBin(bi).burnHat              = [trajBin(bi).burnHat;          impactTrajs{kk,ci_burnHat}];
            if plot_postBurnTrajectories == 1
                trajBin(bi).X_BCRn_postBurn  = [trajBin(bi).X_BCRn_postBurn; padcat(impactTrajs{:,ci_})];
            end
        end
    end
    
    %%% Assign colors to bins
    for bi = 1:binCount_JC
        trajBin(bi).bin_col = colorMatrix(bi,:);
    end
    
    %%% Addings fields denoting certain units
    for bi = 1:binCount_JC
        trajBin(bi).impactSpeeds_mps = trajBin(bi).impactSpeeds .* (vNorm*1000);
    end
    
    % =================================================
    % Pre-bin plotting
    % =================================================
    % -------------------------------------------------
    % Set up plotting options
    % -------------------------------------------------
    markerSize = 13;
    
    % ------------------
    % figure(1), subplot(2,4,[1 2]) - system, apses, and impacts
    % ------------------
    figure(1); 
    % set(gcf,'Position',[266 54 986 678]); % big
    set(gcf,'Position',[1 1 1440 804]); % fullscreen mac
    subplot(2,4,[1 2]); hold all
    axis equal
    PlotBoi3('x$_n$','y$_n$','z$_n$',16,'Latex') 

    %%% Title
    title_body = secondaryName;
    title_J0 = ['J$_0$ = J',sprintf('$_{L%1.0d}$ + %2.2e',Lpoint,-dJC_Lp)];
    title_v0 = ['$\mathopen|\Delta V_',sprintf('{L%1.0d}',Lpoint),' \mathclose|$ = ',sprintf('%2.1f',sqrt(dJC_Lp)*vNorm*1000), ' $m/s$'];
    title_R_n = sprintf('$R_n$ = %1.2e',R2_n);
    title_dV = ['$\mathopen|\Delta V \mathclose|$ = ',sprintf('%2.0f',dVMag_mps), ' $m/s$'];

    title({[title_body, ', \hspace{0.5cm}', title_J0, ', \hspace{0.5cm}', title_v0,...
             ', \hspace{0.5cm}', title_R_n, ', \hspace{0.5cm}', title_dV], ''},'Interpreter','Latex');

    %%% Creating ZV contours
    % Contour bounds
    xCont_min = L123(1,1)-2*R2_n;
    xCont_max = L123(2,1)+2*R2_n;
    yCont_min = -6*R2_n;
    yCont_max = 6*R2_n;

    % Creating x-y grid
    x_grid = linspace(xCont_min,xCont_max,600);
    y_grid = linspace(yCont_min,yCont_max,200);
    [X, Y] = meshgrid(x_grid,y_grid);
    
    % Calculating JCs across grid
    Z = zeros(size(X));
    for xk = 1:size(X,1)
        for yk = 1:size(X,2)
            %%% Zero-Velocity Curve
            zv = JacobiConstantCalculator(MR,[X(xk,yk), Y(xk,yk), 0] ,[0, 0, 0]);
            Z(xk,yk) = zv;
        end
    end

    %%% Plot limit
    xlim([xCont_min xCont_max])
    ylim([yCont_min yCont_max])

    %%% Plotting secondary
    plotBody2(R2_n,[1-MR,0,0],[1,1,1], colors.std.black,1.5,0)

    % ------------------
    % Colorbar for whole figure
    % ------------------
    h8 = get(subplot(2,4,8),'Position');
    cbar = colorbar('Position', [h8(1)+h8(3)+0.01  h8(2)  0.025  h8(2)+h8(3)*4.6]);
    colormap(colorMatrix)
    caxis([0 size(colorMatrix,1)]);
    cbar.FontName     = 'Arial';
    cbar.FontSize     = 10;
    cbar.Ticks        = sort(0:binCount_JC);
    TickLabels = num2cell(bins_JC(1:end-1)-JC_Lp);
    for kk = 1:length(TickLabels)
        TickLabels{kk} = sprintf('%1.2e',TickLabels{kk});
    end
    cbar.TickLabels   = [TickLabels, sprintf('J_{L%1.0d} +',Lpoint)];
    cbar.Label.String = {'Post-Burn Jacobi Constant'};
    cbar.Label.Rotation = 0; % 0 = horizontal, 90 = vertical
    cbar.Label.Position = [h8(4)+.2, 4.15, 0];

    % -------------------------------------------------
    % Loop through bins and plot data
    % -------------------------------------------------
    for bi = 1:binCount_JC
        %%% Check if this bin has any impact
        if isempty(trajBin(bi).X_impact) == 1
            continue
        end
        % -------------------------------------------------
        % figure(1), subplot(2,4,[1 2]) - system, apses, and impacts
        % -------------------------------------------------
        figure(1); subplot(2,4,[1 2]); hold all

        %%% Plotting post-burn contours
        orderedPostBurnJCs = sort(trajBin(bi).postBurnJC);
        contoursForPlotting = linspace(orderedPostBurnJCs(1),orderedPostBurnJCs(end),nContoursPerBin);
        for currentContour = contoursForPlotting
            [C2,h2] = contour(X,Y,Z,[currentContour, currentContour],'color',trajBin(bi).bin_col,'linewidth',3);
        end

        % %%% Plotting post-burn FULL trajectory
        if plot_postBurnTrajectories == 1
            plot3(trajBin(bi).X_BCR_n_postBurn(:,1),trajBin(bi).X_BCR_n_postBurn(:,2),trajBin(bi).X_BCR_n_postBurn(:,3),...
                'linewidth',1,'color',trajBin(kk).bin_col)
        end

        %%% Plotting impact locations
        if isempty(trajBin(bi).X_impact) == 0
            plot3(trajBin(bi).X_impact(:,1),trajBin(bi).X_impact(:,2),trajBin(bi).X_impact(:,3),'.','markersize',markerSize,...
                            'color',trajBin(bi).bin_col,'linewidth',2)
        end

        %%% Plotting burn quivers
        if plot_burnQuivers == 1
            qS = 0.002; % quiverScale
            quiver3(trajBin(bi).burnQuivers(:,1),trajBin(bi).burnQuivers(:,2),trajBin(bi).burnQuivers(:,3),trajBin(bi).burnHat(:,1).*qS,...
                    trajBin(bi).burnHat(:,2).*qS,trajBin(bi).burnHat(:,3).*qS,0,'color',trajBin(bi).bin_col,'linewidth',0.75,'maxheadsize',0.4);
        end

        % -------------------------------------------------
        % figure(1), subplot(2,4,3) - impact longitude vs post burn Jacobi
        % value
        % -------------------------------------------------
        subplot(2,4,3); hold all
        xmin3 = -180;
        xmax3 = 180;
        xlim([-180 180])
        ymin3 = min(vertcat(trajBin(:).postBurnJC))*.9999;
        ymax3 = max(vertcat(trajBin(:).postBurnJC))*1.0001;
        ylim([ymin3 ymax3])
        PlotBoi2('Longitude, $^\circ$','Post-Burn Jacobi Value',12,'Latex')
    
        %%% If there was an impact
        if isempty(trajBin(bi).X_impact) == 0
            %%% Plotting impact longitude vs post burn jacobi value
            plot(trajBin(bi).impactLatLons(:,2),trajBin(bi).postBurnJC,'.','markersize',markerSize,...
                            'color',trajBin(bi).bin_col,'linewidth',2)
                        
            %%% Plotting anti-primary line (0-deg longitude)
            [antiPrimaryLine] = plot([0 0],[0 max(vertcat(trajBin(:).postBurnJC))].*(1.1),'color',colors.std.blue,'linewidth',1.5);
            
            %%% Plotting reference info
            [refPlot] = plot(refTraj.endLatLons(2), refTraj.JC, 'kx', 'markersize', 12, 'linewidth', 1.5);
            
            %%% Legend
            legend([antiPrimaryLine, refPlot],{'Anti-Primary', 'Reference'},'Interpreter','latex','location','best')
        end

        % -------------------------------------------------
        % figure(1), subplot(2,4,4) - impact longitude vs impact angle
        % -------------------------------------------------
        subplot(2,4,4); hold all
        xlim([-180 180])
        ylim([0 90])
        PlotBoi2('Longitude, $^\circ$','Impact Angle, $^\circ$',12,'Latex')

        %%% If there was an impact
        if isempty(trajBin(bi).X_impact) == 0
            %%% Plotting impact longitude vs post burn jacobi value
            plot(trajBin(bi).impactLatLons(:,2),trajBin(bi).impactAngles.*(180/pi),'.','markersize',markerSize,...
                            'color',trajBin(bi).bin_col,'linewidth',2)
            %%% Plotting anti-primary line (0-deg longitude)
            [antiPrimaryLine] = plot([0 0],[0 90],'color',colors.std.blue,'linewidth',1.5);
            
            %%% Plotting reference info
            [refPlot] = plot(refTraj.endLatLons(2), refTraj.impactAngle, 'kx', 'markersize', 12, 'linewidth', 1.5);
            
%             %%% legend
%             legend([antiPrimaryLine],{'Anti-Primary'},'Interpreter','latex','location','best')

        end

        % -------------------------------------------------
        % figure(1), subplot(2,4,5) - impact longitude vs impact speed
        % -------------------------------------------------
        subplot(2,4,5); hold all
        xlim([-180 180])
        ylim([min(vertcat(trajBin(:).impactSpeeds_mps))*.999 max(vertcat(trajBin(:).impactSpeeds_mps))*1.001])
        PlotBoi2('Longitude, $^\circ$','Impact Speed, $m/s$',12,'Latex')

        %%% If there was an impact
        if isempty(trajBin(bi).X_impact) == 0
            %%% Plotting impact longitude vs speed
            plot(trajBin(bi).impactLatLons(:,2),trajBin(bi).impactSpeeds_mps,'.','markersize',markerSize,...
                            'color',trajBin(bi).bin_col,'linewidth',2)

            %%% Plotting anti-primary line (0-deg longitude)
            [antiPrimaryLine] = plot([0 0],[0 max(vertcat(trajBin(:).impactSpeeds_mps))].*1.1,'color',colors.std.blue,'linewidth',1.5);
            
            %%% Plotting reference info
            [refPlot] = plot(refTraj.endLatLons(2), refTraj.endSpeed_mps, 'kx', 'markersize', 12, 'linewidth', 1.5);
            
%             %%% legend
%             legend([antiPrimaryLine],{'Anti-Primary'},'Interpreter','latex','location','best')
        end

        % -------------------------------------------------
        % figure(1), subplot(2,4,6) - impact longitude vs impact time
        % -------------------------------------------------
        subplot(2,4,6); hold all
        xlim([-180 180])
        PlotBoi2('Longitude, $^\circ$','Time to Impact, $\pi$',12,'Latex')

        %%% Plotting impact angle vs impact time
        if isempty(trajBin(bi).X_impact) == 0
            %%% Plotting impact longitude vs speed
            plot(trajBin(bi).impactLatLons(:,2),trajBin(bi).t_impact./pi,'.','markersize',markerSize,...
                            'color',trajBin(bi).bin_col,'linewidth',2)

%             %%% Plotting time of reference impact
%             [referenceImpactLine] = plot([-180 180],[time_n_ref(end) time_n_ref(end)]./pi,'color',colors.std.black,'linewidth',1.5);

            %%% Plotting anti-primary line (0-deg longitude)
            [antiPrimaryLine] = plot([0 0],[0 max(vertcat(trajBin(:).t_impact))].*(1.1/pi),'color',colors.std.blue,'linewidth',1.5);
            
            %%% Plotting reference info
            [refPlot] = plot(refTraj.endLatLons(2), refTraj.t(end)/pi, 'kx', 'markersize', 12, 'linewidth', 1.5);
            
%             %%% legend
%             legend([antiPrimaryLine, referenceImpactLine],{'Anti-Primary','Reference End'},'Interpreter','latex','location','best')
        end

        % -------------------------------------------------
        % figure(1), subplot(2,4,7) - impact longitude vs time of dV
        % -------------------------------------------------
        subplot(2,4,7); hold all
        xlim([-180 180])
        ylim([0 max(vertcat(trajBin(:).t_burn))*1.1/pi])
        PlotBoi2('Longitude, $^\circ$','Time of $\Delta V$, $\pi$',12,'Latex')

        %%% Plotting impact angle vs impact time
        if isempty(trajBin(bi).X_impact) == 0
            %%% Plotting impact longitude vs speed
            plot(trajBin(bi).impactLatLons(:,2),trajBin(bi).t_burn./pi,'.','markersize',markerSize,...
                            'color',trajBin(bi).bin_col,'linewidth',2)

            %%% Plotting anti-primary line (0-deg longitude)
            [antiPrimaryLine] = plot([0 0],[0 max(vertcat(trajBin(:).t_impact))].*(1.1),'color',colors.std.blue,'linewidth',1.5);

%             %%% legend
%             legend([antiPrimaryLine],{'Anti-Primary'},'Interpreter','latex','location','best')
        end

        % -------------------------------------------------
        % figure(1), subplot(2,4,8) - impact angle vs impact speed
        % -------------------------------------------------
        subplot(2,4,8); hold all
        xlim([0 90])
        ylim([min(vertcat(trajBin(:).impactSpeeds_mps))*.999 max(vertcat(trajBin(:).impactSpeeds_mps))*1.001])
        PlotBoi2('Impact Angle, $^\circ$', 'Impact Speed, $m/s$',12,'Latex')

        %%% Plotting impact angle vs speed
        if isempty(trajBin(bi).X_impact) == 0
            plot(trajBin(bi).impactAngles.*(180/pi),trajBin(bi).impactSpeeds_mps,'.','markersize',markerSize,...
                            'color',trajBin(bi).bin_col,'linewidth',2)
                        
            %%% Plotting reference info
            [refPlot] = plot(refTraj.impactAngle, refTraj.endSpeed_mps, 'kx', 'markersize', 12, 'linewidth', 1.5);
        end

    end
    
    % -------------------------------------------------
    % Plot things to be on top
    % -------------------------------------------------
    % ------------------
    % figure(1), subplot(2,4,[1 2 3 4]) - system, apses, and impacts
    % ------------------
    figure(1); subplot(2,4,[1 2]); hold all
    %%% Plotting pre-burn trajectory
    plot3(X_BCR_n_ref(:,1),X_BCR_n_ref(:,2),X_BCR_n_ref(:,3),'linewidth',2,'color',colors.std.black)

    %%% Plotting Lagrange Points
    plot3(L123(1:2,1), L123(1:2,2), L123(1:2,3), '^', 'markeredgecolor', colors.std.black,...
                                           'markerfacecolor',colors.std.mag,'markersize',10);

    %%% Plotting pre-burn contour
    [C1,h1] = contour(X,Y,Z,[JC_scInitial, JC_scInitial],'color',colors.std.black,'linewidth',2);

    % -------------------------------------------------
    % Save figure
    % -------------------------------------------------
    if saveFigures == 1
        figure(1)
        figPath = '/Users/lukebury/Documents/MATLAB/CU/Research/Moon_Landing/Figures/';
        figName = sprintf('%s%sL%1d_%1.1fpi_y%1.2fR_dvLp%2.2f_dv%2.1f.png', figPath, secondaryName(1:3), Lpoint, t_f/pi, r0_n(2)/R2_n, dvLp_mps, dVMag_mps);
        % ex: EurL1_1.0pi_y-0.90R_dvLp56.50_dv50.0.png
        saveas(gcf,figName)
        close all
    end

ii
toc(ticLoop)
end


toc(ticWhole)

















