clear
clc
close all
addpath(genpath('/Users/CU_Google_Drive/lukebury/Documents/MATLAB/mbin'))
ticWhole = tic;
% ========================================================================
%%% Run Switches
% ========================================================================
%%% Data switches
store_trajectories     = 0;

%%% Plot switches
plot_initialConditions = 1;
plot_fullSystemContour = 1;
plot_sectionColors    = 1;

%%% Save plots?d
saveFigures            = 0;

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
% Grid-search settings
% -------------------------------------------------
%%% Collinear lagrange point of interest
Lpoint = 2;  % 1 or 2

%%% Acquire Collinear Lagrange points
L123 = EquilibriumPoints(secondary.MR,1:3); % [3x3] of L1, L2, L3 normalized BCR coordinates
Lpoint_x = L123(Lpoint,1);

% %%% Choosing number of elevations and azimuths per starting spot in 3D neck
% n_el = 3;
% n_az = 3;

%%% How fast the SC would be traveling over the Lagrange point
% dvLp_mps = 200; % Meters per second - Europa
% dvLp_mps = 50; % Meters per second - Enceladus
dvLp_mps = 200; % Meters per second

%%% Spacing of initial positions within 3D neck
% r0GridSpacing_km = 100; % km - Europa
% r0GridSpacing_km = 10; % km - Enceladus
r0GridSpacing_km = 1500; % km

%%% Spacing between azimuths and elevations of v0s per r0
% v0AngularSpace_deg = 15; % degrees - Europa
% v0AngularSpace_deg = 15; % degrees - Enceladus
v0AngularSpace_deg = 90; % degrees

%%% Selecting time vector
t_i = 0; % sec
t_f = 4*pi;
dt = t_f/10000;

% -------------------------------------------------
% Trouble with settings?
% -------------------------------------------------
if v0AngularSpace_deg > 90
    warning('Angular spacing too large for v0s')
    return
end

% -------------------------------------------------
% Bin settings
% -------------------------------------------------
%%% Choosing number of impact angle bins and colors for those bins
binCount_ImpactAngles  = 5;
binColors_ImpactAngles = [colors.sch.s5_1(5,:);colors.sch.s5_1(4,:);...
                         colors.sch.s5_1(3,:);colors.sch.s5_1(2,:);...
                         colors.sch.s5_1(1,:)];
                     
%%% Choosing color for bins of neck sections 
binCount_NeckSections  = 4;
binColors_NeckSections = [colors.sch.d4_1(2,:);colors.sch.d4_1(1,:);...
                         colors.sch.d4_1(4,:);colors.sch.d4_1(3,:)];

% -------------------------------------------------
% Finding initial JC of spacecraft
% -------------------------------------------------
%%% Jacobi constant of Lagrange point
[JC_Lp] = JacobiConstantCalculator(secondary.MR,L123(Lpoint,:),[0,0,0]);

dJC_vel_kps = dvLp_mps/1000;
dJC_Lp = (dJC_vel_kps/vNorm)^2;
    
%%% s/c starting JC (JC_scDesired) is lower than JC_Lp because the dJC_Lp is accounted for
%%% in velocity
JC_scInitial = JC_Lp-dJC_Lp;

% -------------------------------------------------
% Finding upper-y-value of neck at L-point
% -------------------------------------------------
%%% Set values
c = JC_scInitial;
u = secondary.MR;
x = L123(Lpoint,1);
z = 0;

%%% JC function equal to zero
f = @(y) x^2 + y.^2 - c + 2*(1-u)./sqrt((x+u)^2+y.^2+z^2) + 2*u./sqrt((x-1+u)^2+y.^2+z^2);

%%% Find the root of this function in the appropriate range
y_neck_upper = fzero(f,[0 6*secondary.R_n]);

%%% Clear variables
clear c u x f z


% -------------------------------------------------
% Finding upper-z-value of neck at L-point
% -------------------------------------------------
%%% Set values
c = JC_scInitial;
u = secondary.MR;
x = L123(Lpoint,1);
y = 0;

%%% JC function equal to zero
f = @(z) x^2 + y.^2 - c + 2*(1-u)./sqrt((x+u)^2+y^2+z.^2) + 2*u./sqrt((x-1+u)^2+y^2+z.^2);

%%% Find the root of this function in the appropriate range
z_neck_upper = fzero(f,[0 6*secondary.R_n]);

%%% Clear variables
clear c u x f y

% -------------------------------------------------
% Create grid of starting locations based on y-z neck
% -------------------------------------------------
r0GridSpacing_n = r0GridSpacing_km/rNorm;

ys = -y_neck_upper*2 : r0GridSpacing_n : y_neck_upper*2;
zs = -z_neck_upper*2 : r0GridSpacing_n : z_neck_upper*2;
[Y_yz,Z_yz] = meshgrid(ys,zs);

% -------------------------------------------------
% Only keep starting positions that are valid for the energy level
% -------------------------------------------------
%%% Calculating JCs across y-z grid
JCs_yz_Lpoint = zeros(size(Y_yz));
for yk = 1:size(Y_yz,1)
    for zk = 1:size(Y_yz,2)
        %%% Zero-Velocity Curve
        zv = JacobiConstantCalculator(secondary.MR,[L123(Lpoint,1), Y_yz(yk,zk), Z_yz(yk,zk)] ,[0, 0, 0]);
        JCs_yz_Lpoint(yk,zk) = zv;
    end
end

%%% Get points of y-z contour in 3D space
[ yzContourPoints ] = getContourPoints( Y_yz, Z_yz, JCs_yz_Lpoint, JC_scInitial );

%%% Search grid for points lying within y-z contour and keep those as r0s
r0s = [];
for yk = 1:size(Y_yz,1)
    for zk = 1:size(Z_yz,2)
        if inpolygon(Y_yz(yk,zk), Z_yz(yk,zk),yzContourPoints(1,:),yzContourPoints(2,:)) == 1
            r0s = [r0s; L123(Lpoint,1), Y_yz(yk,zk), Z_yz(yk,zk)];
        end
    end
end

n_r0s = size(r0s,1);

% -------------------------------------------------
% Creating Bins for impact angle and neck sections
% -------------------------------------------------
%%% Bins for impact angles
% bins_impactAngles = ceil([0, linspace(5, 90, binCount_ImpactAngles)]);
% % Making last bin too big so everything fits in bins(1:end-1)
% bins_impactAngles(end) = bins_impactAngles(end) + 1;
bins_impactAngles = [0, 5, 20, 40, 60, 91];

%%% Bins for neck sections
smallNeck = smallerValue(z_neck_upper, y_neck_upper);
largeNeck = largerValue(z_neck_upper, y_neck_upper);
bins_neckSections = [linspace(0,smallNeck,binCount_NeckSections), largeNeck];

% -------------------------------------------------
% Defining azimuths and elevations to loop through per r0
% -------------------------------------------------
azimuths   = 0 : v0AngularSpace_deg : 360-v0AngularSpace_deg; % deg
elevations = 0 : v0AngularSpace_deg : 90; % deg

n_v0s_per_r0 = (length(azimuths)*length(elevations)-(length(azimuths)-1));
% -------------------------------------------------
% Predefining some things before parfor
% -------------------------------------------------
%%% Setting time vector
time0_n = t_i:dt:t_f;

%%% Choosing ode45 tolerance
tol = 1e-10;

%%% Setting integrator options
options = odeset('RelTol',tol,'AbsTol',tol);
options_Impact = odeset('Events',@event_Impact_CR3Bn,'RelTol',tol,'AbsTol',tol);

% ========================================================================
%%% Loop through conditions, store results in bins, and plot
% ========================================================================
r0Data = {};
% parforSize = 2;
if Lpoint ~= 2
    warning('must reconfigure v0 calculator for this Lpoint')
end
parfor ii = 1:n_r0s
ticLoop = tic;

    % -------------------------------------------------
    % Reducing broadcast variables
    % -------------------------------------------------
    L123mat = L123;
    % -------------------------------------------------
    % Setting current initial position
    % -------------------------------------------------
    r0_i = r0s(ii,:);
    
    % -------------------------------------------------
    % Finding necessary v0 to maintain correct JC0
    % -------------------------------------------------
    %%% With JC0 defined, starting velocity is a function of position. So first
    %%% we must calculate the JC of the stationary starting position
    JC_initialPos = JacobiConstantCalculator(MR,r0_i,[0,0,0]);
    
    %%% Starting velocity is found from difference between s/c JC (JC_scDesired) and the
    %%% JC of the stationary starting position (JC_initialPos)
    dJC_forInitialVelocity = JC_initialPos - JC_scInitial;
    
    %%% Setting variables so they're not considered 'temporary variables'
    vRef  = [];
    
    %%% Find necessary velocity magnitude
    if dJC_forInitialVelocity < 0
        warning('Spacecraft starting in a forbidden zone')
    elseif dJC_forInitialVelocity >= 0
        v0i = sqrt(abs(dJC_forInitialVelocity));
        
        v0Hat = [0, 0, 1];
        
        % Create initial velocity vector
        vRef = v0Hat .* v0i;
    end
    
    % -------------------------------------------------
    % Looping through velocity vectors
    % -------------------------------------------------
    trajCount = 0;
    elevation90HasBeenDone = 0;
    
    %%% Initializing data matrices
    X0s_ii               = zeros(n_v0s_per_r0,6);
    latLons_ii           = zeros(n_v0s_per_r0,2);
    trajs_ii             = zeros(n_v0s_per_r0,3);
    impactSpeeds_ii      = zeros(n_v0s_per_r0,1);
    impactAngles_ii      = zeros(n_v0s_per_r0,1);
    
    bin_impactAngles_ii  = zeros(n_v0s_per_r0,1);
    bin_neckSections_ii  = zeros(n_v0s_per_r0,1);
    
    %%% Running loop
    for az_i = azimuths
        for el_i = elevations
            %%% Making sure elevation-90 isn't redone
            if el_i == 90 && elevation90HasBeenDone == 1
                continue
            end

            %%% Counting IC
            trajCount = trajCount + 1;

            %%% Determining current v0
            v1 = R2(vRef, -el_i*pi/180);
            v0_i = R1(v1, az_i*pi/180);
            
            %%% Setting X0
            X0_n = [r0_i, v0_i];
            

            % ---------------------------------------
            % Integrating
            % ---------------------------------------
            %%% Propagating trajectory
            [time_n, X_BCR_n, time_eventImpact, X_eventImpact, index_eventImpact] = ode45(@Int_CR3Bn,...
                time0_n, X0_n, options_Impact, MR, R2_n);
            
            % ---------------------------------------
            % Preallocating temporary variables
            % ---------------------------------------
            impactLatLon = []
            impactSpeed  = [];
            impactAngle  = [];
            
            % ---------------------------------------
            % Determining lat/lon, speed, and angle of impact
            % ---------------------------------------
            
            if isempty(X_eventImpact) == 0
                % --------------------------
                % Impact lat/lon
                % --------------------------
                %%% Creating SCR position
                rImpact_SCR = X_eventImpact(1,1:3) - [1-MR, 0, 0];

                %%% Finding lat/lon of impact site
                [lat, lon] = ECEF2latlon(rImpact_SCR,'degrees','stupidMoon');
                impactLatLon = [lat, lon];
                
                % --------------------------
                % Impact speed
                % --------------------------
                impactSpeed = norm(X_eventImpact(end,4:6));
                
                % --------------------------
                % Impact Angle
                % --------------------------
                %%% Creating SCR position vector
                rImpact_SCR_n = X_eventImpact(end,1:3) - [1-MR,0,0];
                rImpact_SCR_n = rImpact_SCR_n./norm(rImpact_SCR_n);
                
                %%% Velocity unit vector at impact
                vHatImpact_n = X_eventImpact(end,4:6)./norm(X_eventImpact(end,4:6));
                
                %%% Angle between velocity and surface
                A = R3(rImpact_SCR_n,pi/2);
                B = vHatImpact_n;
                impactAngle = acos(dot(A,B)/(norm(A)*norm(B)));
                if impactAngle > pi/2
                    impactAngle = pi - impactAngle;
                end
                
                %%% Going with degrees
                impactAngle = impactAngle*180/pi;
                
            elseif isempty(X_eventImpact) == 1
                impactLatLon = [NaN, NaN];
                impactSpeed  = NaN;
                impactAngle  = NaN;
                
            end
            
            % =======================================
            % Storing results of trajectory
            % ---------------------------------------
            X0s_ii(trajCount,:)        = X0_n;
            latLons_ii(trajCount,:)    = impactLatLon;
            impactSpeeds_ii(trajCount) = impactSpeed;
            impactAngles_ii(trajCount) = impactAngle;
            
            % ---------------------------------------
            % Assigning bins
            % ---------------------------------------
            %%% Bin - Impact Angle
            bin_impactAngles_ii(trajCount) = discretize(impactAngle, bins_impactAngles);
    
            %%% Bin - Neck Section
            % Calculate radial distance from r0 to LPoint
            r0LPointDist = norm(X0_n(1:3) - L123mat(Lpoint,:));
            
            bin_neckSections_ii(trajCount) = discretize(r0LPointDist, bins_neckSections);
            
%             if     X0_n(2) >  0 && X0_n(3) >= 0
%                 bin_neckSections_ii(trajCount) = 1;
%             elseif X0_n(2) >= 0 && X0_n(3) <  0
%                 bin_neckSections_ii(trajCount) = 2;
%             elseif X0_n(2) <  0 && X0_n(3) <= 0
%                 bin_neckSections_ii(trajCount) = 3;
%             elseif X0_n(2) <= 0 && X0_n(3) >  0
%                 bin_neckSections_ii(trajCount) = 4;
%             end

            if store_trajectories == 1
%                 trajs = [X_BCR_n(:,1:3); NaN(1,3)];
%                 trajs_ii(trajCount,:) = trajs;
            end
            
            % ---------------------------------------
            % Making sure elevation-90 is only done once in for-loop
            % ---------------------------------------
            if el_i == 90
                elevation90HasBeenDone = 1;
            end
        end
    end
    
    % -------------------------------------------------
    % Storing data from all azimuths/elevations
    % -------------------------------------------------
    r0Data{ii}.X0s          = X0s_ii;
    r0Data{ii}.latLons      = latLons_ii;
    r0Data{ii}.impactSpeeds = impactSpeeds_ii;
    r0Data{ii}.impactAngles = impactAngles_ii;
    
    r0Data{ii}.bin_neckSections = bin_neckSections_ii;
    r0Data{ii}.bin_impactAngles = bin_impactAngles_ii;
    
    if store_trajectories == 1
        r0Data{ii}.trajs   = trajs_ii;
    end
    
    
    % -------------------------------------------------
    % Let the people know
    % -------------------------------------------------
    fprintf('%3.1f Minutes Ellapsed\n',toc(ticWhole)/60)
end

% -------------------------------------------------
% Concatenating data
% -------------------------------------------------
%%% Number of trajectories
n_traj = n_r0s * n_v0s_per_r0;

%%% Preallocate general data
X0s          = zeros(n_traj,6);
latLons      = zeros(n_traj,2);
impactSpeeds = zeros(n_traj,1);
impactAngles = zeros(n_traj,1);

%%% Preallocate bins for neck sections
binData_neckSections(length(bins_neckSections)).latLons = [];

%%% Preallocate bins for impact angles
binData_impactAngles(binCount_ImpactAngles).latLons = [];

%%% Store full trajectories
if store_trajectories == 1
    trajs   = zeros(n_traj,3);
end

%%% Bring data together into single variables
tCount = 0;
for kk = 1:length(r0Data)
    for jj = 1:size(r0Data{kk}.X0s,1)
        %%% Count trajectory
        tCount = tCount + 1;
        
        %%% Store
        X0s(tCount,:)        = r0Data{kk}.X0s(jj,:);
        latLons(tCount,:)    = r0Data{kk}.latLons(jj,:);
        impactSpeeds(tCount) = r0Data{kk}.impactSpeeds(jj);
        impactAngles(tCount) = r0Data{kk}.impactAngles(jj);
        
        %%% Binning appropriate data
        binData_neckSections(r0Data{kk}.bin_neckSections(jj)).latLons =...
            [binData_neckSections(r0Data{kk}.bin_neckSections(jj)).latLons; r0Data{kk}.latLons(jj,:)];
        if isnan(r0Data{kk}.bin_impactAngles(jj)) == 0
            binData_impactAngles(r0Data{kk}.bin_impactAngles(jj)).latLons =...
                [binData_impactAngles(r0Data{kk}.bin_impactAngles(jj)).latLons; r0Data{kk}.latLons(jj,:)];
        end
        
        if store_trajectories == 1
            trajs(tCount,:)  = r0Data{kk}.trajs(jj,:);
        end
    end
end

if tCount ~= n_traj
    warning('Houston we''ve got a problem')
end

% -------------------------------------------------
% 
% -------------------------------------------------

figure(1); hold all
quiver3(r0Data{1}.X0s(:,1),r0Data{1}.X0s(:,2),r0Data{1}.X0s(:,3),...
    r0Data{1}.X0s(:,4),r0Data{1}.X0s(:,5),r0Data{1}.X0s(:,6))
PlotBoi3('X','Y','Z',15)
axis equal


figure; hold all
plot(latLons(:,2),latLons(:,1),'r.','linewidth',1.5,'markersize',13)
xlim([-180 180])
ylim([-90 90])
PlotBoi2('Longitude','Latitude',15)

% if store_trajectories == 1
%     figure(3); hold all
%     plot3(trajs(:,1),trajs(:,2),trajs(:,3))
%     plotBodyTexture3(secondary.R_n, [1-secondary.MR, 0, 0], secondary.img)
%     plot3(L123(Lpoint,1),L123(Lpoint,2),L123(Lpoint,3),'m^','markersize',10)
%     PlotBoi3('X','Y','Z',15)
%     axis equal
%     xlim([L123(1,1)-secondary.R_n, L123(2,1)+secondary.R_n])
%     ylim([-secondary.R_n secondary.R_n].*5)
%     zlim([-secondary.R_n secondary.R_n].*5)
% end

if plot_initialConditions == 1
    %%% Creating ZV contours
    % Contour bounds
    xCont_min = L123(1,1)-2*secondary.R_n;
    xCont_max = L123(2,1)+2*secondary.R_n;
    yCont_min = -9*secondary.R_n;
    yCont_max = 9*secondary.R_n;

    % Creating x-y grid
    xs = linspace(xCont_min,xCont_max,600);
    ys = linspace(yCont_min,yCont_max,200);
    [X_xy, Y_xy] = meshgrid(xs,ys);
%     [Y_yz, Z_yz] = meshgrid(ys,zs);
    clear xs ys zs

    % Calculating JCs across x-y grid
    JCs_xy = zeros(size(X_xy));
    for xk = 1:size(X_xy,1)
        for yk = 1:size(X_xy,2)
            %%% Zero-Velocity Curve
            zv = JacobiConstantCalculator(secondary.MR,[X_xy(xk,yk), Y_xy(xk,yk), 0] ,[0, 0, 0]);
            JCs_xy(xk,yk) = zv;
        end
    end
    
    figure; hold all
    plot3(r0s(:,1),r0s(:,2),r0s(:,3),'.')
    plot3(ones(size(yzContourPoints,2),1).*L123(Lpoint,1), yzContourPoints(1,:),yzContourPoints(2,:),'k','linewidth',3)
    plotBodyTexture3(secondary.R_n, [1-secondary.MR, 0, 0], secondary.img)
    quiver3(X0s(:,1),X0s(:,2),X0s(:,3),X0s(:,4),X0s(:,5),X0s(:,6))
    [xyContourPoints,href] = contour(X_xy,Y_xy,JCs_xy,[JC_scInitial, JC_scInitial],...
        'color',colors.std.black,'linewidth',3);
    PlotBoi3('x$_n$','y$_n$','z$_n$',14,'LaTex')
    view(-35,24)
    camva(9)
    axis equal
end

if plot_fullSystemContour == 1
    %%% Creating ZV contours
    % Contour bounds
    xCont_min = L123(3,1)*1.2;
    xCont_max = L123(2,1)*1.2;
    yCont_min = -1.1;
    yCont_max = 1.1;

    % Creating x-y grid
    xs = linspace(xCont_min,xCont_max,750);
    ys = linspace(yCont_min,yCont_max,750);
    [X_xy, Y_xy] = meshgrid(xs,ys);
    clear xs ys

    % Calculating JCs across x-y grid
    JCs_xy = zeros(size(X_xy));
    for xk = 1:size(X_xy,1)
        for yk = 1:size(X_xy,2)
            %%% Zero-Velocity Curve
            zv = JacobiConstantCalculator(secondary.MR,[X_xy(xk,yk), Y_xy(xk,yk), 0] ,[0, 0, 0]);
            JCs_xy(xk,yk) = zv;
        end
    end
    
    figure; hold all
    plotBodyTexture3(secondary.R_n, [1-secondary.MR, 0, 0], secondary.img)
    plotBodyTexture3(primary.R./rNorm, [-secondary.MR, 0, 0], primary.img)
    [xyContourPoints,href] = contourf(X_xy,Y_xy,JCs_xy,[JC_scInitial, JC_scInitial],...
        'color',colors.std.black,'linewidth',1.5);
    colormap(colors.sch.r6(6,:))
    PlotBoi3('X','Y','Z',14)
    view(0,90)
    camva(9)
    axis equal
end

if plot_sectionColors == 1
    figure; hold all
    
    for kk = binCount_NeckSections+1:-1:2
        plotBody2(bins_neckSections(kk),[0,0,0], binColors_NeckSections(kk-1,:), binColors_NeckSections(kk-1,:), 1)
    end
    
    fillout(yzContourPoints(1,:),yzContourPoints(2,:));
    plot(yzContourPoints(1,:), yzContourPoints(2,:),'k','linewidth',3)
% 
%     eps = 1e-17;
%     pypz = yzContourPoints(:,inpolygon(yzContourPoints(1,:), yzContourPoints(2,:),[0, 100, 100, 0],[eps, eps, 100, 100]));
%     pynz = yzContourPoints(:,inpolygon(yzContourPoints(1,:), yzContourPoints(2,:),[eps, 100, 100, eps],[0, 0, -100, -100]));
%     nynz = yzContourPoints(:,inpolygon(yzContourPoints(1,:), yzContourPoints(2,:),[-100 0 0 -100],[-eps -eps -100 -100]));
%     nypz = yzContourPoints(:,inpolygon(yzContourPoints(1,:), yzContourPoints(2,:),[-100 -eps -eps -100],[100 100 0 0]));
% 
%     pypz = [pypz, [0, 0, y_neck_upper; 0, z_neck_upper, 0]];
%     pynz = [pynz, [0, 0, y_neck_upper; 0, -z_neck_upper, 0]];
%     nynz = [nynz, [-y_neck_upper, 0, 0; 0, -z_neck_upper, 0]];
%     nypz = [nypz, [-y_neck_upper, 0, 0; 0, z_neck_upper, 0]];
% 
%     area(pypz(1,:),pypz(2,:),'FaceColor',binColors_NeckSections(1,:))
%     area(pynz(1,:),pynz(2,:),'FaceColor',binColors_NeckSections(2,:))
%     area(nynz(1,:),nynz(2,:),'FaceColor',binColors_NeckSections(3,:))
%     area(nypz(1,:),nypz(2,:),'FaceColor',binColors_NeckSections(4,:))

    axis equal
    PlotBoi3('Y','Z','X',14)
end

figure; hold all
for kk = 1:4
    if isempty(binData_neckSections(kk).latLons) == 0
        plot(binData_neckSections(kk).latLons(:,2),binData_neckSections(kk).latLons(:,1),...
            '.','markersize',15,'color',binColors_NeckSections(kk,:))
    end
end
title('Neck Sections')
PlotBoi2('Longitude, deg','Latitude, deg',14)
xlim([-180 180])
ylim([-90 90])

%%% Colorbar for whole figure
% a = gcf;
% cbar = colorbar('Position',a.Position);
cbar1 = colorbar;
caxis([0 size(binColors_NeckSections,1)]);
cbar1.FontName     = 'Arial';
cbar1.FontSize     = 10;
cbar1.Ticks        = sort(0:binCount_NeckSections-1);
cbar1.TickLabels = num2cell(bins_neckSections(1:end-1).*rNorm);
cbar1.Label.String = {sprintf('|r_0| from L_%1.0f, km',Lpoint)};
cbar1.Label.Rotation = 0; % 0 = horizontal, 90 = vertical
cbar1.Label.Position = [.2, 4.2, 0];
colormap(binColors_NeckSections)



figure; hold all
for kk = binCount_ImpactAngles:-1:1
if isempty(binData_impactAngles(kk).latLons) == 0
plot(binData_impactAngles(kk).latLons(:,2),binData_impactAngles(kk).latLons(:,1),...
'.','markersize',15,'color',binColors_ImpactAngles(kk,:))
end
end
PlotBoi2('Longitude, $^\circ$','Latitude, $^\circ$',14,'LaTex')
xlim([-180 180])
ylim([-90 90])
%%% Colorbar for whole figure
cbar2 = colorbar;
caxis([0, 90]);
cbar2.FontName     = 'Arial';
cbar2.FontSize     = 10;
cbar2.Ticks        = [bins_impactAngles(1:end-1), 90];
cbar2.TickLabels = num2cell([bins_impactAngles(1:end-1), 90]);
cbar2.Label.String = {'Impact Angle �'};
cbar2.Label.Rotation = 0; % 0 = horizontal, 90 = vertical
cbar2.Label.Position = [.2, 95, 0];
binColors_ImpactAngles_spaced = [...
    repmat(binColors_ImpactAngles(1,:),1,1);...
    repmat(binColors_ImpactAngles(2,:),3,1);...
    repmat(binColors_ImpactAngles(3,:),4,1);...
    repmat(binColors_ImpactAngles(4,:),4,1);...
    repmat(binColors_ImpactAngles(5,:),6,1)];
colormap(binColors_ImpactAngles_spaced)


%     % -------------------------------------------------
%     % Save figure
%     % -------------------------------------------------
%     if saveFigures == 1
%         figure(1)
%         figPath = '/Users/lukebury/Documents/MATLAB/CU/Research/Moon_Landing/Figures/';
%         figName = sprintf('3DMC_%s%sL%1d_%1.1fpi_dvLp%2.2f.png', figPath,...
%             secondaryName(1:3), Lpoint, t_f/pi, dvLp_mps);
%         % ex: 3DMC_EurL1_1.0pi_dvLp56.50.0.png
%         saveas(gcf,figName)
%         close all
%     end

toc(ticWhole)

















