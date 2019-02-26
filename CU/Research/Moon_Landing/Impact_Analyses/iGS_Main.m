clear
% clc
close all

ticWhole = tic;
% ========================================================================
%%% Run Switches
% ========================================================================
%%% Testing?
testCaseOn             = 0;

%%% Zonal harmonics (J21)
on_J21                 = 0;


%%% Set paths based on computer
if isequal(computer,'MACI64')      % Mac
    on_Fortuna         = 0;
    mbinPath = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/mbin';
    moonFuncsPath = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/CU/Research/Moon_Landing/Moon_Landing_funcs';
    savepath = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs';
    computerTag = 'M';
    
elseif isequal(computer,'GLNXA64') % Fortuna
    on_Fortuna         = 1;
    mbinPath = '/home/lubu8198/MatGit/mbin';
    moonFuncsPath = '/home/lubu8198/MatGit/CU/Research/Moon_Landing/Moon_Landing_funcs';
    savepath = '/home/lubu8198/MatGit/MatlabOutputs';
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
addpath(genpath(moonFuncsPath))

% ========================================================================
%%% Setup
% ========================================================================
% -------------------------------------------------
% Choose system
% -------------------------------------------------
%%% General data on solar system bodies
bodies = getBodyData(mbinPath);
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
R1_n = primary.R/rNorm;
R2_n = secondary.R_n;
MR = secondary.MR;
secondaryName = secondary.name;
J21 = primary.J2;
J22 = 0;

% -------------------------------------------------
% Grid-search settings
% -------------------------------------------------
%%% Collinear lagrange point of interest
Lpoint = 2;  % 1 or 2

%%% Acquire Collinear Lagrange points
if on_J21 == 0
    L123 = EquilibriumPoints(secondary.MR,1:3); % [3x3] of L1, L2, L3 normalized BCR coordinates
elseif on_J21 == 1
    L123 = EquilibriumPoints_J2(secondary.MR,J21,J22,R1_n,R2_n,1:3);
end
Lpoint_x = L123(Lpoint,1);

%%% How fast the SC would be traveling over the Lagrange point
% dvLp_mps = 200; % Meters per second - Europa
% dvLp_mps = 50; % Meters per second - Enceladus
% dvLp_mps = 50; % Meters per second

for dvLp_mps = [150] % Europa
% for dvLp_mps = [50, 100, 150, 200, 250, 300, 350] % Europa
% for dvLp_mps = [13, 26, 39, 52, 65, 78, 91] % Enceladus
% for dvLp_mps = [58, 116, 174, 232, 290, 348, 406] % Titan








if isequal(secondary.name,'europa')
    %%% Spacing of initial positions within 3D neck
    r0GridSpacing_km = 50; % km
    %%% Spacing between azimuths and elevations of v0s per r0
    n_v0s_per_r0_target = 145;
elseif isequal(secondary.name,'enceladus')
    %%% Spacing of initial positions within 3D neck
    r0GridSpacing_km = 4; % km
    %%% Spacing between azimuths and elevations of v0s per r0
    n_v0s_per_r0_target = 145;
elseif isequal(secondary.name,'titan')
    %%% Spacing of initial positions within 3D neck
    r0GridSpacing_km = 50*5; % km
    %%% Spacing between azimuths and elevations of v0s per r0
    n_v0s_per_r0_target = 145;
end
% n_target = 1297; % ~ old  5 deg spacing
% n_target = 352;  % ~ old 10 deg spacing
% n_target = 145;  % ~ old 15 deg spacing - should be used
% n_target = 37;   % ~ old 30 deg spacing - decent
% n_target = 17;   % ~ old 45 deg spacing
% n_target = 5;    % ~ old 90 deg spacing

if testCaseOn == 1
%     %%% To get a couple of low impacts ... 
%     r0GridSpacing_km = 2000;
%     n_v0s_per_r0_target = 350;
    
    %%% To get a couple of low impacts ... 
    if isequal(secondary.name,'europa')
        r0GridSpacing_km = 500;
        n_v0s_per_r0_target = 5;
    elseif isequal(secondary.name,'enceladus')
        r0GridSpacing_km = 75;
        n_v0s_per_r0_target = 30;
    elseif isequal(secondary.name,'titan')
        r0GridSpacing_km = 3500;
        n_v0s_per_r0_target = 10;
    end
end

%%% Selecting time vector
t_i = 0; % sec
t_f = 4*pi; % Long bc events are watching for impact or escape
% changing to t_f = 6*pi give 0.1% more impacts out of 5940928
% trajectories.... so we're gonna ignore it for now
dt = t_f/10000;

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
if on_J21 == 0
    [JC_Lp] = JacobiConstantCalculator(secondary.MR,L123(Lpoint,:),[0,0,0]);
elseif on_J21 == 1
    [JC_Lp] = JacobiConstantCalculator_J2(secondary.MR,L123(Lpoint,:),[0,0,0], R1_n, R2_n, J21, J22);
end

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
r1 = @(x,y,z) sqrt((x+u)^2+y.^2+z^2);
r2 = @(x,y,z) sqrt((x-1+u)^2+y.^2+z^2);
if on_J21 == 0
    f = @(y) x^2 + y.^2 - c + 2*(1-u)./r1(x,y,z) + 2*u./r2(x,y,z);
elseif on_J21 == 1
    f = @(y) x^2 + y.^2 - c + 2*(1-u)./r1(x,y,z) + 2*u./r2(x,y,z) + R1_n*R1_n*J21*(1-u)*(r1(x,y,z)^2 - 3*z*z)/(r1(x,y,z)^5);
end


%%% Find the root of this function in the appropriate range
y_neckRange = 15*secondary.R_n;
z_neckRange = 15*secondary.R_n;
if isequal(secondary.name,'titan') == 1
    y_neckRange = 60*secondary.R_n;
    z_neckRange = 30*secondary.R_n;
end
    y_neck_upper = fzero(f,[0 y_neckRange]);



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
r1 = @(x,y,z) sqrt((x+u)^2+y.^2+z^2);
r2 = @(x,y,z) sqrt((x-1+u)^2+y.^2+z^2);

if on_J21 == 0
    f = @(z) x^2 + y.^2 - c + 2*(1-u)./r1(x,y,z) + 2*u./r2(x,y,z);
elseif on_J21 == 1
    f = @(z) x^2 + y.^2 - c + 2*(1-u)./r1(x,y,z) + 2*u./r2(x,y,z) + R1_n*R1_n*J21*(1-u)*(r1(x,y,z)^2 - 3*z*z)/(r1(x,y,z)^5);
end

%%% Find the root of this function in the appropriate range
% % if isequal(secondary.name,'titan') == 0
z_neck_upper = fzero(f,[0 z_neckRange]);
% % end

%%% Clear variables
clear c u x f y

% -------------------------------------------------
% Create grid of starting locations based on y-z neck to find contour
% points
% -------------------------------------------------
if isequal(secondary.name,'titan') == 0
    ys = linspace(-2*y_neck_upper, 2*y_neck_upper, 400);
    zs = linspace(-2*z_neck_upper, 2*z_neck_upper, 400);
elseif isequal(secondary.name,'titan') == 1
    scale = 1.5;
    ys = linspace(-y_neckRange*scale, y_neckRange*scale, 400);
    zs = linspace(-z_neckRange*scale, z_neckRange*scale, 400);
end
[Y_yz,Z_yz] = meshgrid(ys,zs);

% -------------------------------------------------
% Only keep starting positions that are valid for the energy level
% -------------------------------------------------
%%% Calculating JCs across y-z grid
JCs_yz_Lpoint = zeros(size(Y_yz));
for yk = 1:size(Y_yz,1)
    for zk = 1:size(Y_yz,2)
        %%% Zero-Velocity Curve
        if on_J21 == 0
            zv = JacobiConstantCalculator(secondary.MR,[L123(Lpoint,1), Y_yz(yk,zk), Z_yz(yk,zk)] ,[0, 0, 0]);
        elseif on_J21 == 1
            zv = JacobiConstantCalculator_J2(secondary.MR,[L123(Lpoint,1), Y_yz(yk,zk), Z_yz(yk,zk)],[0,0,0], R1_n, R2_n, J21, J22);
        end
        JCs_yz_Lpoint(yk,zk) = zv;
    end
end

%%% Get points of y-z contour in 3D space
[ yzContourPoints4 ] = getContourPoints( Y_yz, Z_yz, JCs_yz_Lpoint, JC_scInitial );

% -------------------------------------------------
% Create grid of starting locations to sample r0s from
% -------------------------------------------------
if isequal(secondary.name,'titan') == 0
    n_r0Grid_y = 2*(2*y_neck_upper) / (r0GridSpacing_km/rNorm);
    n_r0Grid_z = 2*(2*z_neck_upper) / (r0GridSpacing_km/rNorm);
    ys = linspace(-2*y_neck_upper, 2*y_neck_upper, n_r0Grid_y);
    zs = linspace(-2*z_neck_upper, 2*z_neck_upper, n_r0Grid_z);
elseif isequal(secondary.name,'titan') == 1
    n_r0Grid_y = 2*(y_neckRange) / (r0GridSpacing_km/rNorm);
    n_r0Grid_z = 2*(2*z_neckRange) / (r0GridSpacing_km/rNorm);
    ys = linspace(-2*y_neckRange, 2*y_neckRange, n_r0Grid_y);
    zs = linspace(-2*z_neckRange, 2*z_neckRange, n_r0Grid_z);
end
[Y_yz_r0s,Z_yz_r0s] = meshgrid(ys,zs);

%%% Search grid for points lying within y-z contour and keep those as r0s
r0s = [];
for yk = 1:size(Y_yz_r0s,1)
    for zk = 1:size(Z_yz_r0s,2)
        if inpolygon(Y_yz_r0s(yk,zk), Z_yz_r0s(yk,zk),yzContourPoints4(1,:),yzContourPoints4(2,:)) == 1
            if Z_yz_r0s(yk,zk) >= 0 % dynamics are symmetric about z=0, so only need to check one side of that
                r0s = [r0s; L123(2,1), Y_yz_r0s(yk,zk), Z_yz_r0s(yk,zk)];
            end
        end
    end
end


%%% Total number of r0s
n_r0s = size(r0s,1);

% -------------------------------------------------
% Creating Bins for impact angle and neck sections
% -------------------------------------------------
%%% Bins for impact angles
% % Making last bin too big so everything fits in bins(1:end-1)
bins_impactAngles = [0, 5, 20, 40, 60, 91];

%%% Bins for neck sections
% smallNeck = smallerValue(z_neck_upper, y_neck_upper);
% largeNeck = largerValue(z_neck_upper, y_neck_upper);
% bins_neckSections = [linspace(0,smallNeck,binCount_NeckSections), largeNeck];
bins_neckSectionScalars = [0.25, 0.5, 0.75, 1];

% -------------------------------------------------
% Making polygons to bin initial neck position
% -------------------------------------------------
yzContourPoints3 = yzContourPoints4.*bins_neckSectionScalars(3);
yzContourPoints2 = yzContourPoints4.*bins_neckSectionScalars(2);
yzContourPoints1 = yzContourPoints4.*bins_neckSectionScalars(1);

% -------------------------------------------------
% Velocity vectors at each point
% -------------------------------------------------
%%% Acquiring matrix of vHat vectors
if Lpoint == 1
    [vHats2] = vHatHemisphere(n_v0s_per_r0_target,'x');
elseif Lpoint == 2
    [vHats2] = vHatHemisphere(n_v0s_per_r0_target,'-x');
end

%%% For reference, total number of trajectories to simulate
n_v0s_per_r0 = size(vHats2,1);
n_traj = n_r0s * n_v0s_per_r0;

% -------------------------------------------------
% Predefining some things before parfor
% -------------------------------------------------
%%% Setting time vector
time0_n = t_i:dt:t_f;

%%% Choosing ode tolerance
tol = 1e-13;

%%% Setting integrator options
options = odeset('RelTol',tol,'AbsTol',tol);
options_ImpactEscape = odeset('Events',@event_ImpactEscape_CR3Bn,'RelTol',tol,'AbsTol',tol);

% ========================================================================
%%% Loop through conditions, store results in bins, and plot
% ========================================================================
r0Data = {};

parfor ii = 1:n_r0s
    
    ticLoop = tic;
    % -------------------------------------------------
    % Reducing broadcast variables
    % -------------------------------------------------
    L123mat = L123;
    vHats = vHats2;
    yzContourPoints11 = yzContourPoints1;
    yzContourPoints22 = yzContourPoints2;
    yzContourPoints33 = yzContourPoints3;
    yzContourPoints44 = yzContourPoints4;
    
    % -------------------------------------------------
    % Setting parameters for integration
    % -------------------------------------------------
    prms = struct();
    prms.u = MR;
    prms.R2_n = R2_n;
    prms.L1x = L123mat(1,1);
    prms.L2x = L123mat(2,1);
    
    if on_J21 == 1
        prms.R1_n = R1_n;
        prms.J21  = J21;
    end

    % -------------------------------------------------
    % Setting current initial position
    % -------------------------------------------------
    r0_i = r0s(ii,:);

    % -------------------------------------------------
    % Finding necessary v0 to maintain correct JC0
    % -------------------------------------------------
    %%% With JC0 defined, starting velocity is a function of position. So first
    %%% we must calculate the JC of the stationary starting position
    JC_initialPos = [];
    if on_J21 == 0
        JC_initialPos = JacobiConstantCalculator(MR,r0_i,[0,0,0]);
    elseif on_J21 == 1
        JC_initialPos = JacobiConstantCalculator_J2(MR,r0_i, [0,0,0], R1_n, R2_n, J21, J22);
    end
    
    %%% Starting velocity is found from difference between s/c JC (JC_scDesired) and the
    %%% JC of the stationary starting position (JC_initialPos)
    dJC_forInitialVelocity = JC_initialPos - JC_scInitial;
    
    %%% Setting variables so they're not considered 'temporary variables'
    v0i_mag = [];
    
    %%% Find necessary velocity magnitude
    if dJC_forInitialVelocity < 0
        warning('Spacecraft starting in a forbidden zone')
    elseif dJC_forInitialVelocity >= 0
        v0i_mag = sqrt(abs(dJC_forInitialVelocity));
    end
    
    % -------------------------------------------------
    % Looping through velocity vectors
    % -------------------------------------------------    
    %%% Initializing data matrices
    X0s_ii              = zeros(n_v0s_per_r0,6);
    latLons_ii          = zeros(n_v0s_per_r0,2);
    trajs_ii            = zeros(n_v0s_per_r0,3);
    impactSpeeds_ii     = zeros(n_v0s_per_r0,1);
    impactAngles_ii     = zeros(n_v0s_per_r0,1);
    endTimes_ii         = zeros(n_v0s_per_r0,1);
    bin_impactAngles_ii = zeros(n_v0s_per_r0,1);
    bin_neckSections_ii = zeros(n_v0s_per_r0,1);
    trajIDs_ii          = zeros(n_v0s_per_r0,1);
    v0AzEls_ii          = zeros(n_v0s_per_r0,2);
    maxLats_ii          = zeros(n_v0s_per_r0,1);

    %%% Running loop
    for vi = 1:n_v0s_per_r0
        % ---------------------------------------
        % Preallocating temporary variables
        % ---------------------------------------
        trajID       = [];
        impactLatLon = [];
        impactSpeed  = [];
        impactAngle  = [];
        endTime      = [];
        v0AzEl       = [];
        maxLat       = [];
        
        %%% Compute unique trajectory identification
        trajID = (ii-1)*n_v0s_per_r0 + vi;

        %%% Setting v0
        v0_i = vHats(vi,:) .* v0i_mag;
        
        %%% Setting X0
        X0_n = [r0_i, v0_i];

        % ---------------------------------------
        % Integrating
        % ---------------------------------------
        %%% Propagating trajectory
        time_n            = [];
        X_BCR_n           = [];
        time_eventImpact  = [];
        X_eventImpact     = [];
        index_eventImpact = [];
        if on_J21 == 0
            [time_n, X_BCR_n, time_eventImpact, X_eventImpact, index_eventImpact] = ode113(@Int_CR3Bn,...
                time0_n, X0_n, options_ImpactEscape, prms);
        elseif on_J21 == 1
            [time_n, X_BCR_n, time_eventImpact, X_eventImpact, index_eventImpact] = ode113(@Int_CR3Bn_ZH,...
                time0_n, X0_n, options_ImpactEscape, prms);
        end
        
        % ---------------------------------------
        % Determining end conditions
        % ---------------------------------------
        impact1_escape2_orbit3 = []; % 1-impact, 2-escape, 3-still in orbit
        if isempty(X_eventImpact) == 0
            if abs(norm(X_eventImpact(1:3)-[1-MR,0,0])-R2_n) < 1e-9 % impact
                impact1_escape2_orbit3 = 1;
            elseif X_eventImpact(1) <= prms.L1x || X_eventImpact(1) >= prms.L2x % escape
                impact1_escape2_orbit3 = 2;
            end
        elseif isempty(X_eventImpact) == 1 % No event yet, still in orbit
            impact1_escape2_orbit3 = 3;
        end
        
        if impact1_escape2_orbit3 == 2
            if X_eventImpact(1) <= prms.L1x
                format long
                time_n(end)
                X0_n
                fprintf('ttt')
            end
        end
        
        % ---------------------------------------
        % Determining azimuth and elevation of v0
        % ---------------------------------------
        [Az,El,Rad] = cart2sph(X0_n(4), X0_n(5), X0_n(6));
        v0AzEl = [Az, El];
        
        % ---------------------------------------
        % Determining lat/lon, speed, and angle of impact
        % ---------------------------------------
        if impact1_escape2_orbit3 == 1
            % --------------------------
            % Impact lat/lon
            % --------------------------
            %%% Creating SCR position
            rImpact_SCR = X_eventImpact(1,1:3) - [1-MR, 0, 0];

            %%% Finding lat/lon of impact site
%             [lat_deg, lon_deg] = surfCR3BP2latlon(rImpact_SCR, 'secondary', MR)
            [lat_deg, lon_deg] = BCR2latlon(X_eventImpact(1,1:3), 'secondary', MR)
%             [lat, lon] = ECEF2latlon(rImpact_SCR,'degrees','stupidMoon');
            impactLatLon = [lat_deg, lon_deg];
            
            
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
            [impactAngle] = calcImpactAngle(rImpact_SCR_n,vHatImpact_n,'degrees');
            
            %%% Storing endTime
            endTime = time_n(end);
            
            % --------------------------
            % Maximum latitude during trajectory
            % --------------------------
            maxLat = 0;
            for kk = 1:size(X_BCR_n,1)
%                 r_kk_SCR = X_BCR_n(kk,1:3) - [1-MR, 0, 0];
%                 [lat_kk, lon_kk] = ECEF2latlon(r_kk_SCR,'degrees','stupidMoon');
                [lat_kk_deg, lon_kk_deg] = BCR2latlon(X_BCR_n(kk,1:3), 'secondary', MR);
                if abs(lat_kk_deg) > maxLat
                    maxLat = abs(lat_kk_deg);
                end
            end
            
            
        elseif impact1_escape2_orbit3 == 2 % escape
            impactLatLon = [NaN, NaN];
            impactSpeed  = NaN;
            impactAngle  = NaN;     
            endTime      = NaN;
            maxLat       = NaN;
        elseif impact1_escape2_orbit3 == 3 % orbit
            impactLatLon = [NaN, NaN];
            impactSpeed  = NaN;
            impactAngle  = NaN; 
            endTime      = NaN;
            maxLat       = NaN;
        end
            
        % =======================================
        % Storing results of trajectory
        % ---------------------------------------
        trajIDs_ii(vi)      = trajID;
        X0s_ii(vi,:)        = X0_n;
        latLons_ii(vi,:)    = impactLatLon;
        impactSpeeds_ii(vi) = impactSpeed;
        impactAngles_ii(vi) = impactAngle;
        endTimes_ii(vi)     = endTime;
        v0AzEls_ii(vi,:)    = v0AzEl;
        maxLats_ii(vi)      = maxLat;

        % ---------------------------------------
        % Assigning impact-angle bin
        % ---------------------------------------
        %%% Bin - Impact Angle
        bin_impactAngles_ii(vi) = discretize(impactAngle, bins_impactAngles);
        
        % ---------------------------------------
        % Assigning neck-section bin
        % ---------------------------------------
        %%% Bin - Neck Section
        if inpolygon(X0_n(2),X0_n(3),yzContourPoints11(1,:),yzContourPoints11(2,:)) == 1
            bin_neckSections_ii(vi) = 1;
        elseif inpolygon(X0_n(2),X0_n(3),yzContourPoints22(1,:),yzContourPoints22(2,:)) == 1
            bin_neckSections_ii(vi) = 2;
        elseif inpolygon(X0_n(2),X0_n(3),yzContourPoints33(1,:),yzContourPoints33(2,:)) == 1
            bin_neckSections_ii(vi) = 3;
        elseif inpolygon(X0_n(2),X0_n(3),yzContourPoints44(1,:),yzContourPoints44(2,:)) == 1
            bin_neckSections_ii(vi) = 4;
        else % In case it was just missed by the resolution of the outer-most polygon
            bin_neckSections_ii(vi) = 4;
        end
            
%         if inpolygon(Y_yz(yk,zk), Z_yz(yk,zk),yzContourPoints4(1,:),yzContourPoints4(2,:)) == 1
    end
    
    % -------------------------------------------------
    % Storing data from all azimuths/elevations
    % -------------------------------------------------
    r0Data{ii}.trajIDs          = trajIDs_ii;
    r0Data{ii}.X0s              = X0s_ii;
    r0Data{ii}.latLons          = latLons_ii;
    r0Data{ii}.impactSpeeds     = impactSpeeds_ii;
    r0Data{ii}.impactAngles     = impactAngles_ii;
    r0Data{ii}.endTimes         = endTimes_ii;
    r0Data{ii}.bin_neckSections = bin_neckSections_ii;
    r0Data{ii}.bin_impactAngles = bin_impactAngles_ii;
    r0Data{ii}.v0AzEls          = v0AzEls_ii; 
    r0Data{ii}.maxLats          = maxLats_ii;
    
end


% -------------------------------------------------
% Writing data to CSV
% -------------------------------------------------
%%% File Names
if testCaseOn == 0
    if on_J21 == 0
        %%% All impact data
        filename_allTraj = fullfile(savepath, sprintf('%s.iGS_%sL%1.0f_%1.0fmps_%1.0fkm_%1.0fv0s_data.txt',...
            computerTag,secondary.name(1:3),Lpoint,dvLp_mps,r0GridSpacing_km,n_v0s_per_r0));
        %%% Low-impact-angle data
        filename_landingTraj = fullfile(savepath, sprintf('%s.iGS_%sL%1.0f_%1.0fmps_%1.0fkm_%1.0fv0s_land.txt',...
            computerTag,secondary.name(1:3),Lpoint,dvLp_mps,r0GridSpacing_km,n_v0s_per_r0));
        %%% Log file
        filename_runData = fullfile(savepath, sprintf('%s.iGS_%sL%1.0f_%1.0fmps_%1.0fkm_%1.0fv0s_log.txt',...
            computerTag,secondary.name(1:3),Lpoint,dvLp_mps,r0GridSpacing_km,n_v0s_per_r0));
    elseif on_J21 == 1
        %%% All impact data
        filename_allTraj = fullfile(savepath, sprintf('%s.iGS_%sL%1.0f_J21_%1.0fmps_%1.0fkm_%1.0fv0s_data.txt',...
            computerTag,secondary.name(1:3),Lpoint,dvLp_mps,r0GridSpacing_km,n_v0s_per_r0));
        %%% Low-impact-angle data
        filename_landingTraj = fullfile(savepath, sprintf('%s.iGS_%sL%1.0f_J21_%1.0fmps_%1.0fkm_%1.0fv0s_land.txt',...
            computerTag,secondary.name(1:3),Lpoint,dvLp_mps,r0GridSpacing_km,n_v0s_per_r0));
        %%% Log file
        filename_runData = fullfile(savepath, sprintf('%s.iGS_%sL%1.0f_J21_%1.0fmps_%1.0fkm_%1.0fv0s_log.txt',...
            computerTag,secondary.name(1:3),Lpoint,dvLp_mps,r0GridSpacing_km,n_v0s_per_r0));
    end
elseif testCaseOn == 1
    if on_J21 == 0
        %%% All impact data
        filename_allTraj     = fullfile(savepath,'testFile_data.txt');
        %%% Low-impact-angle data
        filename_landingTraj = fullfile(savepath,'testFile_land.txt');
        %%% Log file
        filename_runData     = fullfile(savepath,'testFile_log.txt');
    elseif on_J21 == 1
        %%% All impact data
        filename_allTraj     = fullfile(savepath,'testFile_J21_data.txt');
        %%% Low-impact-angle data
        filename_landingTraj = fullfile(savepath,'testFile_J21_land.txt');
        %%% Log file
        filename_runData     = fullfile(savepath,'testFile_J21_log.txt');
    end
end

% %%% Opening files and writing header
% f_allTraj = fopen(filename_allTraj, 'wt');
% fprintf(f_allTraj,'bin_impactAngle,bin_neckSection,latitude,longitude,impactAngle,endTime\n');  % header
% 
% f_landingTraj = fopen(filename_landingTraj, 'wt');
% fprintf(f_landingTraj,'x0_n,y0_n,z0_n,dx0_n,dy0_n,dz0_n,bin_neckSection,latitude,longitude,endTime\n');  % header

%%% Opening files and writing header
f_allTraj = fopen(filename_allTraj, 'wt');
fprintf(f_allTraj,'trajID,bin_impactAngle,bin_neckSection,latitude,longitude,impactAngle,endTime,y0_n,z0_n,v0_Azimuth,v0_Elevation,maxLatitude\n');  % header

f_landingTraj = fopen(filename_landingTraj, 'wt');
fprintf(f_landingTraj,'trajID,x0_n,y0_n,z0_n,dx0_n,dy0_n,dz0_n,bin_neckSection,latitude,longitude,endTime,maxLatitude\n');  % header

%%% Writing data (and finding maximum latitudes)
maxLat    = 0;
maxLowLat = 0;
lowImpactAngleCounter = 0;
for kk = 1:n_r0s
    %%% If there were impacts from this r0
    if isempty(r0Data{kk}) == 0
        for jj = 1:n_v0s_per_r0
            %%% If this trajectory impacted, write impact data
            if isnan(r0Data{kk}.bin_impactAngles(jj)) == 0
                fprintf(f_allTraj,'%1d,%1d,%1d,%2.7f,%2.7f,%2.1f,%1.5f,%1.16f,%1.16f,%1.6f,%1.6f,%2.7f\n',r0Data{kk}.trajIDs(jj),r0Data{kk}.bin_impactAngles(jj),...
                    r0Data{kk}.bin_neckSections(jj),r0Data{kk}.latLons(jj,1),r0Data{kk}.latLons(jj,2),...
                    r0Data{kk}.impactAngles(jj),r0Data{kk}.endTimes(jj),r0Data{kk}.X0s(jj,2:3),r0Data{kk}.v0AzEls(jj,1:2),...
                    r0Data{kk}.maxLats(jj));
                if abs(maxLat) < abs(r0Data{kk}.latLons(jj,1))
                    maxLat = abs(r0Data{kk}.latLons(jj,1));
                end
                %%% If this was a low-impact-angle trajectory
                if r0Data{kk}.bin_impactAngles(jj) == 1
                    lowImpactAngleCounter = lowImpactAngleCounter + 1;
                    fprintf(f_landingTraj,'%1d,%1.16f,%1.16f,%1.16f,%1.16f,%1.16f,%1.16f,%1d,%2.7f,%2.7f,%1.5f,%2.7f\n',r0Data{kk}.trajIDs(jj),r0Data{kk}.X0s(jj,:),...
                        r0Data{kk}.bin_neckSections(jj),r0Data{kk}.latLons(jj,1),r0Data{kk}.latLons(jj,2),r0Data{kk}.endTimes(jj),...
                        r0Data{kk}.maxLats(jj));
                    if abs(maxLowLat) < abs(r0Data{kk}.latLons(jj,1))
                        maxLowLat = abs(r0Data{kk}.latLons(jj,1));
                    end
                end

            end
        end
    end

end
%%% Just making the file not empty
if lowImpactAngleCounter == 0
    fprintf(f_landingTraj,',\n');
end

%%% Close file
fclose(f_allTraj);
fclose(f_landingTraj);

%%% Preparing to write log file
clear time0_n 
finalToc = toc(ticWhole);
coreNumber = feature('numcores');
% -------------------------------------------------
% Writing log info (basically extra necessary data)
% -------------------------------------------------
%%% Writing log file
f_runData = fopen(filename_runData, 'wt');
fprintf(f_runData,'Run Time (sec):\t\t%1.1f\n',finalToc);
fprintf(f_runData,'Run Time (min):\t\t%1.1f\n',finalToc/60);
fprintf(f_runData,'number of r0s:\t\t%1.0f\n',n_r0s);
fprintf(f_runData,'v0s per r0s:\t\t%1.0f\n',n_v0s_per_r0);
fprintf(f_runData,'total trajs:\t\t%1.0f\n',n_traj);
fprintf(f_runData,'primary:%s\n',primary.name);
fprintf(f_runData,'secondary:%s\n',secondary.name);
fprintf(f_runData,'bins_impactAngles:%1.0f,%1.0f,%1.0f,%1.0f,%1.0f,%1.0f\n',bins_impactAngles);
fprintf(f_runData,'bins_neckSectionScalars:%1.15f,%1.15f,%1.15f,%1.15f\n',bins_neckSectionScalars);
fprintf(f_runData,'dvLp_mps:%1.0f\n',dvLp_mps);
fprintf(f_runData,'JC_scInitial:%1.15f\n',JC_scInitial);
fprintf(f_runData,'Lpoint:%1.0d\n',Lpoint);
fprintf(f_runData,'y_neck_upper:%1.15f\n',y_neck_upper);
fprintf(f_runData,'z_neck_upper:%1.15f\n',z_neck_upper);
fprintf(f_runData,'maxLat:%1.2f\n',maxLat);
fprintf(f_runData,'maxLowLat:%1.2f\n',maxLowLat);

%%% Close file
fclose(f_runData);

finalToc = toc(ticWhole);

end % dvLp_mps = [50, 100, 150, 200, 250, 300, 350]
















