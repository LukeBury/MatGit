clear
clc
close all
mbinPath = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/mbin';
moonFuncsPath = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/CU/Research/Moon_Landing/Moon_Landing_funcs';
polarTrajDataPath = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/CU/Research/Moon_Landing/Impact_Analyses/polarTrajectoryConditions';
addpath(genpath(mbinPath))
addpath(genpath(moonFuncsPath))
tic

% ========================================================================
%%% Run Switches
% ========================================================================
run_newPolarTrajs = 0;

% ========================================================================
%%% Free variables
% ========================================================================
%%% How fast the SC would be traveling over L2
dvLp_mps = 350; % Meters per second

%%% Latitude cutoff (this and above)
latCutoff_deg = 80; % degrees

%%% Timing info
n_dt = 10000;
t_f = 4*pi;

%%% Spacing of grid in neck
r0GridSpacing_km = 1000; % km

%%% Spacing between azimuths and elevations of v0s per r0
n_v0s_per_r0_target = 75;

%%% Load file
polarTrajFile = [polarTrajDataPath,''];
% ========================================================================
%%% Importing Data
% ========================================================================
%%% General data on solar system bodies
bodies = getBodyData(mbinPath);

%%% Color options/schemes
colors = get_colors();

% ========================================================================
%%% Setup
% ========================================================================
% -------------------------------------------------
% Choose system
% -------------------------------------------------
%%% 3B system
primary   = bodies.jupiter;
secondary = bodies.europa;

%%% Normalizing constants
rNorm = secondary.a;         % n <-> km
tNorm = 1/secondary.meanMot; % n <-> sec
vNorm = rNorm / tNorm;       % n <-> km/sec

%%% Lagrange Points
L123 = EquilibriumPoints(secondary.MR,1:3); % [3x3] of L1, L2, L3 normalized BCR coordinates

%%% Setting new variables because parfor is picky
R2_n = secondary.R_n;
MR = secondary.MR;
secondaryName = secondary.name;

% -------------------------------------------------
% Finding initial JC of spacecraft
% -------------------------------------------------
%%% Jacobi constant of Lagrange point
[JC_Lp] = JacobiConstantCalculator(secondary.MR,L123(2,:),[0,0,0]);

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
x = L123(2,1);
z = 0;

%%% JC function equal to zero
f = @(y) x^2 + y.^2 - c + 2*(1-u)./sqrt((x+u)^2+y.^2+z^2) + 2*u./sqrt((x-1+u)^2+y.^2+z^2);

%%% Find the root of this function in the appropriate range
y_neck_upper = fzero(f,[0 10*secondary.R_n]);

%%% Clear variables
clear c u x f z

% -------------------------------------------------
% Finding upper-z-value of neck at L-point
% -------------------------------------------------
%%% Set values
c = JC_scInitial;
u = secondary.MR;
x = L123(2,1);
y = 0;

%%% JC function equal to zero
f = @(z) x^2 + y.^2 - c + 2*(1-u)./sqrt((x+u)^2+y^2+z.^2) + 2*u./sqrt((x-1+u)^2+y^2+z.^2);

%%% Find the root of this function in the appropriate range
z_neck_upper = fzero(f,[0 10*secondary.R_n]);

%%% Clear variables
clear c u x f y

% -------------------------------------------------
% Create grid of starting locations based on y-z neck to find contour
% points
% -------------------------------------------------
ys = linspace(-2*y_neck_upper, 2*y_neck_upper, 400);
zs = linspace(-2*z_neck_upper, 2*z_neck_upper, 400);
[Y_yz,Z_yz] = meshgrid(ys,zs);

% -------------------------------------------------
% Only keep starting positions that are valid for the energy level
% -------------------------------------------------
%%% Calculating JCs across y-z grid
JCs_yz_Lpoint = zeros(size(Y_yz));
for yk = 1:size(Y_yz,1)
    for zk = 1:size(Y_yz,2)
        %%% Zero-Velocity Curve
        zv = JacobiConstantCalculator(secondary.MR,[L123(2,1), Y_yz(yk,zk), Z_yz(yk,zk)] ,[0, 0, 0]);
        JCs_yz_Lpoint(yk,zk) = zv;
    end
end

%%% Get points of y-z contour in 3D space
[ yzContourPoints4 ] = getContourPoints( Y_yz, Z_yz, JCs_yz_Lpoint, JC_scInitial );

% -------------------------------------------------
% Create grid of starting locations to sample r0s from
% -------------------------------------------------
n_r0Grid_y = 2*(2*y_neck_upper) / (r0GridSpacing_km/rNorm);
n_r0Grid_z = 2*(2*z_neck_upper) / (r0GridSpacing_km/rNorm);
ys = linspace(-2*y_neck_upper, 2*y_neck_upper, n_r0Grid_y);
zs = linspace(-2*z_neck_upper, 2*z_neck_upper, n_r0Grid_z);
[Y_yz_r0s,Z_yz_r0s] = meshgrid(ys,zs);

%%% Search grid for points lying within y-z contour and keep those as r0s
r0s = [];
for yk = 1:size(Y_yz_r0s,1)
    for zk = 1:size(Z_yz_r0s,2)
        if inpolygon(Y_yz_r0s(yk,zk), Z_yz_r0s(yk,zk),yzContourPoints4(1,:),yzContourPoints4(2,:)) == 1
            r0s = [r0s; L123(2,1), Y_yz_r0s(yk,zk), Z_yz_r0s(yk,zk)];
        end
    end
end

%%% Total number of r0s
n_r0s = size(r0s,1);

% -------------------------------------------------
% Velocity vectors at each point
% -------------------------------------------------
%%% Acquiring matrix of vHat vectors
[vHats2] = vHatHemisphere(n_v0s_per_r0_target,'-x');

%%% For reference, total number of trajectories to simulate
n_v0s_per_r0 = size(vHats2,1);
n_traj = n_r0s * n_v0s_per_r0;

% -------------------------------------------------
% Predefining some things before parfor
% -------------------------------------------------
%%% Setting time vector
t_i = 0;
time0_n = linspace(t_i,t_f,n_dt);

%%% Choosing ode tolerance
tol = 1e-13;

%%% Setting integrator options
% options = odeset('RelTol',tol,'AbsTol',tol);
options_ImpactEscape = odeset('Events',@event_ImpactEscape_CR3Bn,'RelTol',tol,'AbsTol',tol);


% ========================================================================
%%% Loop through conditions, store results - or load results
% ========================================================================
polarR0Data = {};

parfor ii = 1:n_r0s
    % -------------------------------------------------
    % Reducing broadcast variables
    % -------------------------------------------------
    L123mat = L123;
    vHats = vHats2;
    
    % -------------------------------------------------
    % Setting parameters for integration
    % -------------------------------------------------
    prms = struct();
    prms.u = MR;
    prms.R2_n = R2_n;
    prms.L1x = L123mat(1,1);
    prms.L2x = L123mat(2,1);

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
    Xfs_ii              = zeros(n_v0s_per_r0,6);
    latLons_ii          = zeros(n_v0s_per_r0,2);
    endTimes_ii         = zeros(n_v0s_per_r0,1);
        
    %%% Running loop
    for vi = 1:n_v0s_per_r0

        %%% Setting v0
        v0_i = vHats(vi,:) .* v0i_mag;
        
        %%% Setting X0
        X0_n = [r0_i, v0_i];

        % ---------------------------------------
        % Integrating
        % ---------------------------------------
        %%% Propagating trajectory
        [time_n, X_BCR_n, time_eventImpact, X_eventImpact, index_eventImpact] = ode113(@Int_CR3Bn,...
            time0_n, X0_n, options_ImpactEscape, prms);
        
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
        
        
        % ---------------------------------------
        % Preallocating temporary variables
        % ---------------------------------------
        impactLatLon = [];
        endTime      = [];

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
            [lat, lon] = ECEF2latlon(rImpact_SCR,'degrees','stupidMoon');
            impactLatLon = [lat, lon];

            %%% End time
            endTime = time_eventImpact(end);

        elseif impact1_escape2_orbit3 == 2 % escape
            impactLatLon = [NaN, NaN];
            endTime      = NaN;
        elseif impact1_escape2_orbit3 == 3 % orbit
            impactLatLon = [NaN, NaN];
            endTime      = NaN;
        end
            
        % =======================================
        % Storing results of trajectory
        % ---------------------------------------
        if impactLatLon(1,1) >= latCutoff_deg
            X0s_ii(vi,:)        = X0_n;
            Xfs_ii(vi,:)        = X_BCR_n(end,:);
            latLons_ii(vi,:)    = impactLatLon;
            endTimes_ii(vi)     = endTime;
        else
            X0s_ii(vi,:)        = [NaN, NaN, NaN, NaN, NaN, NaN];
            Xfs_ii(vi,:)        = [NaN, NaN, NaN, NaN, NaN, NaN];
            latLons_ii(vi,:)    = [NaN, NaN];
            endTimes_ii(vi)     = NaN;
        end

    end
    
    % -------------------------------------------------
    % Storing data from all azimuths/elevations
    % -------------------------------------------------
    polarR0Data{ii}.X0s          = X0s_ii;
    polarR0Data{ii}.Xfs          = Xfs_ii;
    polarR0Data{ii}.latLons      = latLons_ii;
    polarR0Data{ii}.endTimes     = endTimes_ii;

    
end
% -------------------------------------------------
% Saving data
% -------------------------------------------------
%%% Filename
filename_polarTraj = fullfile(polarTrajDataPath, sprintf('PolarXs_%s_%1.0fmps_%1.0fkm_%1.0fv0s.txt',...
    secondary.name(1:3),dvLp_mps,r0GridSpacing_km,n_v0s_per_r0));
    
%%% Opening files and writing header
f_polarTraj = fopen(filename_polarTraj, 'wt');
fprintf(f_polarTraj,'t0/tf,x,y,z,dx,dy,dz\n');  % header

%%% Writing data
for kk = 1:n_r0s
    for jj = 1:n_v0s_per_r0
        %%% If this trajectory impacted, write impact data
        if isnan(polarR0Data{kk}.endTimes(jj)) == 0
            %%% X0 (t0 = 0)
            fprintf(f_polarTraj,'0,%1.16f,%1.16f,%1.16f,%1.16f,%1.16f,%1.16f\n',polarR0Data{kk}.X0s(jj,:));
            %%% Xf (tf = tf);
            fprintf(f_polarTraj,'%1.8f,%1.16f,%1.16f,%1.16f,%1.16f,%1.16f,%1.16f\n',polarR0Data{kk}.endTimes(jj),polarR0Data{kk}.Xfs(jj,:));
        end
    end

end

%%% Close file
fclose(f_polarTraj);

toc






