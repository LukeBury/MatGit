clear
clc
close all
mbinPath = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/mbin';
addpath(genpath(mbinPath))
ticWhole = tic;
% ========================================================================
%%% Importing Data
% ========================================================================
%%% General data on solar system bodies
bodies = getBodyData(mbinPath);

%%% Color options/schemes
colors = get_colors();

% ========================================================================
%%% Run Switches
% ========================================================================
testCaseOn = 0; % Lowers the resolution

generateData = 1; % Run the grid search
plotResults  = 0; % Look at data from grid search
    plotAllLatLons = 1;
storeAllLatLons = 1; % Store lat/lon at each time step for each trajectory
% ========================================================================
%%% Setup
% ========================================================================
% -------------------------------------------------
%%% Paths and Data
% -------------------------------------------------
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
else 
    warning('This computer will explode in 5 seconds')
end

%%% Add the function paths to matlab
addpath(genpath(mbinPath))
addpath(genpath(moonFuncsPath))

%%% General data on solar system bodies
bodies = getBodyData(mbinPath);

%%% Color options/schemes
colors = get_colors();

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
%%% Grid Search Settings
% -------------------------------------------------
%%% Energy levels to run in terms of L2-flyover-velocity
L2FlyoverVelocity_mps_vec = [50,100]; % m/s

if isequal(secondary.name,'europa')
    %%% Spacing of initial positions within y-z-plane neck
    r0GridSpacing_km = 50; % km
    
    %%% Number of v0s per r0
    n_v0s_per_r0_target = 145;
end


%%% Lowering resolution for test cases
if testCaseOn == 1
    if isequal(secondary.name,'europa')
        r0GridSpacing_km = 500;
        n_v0s_per_r0_target = 5;
    end
end

% -------------------------------------------------
%%% Plot settings
% -------------------------------------------------
LatStudyEventFile = [savepath,'/LatStudy_testFile.txt'];
LatStudyAllLatLonsFile = [savepath,'/LatStudy_allLatLons_testFile.txt'];

% -------------------------------------------------
%%% Integration settings
% -------------------------------------------------
%%% Time vector
t0 = 0;
tf = 4*pi;
n_t = 1000;
time0_n = linspace(t0, tf, n_t);

%%% Choosing ode tolerance
tol = 1e-13;

%%% Setting integrator options
options_ImpactNoStop_EscapeStop = odeset('Events',@event_ImpactNoStop_EscapeStop_CR3Bn,'RelTol',tol,'AbsTol',tol);

% ========================================================================
%%% Preparing for Grid Search - Looping through dvLp values
% ========================================================================
if generateData == 1
for L2FlyoverVelocity_mps = [L2FlyoverVelocity_mps_vec]
% -------------------------------------------------
%%% Initial JC of spacecraft
% -------------------------------------------------
[JC_scInitial] = L2FlyoverVelocity_2_JC(L2FlyoverVelocity_mps, MR, L123(LPoint,:), vNorm);

% -------------------------------------------------
% Finding upper-y-value of neck at L-point
% -------------------------------------------------
%%% Set values
c = JC_scInitial;
u = secondary.MR;
x = L123(LPoint,1);
z = 0;

%%% JC function equal to zero
r1 = @(x,y,z) sqrt((x+u)^2+y.^2+z^2);
r2 = @(x,y,z) sqrt((x-1+u)^2+y.^2+z^2);
f = @(y) x^2 + y.^2 - c + 2*(1-u)./r1(x,y,z) + 2*u./r2(x,y,z);

%%% Find the root of this function in the appropriate range
y_neckRange = 15*secondary.R_n;
z_neckRange = 15*secondary.R_n;

y_neck_upper = fzero(f,[0 y_neckRange]);

%%% Clear variables
clear c u x f z

% -------------------------------------------------
% Finding upper-z-value of neck at L-point
% -------------------------------------------------
%%% Set values
c = JC_scInitial;
u = secondary.MR;
x = L123(LPoint,1);
y = 0;

%%% JC function equal to zero
r1 = @(x,y,z) sqrt((x+u)^2+y.^2+z^2);
r2 = @(x,y,z) sqrt((x-1+u)^2+y.^2+z^2);

f = @(z) x^2 + y.^2 - c + 2*(1-u)./r1(x,y,z) + 2*u./r2(x,y,z);

%%% Find the root of this function in the appropriate range
z_neck_upper = fzero(f,[0 z_neckRange]);

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
JCs_yz_LPoint = zeros(size(Y_yz));
for yk = 1:size(Y_yz,1)
    for zk = 1:size(Y_yz,2)
        %%% Zero-Velocity Curve
        zv = JacobiConstantCalculator(secondary.MR,[L123(LPoint,1), Y_yz(yk,zk), Z_yz(yk,zk)] ,[0, 0, 0]);
        
        JCs_yz_LPoint(yk,zk) = zv;
    end
end

%%% Get points of y-z contour in 3D space
[ yzNeckContourPoints ] = getContourPoints( Y_yz, Z_yz, JCs_yz_LPoint, JC_scInitial );

% -------------------------------------------------
% Create grid of starting locations to sample r0s from
% -------------------------------------------------
n_r0Grid_y = 2*(2*y_neck_upper) / (r0GridSpacing_km/rNorm);
n_r0Grid_z = 2*(2*z_neck_upper) / (r0GridSpacing_km/rNorm);
ys = linspace(-2*y_neck_upper, 2*y_neck_upper, n_r0Grid_y);
zs = linspace(-2*z_neck_upper, 2*z_neck_upper, n_r0Grid_z);

[Y_yz_r0s,Z_yz_r0s] = meshgrid(ys,zs);

%%% Search grid for points lying within y-z contour and keep those as r0s
%%% *** Only keeping those with z>=0 so the results can later be flipped
r0s = [];
for yk = 1:size(Y_yz_r0s,1)
    for zk = 1:size(Z_yz_r0s,2)
        if inpolygon(Y_yz_r0s(yk,zk), Z_yz_r0s(yk,zk),yzNeckContourPoints(1,:),yzNeckContourPoints(2,:)) == 1
            if Z_yz_r0s(yk,zk) >= 0 % dynamics are symmetric about z=0, so only need to check one side of that
                r0s = [r0s; L123(2,1), Y_yz_r0s(yk,zk), Z_yz_r0s(yk,zk)];
            end
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

% ========================================================================
%%% Loop through conditions in parallel and store results
% ========================================================================
%%% Create data structure
r0Data = {};

%%% Loop through positions in parallel
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
    JC_initialPos = [];
    JC_initialPos = JacobiConstantCalculator(MR,r0_i,[0,0,0]);
    
    %%% Starting velocity is found from difference between s/c JC (JC_scDesired) and the
    %%% JC of the stationary starting position (JC_initialPos)
    dJC_forInitialVelocity = JC_initialPos - JC_scInitial;
    
    %%% Setting variables so they're not considered 'temporary variables'
    v0i_mag = [];
    
    %%% Find necessary velocity magnitude
%     if dJC_forInitialVelocity < 0
%         warning('Spacecraft starting in a forbidden zone')
%     elseif dJC_forInitialVelocity >= 0
%         v0i_mag = sqrt(abs(dJC_forInitialVelocity));
%     end
    [v0i_mag] = dJC_2_dv(dJC_forInitialVelocity);
    
    % -------------------------------------------------
    % Preallocating data storage for each velocity vector
    % -------------------------------------------------    
    %%% Initializing data matrices
    X0s_ii              = zeros(n_v0s_per_r0,6);
    eventLatLons_ii     = zeros(n_v0s_per_r0,2);
    maxLats_ii          = zeros(n_v0s_per_r0,1);
    if storeAllLatLons == 1
        allLatLons_ii       = zeros(n_v0s_per_r0*n_t, 2);
    end
    
    % -------------------------------------------------
    % Looping through all velocities at the position
    % ------------------------------------------------- 
    for vi = 1:n_v0s_per_r0
        % ---------------------------------------
        % Preallocating temporary variables
        % ---------------------------------------
        trajID       = [];
        eventLatLon  = [];
        maxLat       = 0;
        if storeAllLatLons == 1
            allLatLons = NaN(n_t,2);
        end
        
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

        [time_n, X_BCR_n, time_eventImpact, X_eventImpact, index_eventImpact] = ode113(@Int_CR3Bn,...
            time0_n, X0_n, options_ImpactNoStop_EscapeStop, prms);
        
        % ---------------------------------------
        % Store each surface crossing
        % ---------------------------------------
        if isempty(time_eventImpact) == 0 % if there was at least one event
            for event_i = 1:length(time_eventImpact)
                %%% Finding lat/lon of event site
                [lat_deg, lon_deg] = BCR2latlon(X_eventImpact(event_i,1:3), 'secondary', MR);
                eventLatLon = [lat_deg, lon_deg];
            end
        else
            eventLatLon = [NaN, NaN];
        end
        
        % --------------------------
        % Maximum latitude during trajectory
        % --------------------------
        for kk = 1:size(X_BCR_n,1)
            [lat_kk_deg, lon_kk_deg] = BCR2latlon(X_BCR_n(kk,1:3), 'secondary', MR);
            if abs(lat_kk_deg) > maxLat
                maxLat = abs(lat_kk_deg);
            end
            if storeAllLatLons == 1
                allLatLons(kk,:) = [lat_kk_deg, lon_kk_deg];
            end
        end
            
        
        % =======================================
        % Storing results of trajectory
        % ---------------------------------------
        X0s_ii(vi,:)          = X0_n;
        eventLatLons_ii(vi,:) = eventLatLon;
        maxLats_ii(vi)        = maxLat;
        if storeAllLatLons == 1
%             allLatLons_ii{vi} = [allLatLons; NaN, NaN];
            allLatLons_ii(((vi-1)*n_t+1):(vi*n_t),:) = allLatLons;
        end
        
        %%% Cleaning up allLatLons_ii (getting rid of nans and zeros
        nonNanIndices = ~isnan(allLatLons_ii(:,1));
        allLatLons_ii = allLatLons_ii(nonNanIndices,:);
        nonPlanarIndices = find(allLatLons_ii(:,1) ~= 0);
        allLatLons_ii = allLatLons_ii(nonPlanarIndices,:);
        
    end % vi = 1:n_v0s_per_r0
    
    % -------------------------------------------------
    % Storing data from all azimuths/elevations
    % -------------------------------------------------
    r0Data{ii}.X0s              = X0s_ii;
    r0Data{ii}.latLons          = eventLatLons_ii;
    r0Data{ii}.maxLats          = maxLats_ii;
    r0Data{ii}.allLatLons       = allLatLons_ii;
    
end % parfor ii = 1:n_r0s

% -------------------------------------------------
% Writing data to CSV
% -------------------------------------------------
%%% File Names
if testCaseOn == 0
    %%% Latitdue/Longitude data
    filename_LatStudy = fullfile(savepath, sprintf('LatStudy_%s.iGS_%s_%1.0fmps_%1.0fkm_%1.0fv0s.txt',...
        computerTag,secondary.name(1:3),dvLp_mps,r0GridSpacing_km,n_v0s_per_r0));
    
    %%% If storing all lat/lons
    if storeAllLatLons == 1
        filename_allLatLons = fullfile(savepath, sprintf('LatStudy_allLatLons_%s.iGS_%s_%1.0fmps_%1.0fkm_%1.0fv0s.txt',...
        computerTag,secondary.name(1:3),dvLp_mps,r0GridSpacing_km,n_v0s_per_r0));
    end
    
elseif testCaseOn == 1
    %%% Latitude/Longitude data
    filename_LatStudy     = fullfile(savepath,'LatStudy_testFile.txt');
    
    %%% If storing all lat/lons
    if storeAllLatLons == 1
        filename_allLatLons     = fullfile(savepath,'LatStudy_allLatLons_testFile.txt');
    end
end

%%% Opening files and writing header
f_LatStudy = fopen(filename_LatStudy, 'wt');
runTime_min = toc(ticWhole)/60;
fprintf(f_LatStudy,['y0_n,z0_n,maxLatitude(deg),latitude(deg),longitude(deg),',sprintf('runtime=%1.1fmin',runTime_min),'\n']);  % header

if storeAllLatLons == 1
    %%% Opening files and writing header
    f_allLatLons = fopen(filename_allLatLons, 'wt');
    fprintf(f_allLatLons,'latitude(deg),longitude(deg)\n');  % header
end

%%% Writing data (and finding maximum latitudes)
for r0_k = 1:n_r0s
    %%% If there were events from this r0
    if isempty(r0Data{r0_k}) == 0
        for v0_k = 1:n_v0s_per_r0
            fprintf(f_LatStudy,'%1.16f,%1.16f,%2.3f,%2.3f,%2.3f\n',...
                r0Data{r0_k}.X0s(v0_k,2),...
                r0Data{r0_k}.X0s(v0_k,3),...
                r0Data{r0_k}.maxLats(v0_k),...
                r0Data{r0_k}.latLons(v0_k,1),...
                r0Data{r0_k}.latLons(v0_k,2));
            
        end
        
        if storeAllLatLons == 1
            for LatLon_k = 1:size(r0Data{r0_k}.allLatLons(:,1))
                fprintf(f_allLatLons,'%2.1f,%2.1f\n',...
                r0Data{r0_k}.allLatLons(LatLon_k,1),...
                r0Data{r0_k}.allLatLons(LatLon_k,2));
            end
        end
    end

end

%%% Close file
fclose(f_LatStudy);

end % dvLp_mps = [dvLp_mps_vec]

end % if generateData == 1




% ========================================================================
%%% Plotting Results
% ========================================================================
if (plotResults == 1) && isequal(computerTag,'M')
% -------------------------------------------------
% Load data into matrix
% -------------------------------------------------
%%% Load data
LatStudyData = dlmread(LatStudyEventFile,',',1,0);

%%% set column specifiers
c_y0_n      = 1;
c_z0_n      = 2;
c_maxLat    = 3;
c_eventLat  = 4;
c_eventLon  = 5;

% -------------------------------------------------
% Plotting event map
% -------------------------------------------------
figure; hold all
plot(LatStudyData(:,c_eventLon),LatStudyData(:,c_eventLat),'.','markersize',8,'color',colors.std.black)
xlim([-180, 180])
ylim([-90, 90])

figure; hold all
plot(LatStudyData(:,c_maxLat),'.','markersize',8,'color',colors.std.black)

% -------------------------------------------------
% Plotting event map
% -------------------------------------------------
if plotAllLatLons == 1
%%% Load data
LatStudyAllLonLonsData = dlmread(LatStudyAllLatLonsFile,',',1,0);

figure; hold all
plot(LatStudyAllLonLonsData(:,2),LatStudyAllLonLonsData(:,1),'.','markersize',8,'color',colors.std.blue)


end

end % (plotResults == 1) && isequal(computerTag,'M')














