clear
clc
close all

ticWhole = tic;
% ========================================================================
%%% Run Switches
% ========================================================================
%%% Testing?
testCaseOn             = 1;

%%% Plot switches
plot_initialConditions = 1;
plot_fullSystemContour = 1;
plot_sectionColors     = 1;

%%% Computer?
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

%%% How fast the SC would be traveling over the Lagrange point
% dvLp_mps = 200; % Meters per second - Europa
% dvLp_mps = 50; % Meters per second - Enceladus
dvLp_mps = 200; % Meters per second

%%% Spacing of initial positions within 3D neck
% r0GridSpacing_km = 100; % km - Europa
% r0GridSpacing_km = 10;  % km - Enceladus
r0GridSpacing_km = 50; % km

%%% Spacing between azimuths and elevations of v0s per r0
n_v0s_per_r0_target = 145;
% n_target = 1297; % ~ old  5 deg spacing
% n_target = 352;  % ~ old 10 deg spacing
% n_target = 145;  % ~ old 15 deg spacing - should be used
% n_target = 37;   % ~ old 30 deg spacing - decent
% n_target = 17;   % ~ old 45 deg spacing
% n_target = 5;    % ~ old 90 deg spacing

if testCaseOn == 1
    r0GridSpacing_km = 2000;
    n_v0s_per_r0_target = 350;
end

%%% Selecting time vector
t_i = 0; % sec
t_f = 6*pi; % Long bc events are watching for impact or escape
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

%%% Total number of r0s
n_r0s = size(r0s,1);

% -------------------------------------------------
% Creating Bins for impact angle and neck sections
% -------------------------------------------------
%%% Bins for impact angles
% % Making last bin too big so everything fits in bins(1:end-1)
bins_impactAngles = [0, 5, 20, 40, 60, 91];

%%% Bins for neck sections
smallNeck = smallerValue(z_neck_upper, y_neck_upper);
largeNeck = largerValue(z_neck_upper, y_neck_upper);
bins_neckSections = [linspace(0,smallNeck,binCount_NeckSections), largeNeck];

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
tol = 1e-10;

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
    latLons_ii          = zeros(n_v0s_per_r0,2);
    trajs_ii            = zeros(n_v0s_per_r0,3);
    impactSpeeds_ii     = zeros(n_v0s_per_r0,1);
    impactAngles_ii     = zeros(n_v0s_per_r0,1);
    
    bin_impactAngles_ii = zeros(n_v0s_per_r0,1);
    bin_neckSections_ii = zeros(n_v0s_per_r0,1);
    
%     landingTrajX0s_ii   = zeros(n_v0s_per_r0,6);
    
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
        endCase = []; % 1-impact, 2-escape, 3-still in orbit
        if isempty(X_eventImpact) == 0
            if abs(norm(X_eventImpact(1:3)-[1-MR,0,0])-R2_n) < 1e-9 % impact
                endCase = 1;
            elseif X_eventImpact(1) <= L123mat(1,1) || X_eventImpact(1) >= L123mat(2,1) % escape
                endCase = 2;
            end
        elseif isempty(X_eventImpact) == 1 % No event yet, still in orbit
            endCase = 3;
        end
        
        
        % ---------------------------------------
        % Preallocating temporary variables
        % ---------------------------------------
        impactLatLon = [];
        impactSpeed  = [];
        impactAngle  = [];
%         landingTrajX0  = [];

        % ---------------------------------------
        % Determining lat/lon, speed, and angle of impact
        % ---------------------------------------

        if endCase == 1
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

        elseif endCase == 2 % escape
            impactLatLon = [NaN, NaN];
            impactSpeed  = NaN;
            impactAngle  = NaN;     
        elseif endCase == 3 % orbit
            impactLatLon = [NaN, NaN];
            impactSpeed  = NaN;
            impactAngle  = NaN; 
        end
            
        % =======================================
        % Storing results of trajectory
        % ---------------------------------------
        X0s_ii(vi,:)        = X0_n;
        latLons_ii(vi,:)    = impactLatLon;
        impactSpeeds_ii(vi) = impactSpeed;
        impactAngles_ii(vi) = impactAngle;

        % ---------------------------------------
        % Assigning bins
        % ---------------------------------------
        %%% Bin - Impact Angle
        bin_impactAngles_ii(vi) = discretize(impactAngle, bins_impactAngles);

        %%% Bin - Neck Section
        % Calculate radial distance from r0 to LPoint
        r0LPointDist = norm(X0_n(1:3) - L123mat(Lpoint,:));

        bin_neckSections_ii(vi) = discretize(r0LPointDist, bins_neckSections);

%         % ---------------------------------------
%         % Storing Low-landing-angle trajectories
%         % ---------------------------------------
%         if isempty(X_eventImpact) == 0 && bin_impactAngles_ii(vi) == 1
%             landingTrajX0 = X0_n;
%             landingTrajX0s_ii = [landingTrajX0s_ii ;landingTrajX0];
%         elseif isempty(X_eventImpact) == 1
%             landingTrajX0 = NaN;
%         end

    end
    
    % -------------------------------------------------
    % Storing data from all azimuths/elevations
    % -------------------------------------------------
    r0Data{ii}.X0s          = X0s_ii;
    r0Data{ii}.latLons      = latLons_ii;
    r0Data{ii}.impactSpeeds = impactSpeeds_ii;
    r0Data{ii}.impactAngles = impactAngles_ii;
    
%     r0Data{ii}.landingTrajs = landingTrajX0s_ii;
    
    r0Data{ii}.bin_neckSections = bin_neckSections_ii;
    r0Data{ii}.bin_impactAngles = bin_impactAngles_ii;
    
end


% -------------------------------------------------
% Writing data to CSV
% -------------------------------------------------
%%% File Names
if testCaseOn == 0
    %%% All impact data
    filename_allTraj = fullfile(savepath, sprintf('%s.iGS_%sL%1.0f_%1.0fmps_%1.0fkm_%1.0fv0s_data.txt',...
        computerTag,secondary.name(1:3),Lpoint,dvLp_mps,r0GridSpacing_km,n_v0s_per_r0));
    %%% Low-impact-angle data
    filename_landingTraj = fullfile(savepath, sprintf('%s.iGS_%sL%1.0f_%1.0fmps_%1.0fkm_%1.0fv0s_land.txt',...
        computerTag,secondary.name(1:3),Lpoint,dvLp_mps,r0GridSpacing_km,n_v0s_per_r0));
    %%% Log file
    filename_runData = fullfile(savepath, sprintf('%sLog.iGS_%sL%1.0f_%1.0fmps_%1.0fkm_%1.0fv0s_log.txt',...
        computerTag,secondary.name(1:3),Lpoint,dvLp_mps,r0GridSpacing_km,n_v0s_per_r0));
elseif testCaseOn == 1
    %%% All impact data
    filename_allTraj     = fullfile(savepath,'testFile_data.txt');
    %%% Low-impact-angle data
    filename_landingTraj = fullfile(savepath,'testFile_land.txt');
    %%% Log file
    filename_runData     = fullfile(savepath,'testFile_log.txt');
end

%%% Opening files and writing header
f_allTraj = fopen(filename_allTraj, 'wt');
fprintf(f_allTraj,'bin_impactAngle,bin_neckSection,latitude,longitude,impactAngle\n');  % header

f_landingTraj = fopen(filename_landingTraj, 'wt');
fprintf(f_landingTraj,'x0_n,y0_n,z0_n,dx0_n,dy0_n,dz0_n,bin_neckSection,latitude,longitude\n');  % header

%%% Writing data
lowImpactAngleCounter = 0;
for kk = 1:n_r0s
    for jj = 1:n_v0s_per_r0
        %%% If this trajectory impacted, write impact data
        if isnan(r0Data{kk}.bin_impactAngles(jj)) == 0
            fprintf(f_allTraj,'%1d,%1d,%2.1f,%2.1f,%2.1f\n',r0Data{kk}.bin_impactAngles(jj),...
                r0Data{kk}.bin_neckSections(jj),r0Data{kk}.latLons(jj,1),r0Data{kk}.latLons(jj,2),r0Data{kk}.impactAngles(jj));
            
            %%% If this was a low-impact-angle trajectory
            if r0Data{kk}.bin_impactAngles(jj) == 1
                lowImpactAngleCounter = lowImpactAngleCounter + 1;
                fprintf(f_landingTraj,'%1.11f, %1.11f, %1.11f, %1.10f, %1.10f, %1.10f, %1d, %2.1f, %2.1f\n',r0Data{kk}.X0s(jj,:),...
                    r0Data{kk}.bin_neckSections(jj),r0Data{kk}.latLons(jj,1),r0Data{kk}.latLons(jj,2));
            end
            
        end
    end
%     if isempty(r0Data{kk}.landingTrajs) == 0
%         for rr = 1:size(r0Data{kk}.landingTrajs,1)
%             if isnan(r0Data{kk}.landingTrajs(rr,1)) == 1
%                 fprintf(f_landingTraj,'NaN, NaN, NaN, NaN, NaN, NaN, NaN\n');
%             else
%                 fprintf(f_landingTraj,'%1.6f, %1.8f, %1.8f, %1.8f, %1.8f, %1.8f, %1.8f\n',r0Data{kk}.X0s(rr,:));
%             end
%             
%         end
%     end
end
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
% Creating file of Fortuna info
% -------------------------------------------------
%%% Writing log file
f_runData = fopen(filename_runData, 'wt');
fprintf(f_runData,'Run Time (sec):\t\t%1.1f\n',finalToc);
fprintf(f_runData,'Run Time (min):\t\t%1.1f\n',finalToc/60);
fprintf(f_runData,'number of r0s:\t\t%1.0f\n',n_r0s);
fprintf(f_runData,'v0s per r0s:\t\t%1.0f\n',n_v0s_per_r0);
fprintf(f_runData,'total trajs:\t\t%1.0f\n',n_traj);
% fprintf(f_runData,'bins_impactAngles:\t%1.0f,%1.0f,%1.0f,%1.0f,%1.0f,%1.0f\n',bins_impactAngles);
% fprintf(f_runData,'bins_neckSections:\t%1.0f,%1.0f,%1.0f,%1.0f,%1.0f,\n',bins_neckSections);
% fprintf(f_runData,'JC_scInitial:\t\t%1.15f,,,,,\n',JC_scInitial);
% fprintf(f_runData,'Lpoint:\t\t\t\t%1.0d,,,,,\n',Lpoint);
% fprintf(f_runData,'y_neck_upper:\t\t%1.15f,,,,,\n',y_neck_upper);
% fprintf(f_runData,'z_neck_upper:\t\t%1.15f,,,,,\n',z_neck_upper);

%%% Close file
fclose(f_runData);

finalToc = toc(ticWhole);


















