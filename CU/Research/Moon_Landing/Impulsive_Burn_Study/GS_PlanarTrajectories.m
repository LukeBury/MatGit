% ========================================================================
%%% Description
% ========================================================================
% 

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
testCaseOn    = 0; % Lowers the resolution

plotResults   = 1;
storeAllTrajs = 1;
% ========================================================================
%%% Setup
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
warning('This spacing method will result in different resolutions at different energies')
if isequal(secondary.name,'europa')
    %%% Energy levels to run in terms of L2-flyover-velocity
    % L2FlyoverVelocity_mps_vec = [50, 100, 150, 200, 250, 300, 350]; % m/s
    L2FlyoverVelocity_mps_vec = [50]; % m/s

%     %%% Spacing of initial positions within y-z-plane neck
%     r0GridSpacing_km = 50; % km
    %%% Total number of r0s
    n_r0s = 200;

    %%% Number of velocity orientations per position
    n_v0s_per_r0 = 90;
    
    
end
if isequal(secondary.name,'enceladus')
    %%% Energy levels to run in terms of L2-flyover-velocity
    % L2FlyoverVelocity_mps_vec = [13, 26, 39, 52, 65, 78, 91] % m/s
    L2FlyoverVelocity_mps_vec = [26]; % m/s
    
    %%% Total number of r0s
    n_r0s = 200;

    %%% Number of velocity orientations per position
    n_v0s_per_r0 = 90;
end


%%% Lowering resolution for test cases
if testCaseOn == 1
    if isequal(secondary.name,'europa')
%         r0GridSpacing_km = 200;
        %%% Total number of r0s
        n_r0s = 10;
        
        %%% Number of velocity orientations per position
        n_v0s_per_r0 = 5;
    end
    if isequal(secondary.name,'enceladus')
        %%% Total number of r0s
        n_r0s = 10;
        
        %%% Number of velocity orientations per position
        n_v0s_per_r0 = 5;
    end
end

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
options_ImpactEscape = odeset('Events',@event_ImpactEscape_CR3Bn,'RelTol',tol,'AbsTol',tol);

% ========================================================================
%%% Preparing for Grid Search - Looping through L2FlyoverVelocity_mps values
% ========================================================================
for L2FlyoverVelocity_mps = [L2FlyoverVelocity_mps_vec]
% -------------------------------------------------
% %% Initial JC of spacecraft
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

%%% Populating vector of initial positions
r0s =  [ones(n_r0s,1).*L123(LPoint,1), linspace(-y_neck_upper, y_neck_upper, n_r0s)', zeros(n_r0s,1)];

% -------------------------------------------------
% Velocity vectors at each point
% -------------------------------------------------
[vHats2] = vHatPlanarHemisphere(n_v0s_per_r0,'-x','xy');

%%% For reference, total number of trajectories to simulate
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
    prms      = struct();
    prms.u    = MR;
    prms.R2_n = R2_n;
    prms.L1x  = L123mat(1,1);
    prms.L2x  = L123mat(2,1);

    % -------------------------------------------------
    % Setting current initial position
    % -------------------------------------------------
    r0_i = r0s(ii,:);
    
    
    % -------------------------------------------------
    % Finding necessary V0 to maintain correct JC0
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
    [v0i_mag] = dJC_2_dv(dJC_forInitialVelocity);

    % -------------------------------------------------
    % Preallocating data storage for each velocity vector
    % -------------------------------------------------    
    %%% Initializing data matrices
    trajIDs_ii              = zeros(n_v0s_per_r0,1);
    X0s_ii                  = zeros(n_v0s_per_r0,6);
    impactLongitudes_deg_ii = zeros(n_v0s_per_r0,1);
    impactAngles_deg_ii     = zeros(n_v0s_per_r0,1);
    TOFs_ii                 = zeros(n_v0s_per_r0,1);

    allTrajsAtThisr0_ii      = [];

    % -------------------------------------------------
    % Looping through all velocities at the position
    % ------------------------------------------------- 
    for vi = 1:n_v0s_per_r0
        % ---------------------------------------
        % Preallocating temporary variables
        % ---------------------------------------
%         trajID          = [];
        impactLon_deg   = [];
        impactAngle_deg = [];
        endTime         = [];

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

        [time_n, X_BCR_n, time_eventImpact, X_eventImpact, index_eventImpact] = ode113(@Int_CR3Bn,...
            time0_n, X0_n, options_ImpactEscape, prms);     
        
        % ---------------------------------------
        % Determining end conditions
        % ---------------------------------------
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
        
%         if impact1_escape2_orbit3 == 2
%             if X_eventImpact(1) <= prms.L1x
%                 format long
%                 time_n(end)
%                 X0_n
%                 fprintf('ttt')
%             end
%         end
        
%         % ---------------------------------------
%         % Store each surface crossing
%         % ---------------------------------------
%         if isempty(time_eventImpact) == 0 % if there was at least one event
%             for event_i = 1:length(time_eventImpact)
%                 %%% Finding lat/lon of event site
%                 [lat_deg, lon_deg] = BCR2latlon(X_eventImpact(event_i,1:3), 'secondary', MR);
%                 eventLatLon = [lat_deg, lon_deg];
%             end
%         else
%             eventLatLon = [NaN, NaN];
%         end
%         
%         % --------------------------
%         % Maximum latitude during trajectory
%         % --------------------------
%         for kk = 1:size(X_BCR_n,1)
%             [lat_kk_deg, lon_kk_deg] = BCR2latlon(X_BCR_n(kk,1:3), 'secondary', MR);
%             if abs(lat_kk_deg) > maxTrajLat
%                 maxTrajLat = abs(lat_kk_deg);
%             end
%             if storeAllLatLons == 1
%                 allLatLons(kk,:) = [lat_kk_deg, lon_kk_deg];
%             end
%         end
%         
        % ---------------------------------------
        % Determining impact conditions
        % ---------------------------------------
        if impact1_escape2_orbit3 == 1
            % --------------------------
            % Impact Longitude
            % --------------------------
            %%% Creating SCR position
            rImpact_SCR = X_eventImpact(1,1:3) - [1-MR, 0, 0];

            %%% Finding lat/lon of impact site
            [impactLat_deg, impactLon_deg] = BCR2latlon(X_eventImpact(1,1:3), 'secondary', MR);         
            
            % --------------------------
            % Impact Angle
            % --------------------------
            %%% Creating SCR position vector
            rImpact_SCR_n = X_eventImpact(end,1:3) - [1-MR,0,0];
            rImpact_SCR_n = rImpact_SCR_n./norm(rImpact_SCR_n);

            %%% Velocity unit vector at impact
            vHatImpact_n = X_eventImpact(end,4:6)./norm(X_eventImpact(end,4:6));

            %%% Angle between velocity and surface
            [impactAngle_deg] = calcImpactAngle(rImpact_SCR_n,vHatImpact_n,'degrees');
            
            % --------------------------
            % TOF
            % --------------------------
            %%% Storing endTime
            endTime = time_n(end);
            
            
        elseif impact1_escape2_orbit3 == 2 % escape
            impactLon_deg   = NaN;
            impactAngle_deg = NaN;     
            endTime         = NaN;
        elseif impact1_escape2_orbit3 == 3 % still orbiting
            impactLon_deg   = NaN;
            impactAngle_deg = NaN; 
            endTime         = NaN;
        end
        
        % =======================================
        % Storing results of trajectory
        % ---------------------------------------
        X0s_ii(vi,:)                  = X0_n;
        trajIDs_ii(vi)                = trajID;
        impactLongitudes_deg_ii(vi,:) = impactLon_deg;
        impactAngles_deg_ii(vi)       = impactAngle_deg;
        TOFs_ii(vi)                   = endTime;
        
%         eventLatLons_ii(vi,:) = eventLatLon;
%         maxTrajLats_ii(vi)    = maxTrajLat;

        if storeAllTrajs == 1
            allTrajsAtThisr0_ii = [allTrajsAtThisr0_ii; X_BCR_n(:,1:3); NaN(1,3)];
            
% %             allLatLons_ii{vi} = [allLatLons; NaN, NaN];
%             allTrajsAtThisr0_ii(((vi-1)*n_t+1):(vi*n_t),:) = X_BCR_n(:,1:3);
%         
%             %%% Cleaning up allTrajs_ii (getting rid of nans and zeros
%             nonNanIndices = ~isnan(allTrajsAtThisr0_ii(:,1));
%             allTrajsAtThisr0_ii = allTrajsAtThisr0_ii(nonNanIndices,:);
%             nonPlanarIndices = find(allTrajsAtThisr0_ii(:,1) ~= 0);
%             allTrajsAtThisr0_ii = allTrajsAtThisr0_ii(nonPlanarIndices,:);
        end
        
    end % vi = 1:n_v0s_per_r0
%     
    % -------------------------------------------------
    % Storing data from all azimuths/elevations
    % -------------------------------------------------
    r0Data{ii}.X0s                  = X0s_ii;
    r0Data{ii}.trajIDs              = trajIDs_ii;
    r0Data{ii}.impactLongitudes_deg = impactLongitudes_deg_ii;
    r0Data{ii}.impactAngles_deg     = impactAngles_deg_ii;
    r0Data{ii}.TOFs                 = TOFs_ii;
%     r0Data{ii}.eventLatLons     = eventLatLons_ii;
%     r0Data{ii}.maxTrajLats      = maxTrajLats_ii;

    if storeAllTrajs == 1
        r0Data{ii}.allTrajs       = allTrajsAtThisr0_ii;
    end
end % parfor ii = 1:n_r0s
% % -------------------------------------------------
% % Rearranging certain data
% % -------------------------------------------------
% eventLatLons_full = [];
% maxTrajLats_full  = [];

trajIDs              = [];
X0s                  = [];
impactLongitudes_deg = [];
impactAngles_deg     = [];
TOFs                 = [];
if storeAllTrajs == 1
    allTrajs             = [];
end

for kk = 1:length(r0Data)
%     eventLatLons_full = [eventLatLons_full; r0Data{kk}.eventLatLons];
%     maxTrajLats_full = [maxTrajLats_full; r0Data{kk}.maxTrajLats];
    X0s                  = [X0s; r0Data{kk}.X0s];
    trajIDs              = [trajIDs; r0Data{kk}.trajIDs];
    impactLongitudes_deg = [impactLongitudes_deg; r0Data{kk}.impactLongitudes_deg];
    impactAngles_deg     = [impactAngles_deg; r0Data{kk}.impactAngles_deg];
    TOFs                 = [TOFs; r0Data{kk}.TOFs];
    
    
    if storeAllTrajs == 1
        allTrajs = [allTrajs; r0Data{kk}.allTrajs];
    end
end
% 
% % -------------------------------------------------
% % Writing data to CSV
% % -------------------------------------------------
% %%% File Names
% if testCaseOn == 0
%     %%% Latitdue/Longitude data
%     filename_LatStudy = fullfile(savepath, sprintf('LatStudy_%s.iGS_%s_%1.0fmps_%1.0fkm_%1.0fv0s.txt',...
%         computerTag,secondary.name(1:3),L2FlyoverVelocity_mps,r0GridSpacing_km,n_v0s_per_r0));
%     
%     %%% If storing all lat/lons
%     if storeAllLatLons == 1
%         filename_allLatLons = fullfile(savepath, sprintf('LatStudy_allLatLons_%s.iGS_%s_%1.0fmps_%1.0fkm_%1.0fv0s.txt',...
%         computerTag,secondary.name(1:3),L2FlyoverVelocity_mps,r0GridSpacing_km,n_v0s_per_r0));
%     end
%     
% elseif testCaseOn == 1
%     %%% Latitude/Longitude data
%     filename_LatStudy     = fullfile(savepath,'LatStudy_testFile.txt');
%     
%     %%% If storing all lat/lons
%     if storeAllLatLons == 1
%         filename_allLatLons     = fullfile(savepath,'LatStudy_allLatLons_testFile.txt');
%     end
% end
% 
% %%% Opening files and writing header
% f_LatStudy = fopen(filename_LatStudy, 'wt');
% runTime_min = toc(ticWhole)/60;
% fprintf(f_LatStudy,['y0_n,z0_n,maxLatitude(deg),latitude(deg),longitude(deg)...',...
%     sprintf('runtime=%1.1fmin,',runTime_min),...
%     sprintf('MaxImpactLat=%1.3f,MaxTrajLat=%1.3f',max(abs(eventLatLons_full(:,1))),max(abs(maxTrajLats_full))),...
%     '\n']);  % header
% 
% if storeAllLatLons == 1
%     %%% Opening files and writing header
%     f_allLatLons = fopen(filename_allLatLons, 'wt');
%     fprintf(f_allLatLons,'latitude(deg),longitude(deg)\n');  % header
% end
% 
% %%% Writing data (and finding maximum latitudes)
% for r0_k = 1:n_r0s
%     %%% If there were events from this r0
%     if isempty(r0Data{r0_k}) == 0
%         for v0_k = 1:n_v0s_per_r0
%             fprintf(f_LatStudy,'%1.16f,%1.16f,%1.16f,%1.16f,%1.16f,%2.3f,%2.3f,%2.3f\n',...
%                 r0Data{r0_k}.X0s(v0_k,2),...
%                 r0Data{r0_k}.X0s(v0_k,3),...
%                 r0Data{r0_k}.X0s(v0_k,4),...
%                 r0Data{r0_k}.X0s(v0_k,5),...
%                 r0Data{r0_k}.X0s(v0_k,6),...
%                 r0Data{r0_k}.maxTrajLats(v0_k),...
%                 r0Data{r0_k}.eventLatLons(v0_k,1),...
%                 r0Data{r0_k}.eventLatLons(v0_k,2));
%             
%         end
%         
%         if storeAllLatLons == 1
%             for LatLon_k = 1:size(r0Data{r0_k}.allLatLons(:,1))
%                 fprintf(f_allLatLons,'%2.1f,%2.1f\n',...
%                 r0Data{r0_k}.allLatLons(LatLon_k,1),...
%                 r0Data{r0_k}.allLatLons(LatLon_k,2));
%             end
%         end
%     end
% 
% end
% 
% %%% Close file
% fclose(f_LatStudy);

end % L2FlyoverVelocity_mps = [L2FlyoverVelocity_mps_vec]





% ========================================================================
%%% Plotting Results
% ========================================================================
if (plotResults == 1) && isequal(computerTag,'M')
    figure; hold all
    plot(impactLongitudes_deg,impactAngles_deg,'b.','markersize',10)
    PlotBoi2('Impact Longitude, $^\circ$','Impact Angle, $^\circ$',18,'LaTex')
    
    
    if storeAllTrajs == 1
        figure; hold all
        axis equal
        plotCR3BP_Neck(secondary,L123,JC_scInitial,600,200,colors.std.black,1.5)
        prms.R2_n = secondary.R_n;
        plotBodyTexture3(secondary.R_n,[1-secondary.MR,0,0],secondary.img)
        plot(allTrajs(:,1),allTrajs(:,2),'b','linewidth',1)
        view(0,90)
        PlotBoi2('$X_n$','$Y_n$',20,'LaTex')
    end
end % (plotResults == 1) && isequal(computerTag,'M')
% 
% 


fprintf('Finished: %1.1f Minutes\n',toc(ticWhole)/60)









