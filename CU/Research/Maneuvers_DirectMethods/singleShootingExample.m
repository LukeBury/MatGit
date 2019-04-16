clear
clc
close all
mbinPath = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/mbin';
moonFuncsPath = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/CU/Research/Moon_Landing/Moon_Landing_funcs';
addpath(genpath(mbinPath))
addpath(genpath(moonFuncsPath))

% ========================================================================
%%% Importing Data
% ========================================================================
%%% General data on solar system bodies
bodies = getBodyData(mbinPath);

%%% Color options/schemes
colors = get_colors();

% ========================================================================
%%% Switches
% ========================================================================
%%% The following are set up to be mutually exclusive (1 = on, 0 = off)
run_r0error = 0; % Perturb & correct only initial position
run_v0error = 0; % Perturb & correct only initial velocity
run_X0error = 1; % Perturb & correct full initial state

% ========================================================================
%%% Setup
% ========================================================================
% ------------------------------------
%%% Setup - system
% ------------------------------------
%%% bodies
primary = bodies.jupiter;
secondary = bodies.europa;

%%% Normalizing constants
rNorm = secondary.a;         % n <-> km
tNorm = 1/secondary.meanMot; % n <-> sec
vNorm = rNorm / tNorm;       % n <-> km/sec

% ------------------------------------
%%% Setup - options
% ------------------------------------
%%% Weight the dX0 correction
dX0_scale = 1;

%%% Number of iterations
n_iterations = 30;

%%% Set introduced error to X0
introducedPositionError = [0; 1; 1].*(2.0)./rNorm;
introducedVelocityError = [0; 1; 1].*(1e-9)./vNorm;
introducedStateError = [[0; 1; 1].*(1e-2)./rNorm; [0; 1; 1].*(1e-7)./vNorm];

% ------------------------------------
%%% Setup - initial integration options
% ------------------------------------
%%% Time
t_i = 0;
t_f = 4*pi;
n_dt = 10000;
time0_n = linspace(t_i,t_f,n_dt); % normalized time

%%% Choosing ode tolerance
tol = 1e-13;

%%% Setting integrator options
options = odeset('Event',@event_Impact_CR3Bn, 'RelTol',tol,'AbsTol',tol);

%%% Setting necessary parameters for integration
prms.u = secondary.MR;
prms.R2_n = secondary.R_n;

% ------------------------------------
%%% Generate 'truth' trajectory
% ------------------------------------
%%% Set an initial truth state - this initial condition set impacts Europa
%%% at (LatLon [11.713313624820113, -17.884817025059917]).
X0_X_truth = [1.0204617015266166,-0.0019850725823437,0.0000000000000000,-0.0010976413130697,-0.0030157447223195,0.0098771743854318]';

%%% Integrate truth
[time_n, X_BCR_n_truth] = ode113(@Int_CR3Bn, time0_n, X0_X_truth, options, prms);

%%% Finding truth impact site of trajectory
[targetLat_truth, targetLon_truth] = BCR2latlon(X_BCR_n_truth(end,1:3), 'secondary', secondary.MR);

%%% Hard-coding target state
rTarget     = X_BCR_n_truth(end,1:3)';
vTarget     = X_BCR_n_truth(end,4:6)';
stateTarget = X_BCR_n_truth(end,1:6)';

% ------------------------------------
%%% Add error to the truth state and prepare for shooting
% ------------------------------------
%%% Introduce an error to the initial state
if (run_r0error == 1) && (run_v0error == 0) && (run_X0error == 0)
    X0_X_error = X0_X_truth + [introducedPositionError; 0; 0; 0];
    
elseif (run_r0error == 0) && (run_v0error == 1) && (run_X0error == 0)
    X0_X_error = X0_X_truth + [0; 0; 0; introducedVelocityError];
    
elseif (run_r0error == 0) && (run_v0error == 0) && (run_X0error == 1)
    X0_X_error = X0_X_truth + introducedStateError;
end

%%% Initial stm is identity (t0 to t0)
stm0 = reshape(eye(6),36,1);

%%% Set full initial state that includes error
X0_error = [X0_X_error; stm0];

% ------------------------------------
%%% Initial integration
% ------------------------------------
%%% Integrate trajectory
[time_n, X_BCR_n_error] = ode113(@Int_CR3BnSTM, time0_n, X0_error, options, prms);

%%% Check if the provided error still led to an impact of the moon
%%% (basically making sure the end-error isn't absurdly large)
if (norm(X_BCR_n_error(end,1:3) - [1-secondary.MR,0,0]) - secondary.R_n) > (1/rNorm)
    warning('no impact')
    return
end

%%% Finding initial impact site of trajectory
[lat0, lon0] = BCR2latlon(X_BCR_n_error(end,1:3), 'secondary', secondary.MR);
missLatLon0 = [lat0 - targetLat_truth, lon0 - targetLon_truth];

% ------------------------------------
%%% Setting up for shooting
% ------------------------------------
%%% Don't let future trajectories keep integrating way past the miss point
%%% if they miss Europa
timeLimit_n = time_n(end)*1.05;

%%% Time
t_i = 0; % sec
t_f = timeLimit_n; % Long bc events are watching for impact or escape
n_dt = 10000;
time0_n = linspace(t_i,t_f,n_dt);

%%% Initialize new variables to iterate over
X0_new = X0_X_error;
X_BCR_n_new = X_BCR_n_error;

%%% Preallocate space to track lat/lon impact sites and norm of lat/lon
%%% error
latLons_iter = zeros(n_iterations,2);
missLatLons_iter = zeros(n_iterations,2);

%%% Iterate
for kk = 1:n_iterations
    % ------------------------------------
    %%% Correct state
    % ------------------------------------
    %%% Set error in final state
    if (run_r0error == 1) && (run_v0error == 0) && (run_X0error == 0)
        dXf = X_BCR_n_new(end,1:3)' - rTarget;

    elseif (run_r0error == 0) && (run_v0error == 1) && (run_X0error == 0)
        dXf = X_BCR_n_new(end,4:6)' - vTarget;

    elseif (run_r0error == 0) && (run_v0error == 0) && (run_X0error == 1)
        dXf = X_BCR_n_new(end,1:6)' - stateTarget;
    end

    %%% Grab STM(tf, t0) and invert it to get STM(t0, tf)
    stm_tf_t0 = reshape(X_BCR_n_new(end,7:42),6,6);
    if (run_r0error == 1) && (run_v0error == 0) && (run_X0error == 0)
        stm_tf_t0 = stm_tf_t0(1:3,1:3); % position

    elseif (run_r0error == 0) && (run_v0error == 1) && (run_X0error == 0)
        stm_tf_t0 = stm_tf_t0(4:6,4:6); % velocity
    end
    stm_t0_tf = inv(stm_tf_t0); % state

    %%% Map dXf to t0
    dX0 = stm_t0_tf * dXf;

    %%% Correct X0 with dX0
    if (run_r0error == 1) && (run_v0error == 0) && (run_X0error == 0)
        X0_new = [X0_new(1:6) - [dX0; 0; 0; 0].*dX0_scale; stm0]; % position only
    elseif (run_r0error == 0) && (run_v0error == 1) && (run_X0error == 0)
        X0_new = [X0_new(1:6) - [0; 0; 0; dX0].*dX0_scale; stm0]; % velocity only
    elseif (run_r0error == 0) && (run_v0error == 0) && (run_X0error == 1)
        X0_new = [X0_new(1:6) - dX0.*dX0_scale; stm0]; % full state
    end

    % ------------------------------------
    %%% Integrate new correction
    % ------------------------------------
    %%% Integrate new state
    [time_n, X_BCR_n_new] = ode113(@Int_CR3BnSTM, time0_n, X0_new, options, prms);

    %%% Find new impact conditions
    [lat_new, lon_new] = BCR2latlon(X_BCR_n_new(end,1:3), 'secondary', secondary.MR);
    missLatLons_iter(kk,:) = [lat_new - targetLat_truth, lon_new - targetLon_truth];
    latLons_iter(kk,:) = [lat_new, lon_new];

    %%% Print results from current iteration
    fprintf('-------------------------------------------- %1d\n',kk)
    % fprintf('Latitude Miss: \t%1.4f\n',missLatLonf(1))
    % fprintf('Longitude Miss:\t%1.4f\n',missLatLonf(2))
    fprintf('Norm of Miss:  \t%1.4e\n',norm(missLatLons_iter(kk,:)))
    fprintf('Norm of dXf:  \t%1.4e\n',norm(dXf))

end

% ------------------------------------
%%% Plotting of results
% ------------------------------------
%%% Set color gradient for iterations
[ colorMatrix ] = colorScale([colors.std.mag; colors.std.cyan],n_iterations );

%%% Plot magnitude of lat/lon miss of impact site
figure; hold all
scatter(linspace(1,n_iterations,n_iterations),rowNorm(missLatLons_iter),150,colorMatrix,'filled','marker','o');
p2 = plot([0 n_iterations],[norm(missLatLon0) norm(missLatLon0)],'b','linewidth',2);
legend([p2], 'Norm of initial Lat/Lon miss','location','best')
PlotBoi2('Iteration','Norm of Lat/Lon miss',18,'LaTex')    
setLogPlot('y')
colorbar('Ticks',[0, 1], 'TickLabels',{'First Iter','Last Iter'})
colormap(colorMatrix)


%%% Plot Lat/Lon sites of each iteration
figure; hold all
scatter(latLons_iter(:,2),latLons_iter(:,1),30,colorMatrix,'filled','marker','o')
plot(targetLon_truth, targetLat_truth, 'bo','markersize',15)
PlotBoi2('Longitude, $^\circ$','Latitude, $^\circ$',18,'LaTex')
colorbar('Ticks',[0, 1], 'TickLabels',{'First Iter','Last Iter'})
colormap(colorMatrix)

%%% Plot 3D system with truth and final trajectories
figure; hold all 
p1 = plot3(X_BCR_n_truth(:,1),X_BCR_n_truth(:,2),X_BCR_n_truth(:,3),'b','linewidth',2);
p2 = plot3(X_BCR_n_error(:,1),X_BCR_n_error(:,2),X_BCR_n_error(:,3),'--r','linewidth',2);
p3 = plot3(X_BCR_n_new(:,1),X_BCR_n_new(:,2),X_BCR_n_new(:,3),'color',colorMatrix(end,:),'linewidth',2);
plotBodyTexture3(secondary.R_n, [1-secondary.MR, 0, 0], secondary.img)
axis equal
PlotBoi3('$X_n$','$Y_n$','$Z_n$',18,'LaTex')
legend([p1 p2 p3],'Truth','Perturbed','Final Iteration','location','best')





