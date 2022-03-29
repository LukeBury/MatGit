clear
clc
close all
mbinPath = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/mbin';
addpath(genpath(mbinPath))
tic
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


% ========================================================================
%%% Setup
% ========================================================================
% ------------------------------------
%%% System
% ------------------------------------
%%% Load
primary = bodies.earth;
secondary = bodies.moon;

%%% Replace desired values
secondary.MR = 0.012150585609624;
secondary.a = 384747.962856037;  % km

%%% Normalizing constants
rNorm = secondary.a;         % n <-> km
tNorm = 1/secondary.meanMot; % n <-> sec
vNorm = rNorm / tNorm;       % n <-> km/sec

% t_i = 0; % sec
% t_f = 4*pi; % Long bc events are watching for impact or escape
% n_dt = 10000;
% time0_n = linspace(t_i,t_f,n_dt);
% 

% ------------------------------------
%%% Initial Conditions
% ------------------------------------
X0 = [1.14; 0; -0.16; 0; -0.223; 0];
% X0 = [1.14; 0; -0.163379488757859; 0; -0.223429612636086; 0];

% ------------------------------------
%%% PO Options
% ------------------------------------
%%% Tolerance on PO
tol_PO = 1e-12;

iterLimit = 100;

% ------------------------------------
%%% Integrating initial PO for reference
% ------------------------------------
t_i = 0; % sec
t_f = 2*pi; % Long bc events are watching for impact or escape
n_dt = 100;
time0_n = linspace(t_i,t_f,n_dt);

%%% Choosing ode tolerance
tol = 1e-13;

%%% Initializing stm
stm0 = eye(6);
stm0_vec = reshape(stm0,36,1);
X0 = [X0; stm0_vec];

%%% Setting integrator options
prms.u = secondary.MR;
options = odeset('RelTol',tol,'AbsTol',tol);
options_event = odeset('event',@event_xAxisNegYCrossing,'RelTol',tol,'AbsTol',tol);

%%% Integrating
[time_n, X_BCR_n] = ode113(@Int_CR3BnSTM, time0_n, X0, options_event, prms);

% figure; hold all
% % plot3(X_BCR_n(:,1),X_BCR_n(:,2),X_BCR_n(:,3),'b','linewidth',2)
% plot3(X_BCR_n(:,1),X_BCR_n(:,2),X_BCR_n(:,3),'r','linewidth',2)
% % plotBodyTexture3(primary.R/rNorm,[-secondary.MR, 0, 0],primary.img)
% plotBodyTexture3(secondary.R_n, [1-secondary.MR, 0, 0],secondary.img)
% PlotBoi3('$X_n$','$Y_n$','$Z_n$',18,'LaTex')
% axis equal
% return

% ========================================================================
%%% Iterating to Find PO
% ========================================================================
%%% Setting up variables to be iterated over with their first value
X0_new = X0;
X_new = X_BCR_n;
time0_new = time0_n;
T_new = time_n(end);
error = 100;
iteration = -1;

%%% Preallocate
errors = [];

while (error > tol_PO)
% % %     %%% Set correction at final time (dXf)
% % %     dXf = X_new(end,1:6)' - [X_new(end,1:3)'; 0; X_new(end,5); 0];
% % % 
% % %     %%% Check error and count the iteration
% % %     error = norm([dXf(4), dXf(6)]);
% % %     iteration = iteration + 1;
% % %     
% % %     %%% Grab STM(tf,t0) and invert it to get STM(t0,tf)
% % %     stm_tf_t0 = reshape(X_new(end,7:42),6,6);
% % % 
% % %     stm_t0_tf = inv(stm_tf_t0);
% % %     
% % %     %%% Map dXf to dX0
% % %     dX0 = stm_t0_tf * dXf;
% % %     
% % %     %%% Apply correction to X0
% % %     X0_new = [X0_new(1:6) - [0; 0; dX0(3); 0; dX0(5); 0]; stm0_vec];
    
    
    
% % %     %%% Set correction at final time (dXf)
% % %     dXf = X_new(end,3:6)';
% % %     
% % %     %%% Check error and count the iteration
% % %     error = norm(dXf);
% % %     iteration = iteration + 1;
% % %     
% % %     %%% Grab STM(tf,t0) and invert it to get STM(t0,tf)
% % %     stm_tf_t0 = reshape(X_new(end,7:42),6,6);
% % %     
% % %     stm_tf_t0 = stm_tf_t0(3:6,3:6);
% % %     
% % %     stm_t0_tf = inv(stm_tf_t0);
% % %     
% % %     %%% Map dXf to dX0
% % %     dX0 = stm_t0_tf * dXf;
% % %     
% % %     %%% Apply correction to X0
% % %     X0_new = [X0_new(1:6) - [0; 0; dX0(1); 0; dX0(3); 0]; stm0_vec];
    

    %%% Set correction at final time (dXf)
    dXf = X_new(end,1:6)' - [X_new(end,1:3)'; 0; X_new(end,5); 0];
    
    
    
    
    %%% Check error and count the iteration
    error = norm([dXf(4), dXf(6)]);
    iteration = iteration + 1;
    
    %%% Grab STM(tf,t0) and invert it to get STM(t0,tf)
    stm_tf_t0 = reshape(X_new(end,7:42),6,6);

%     stm_t0_tf = inv(stm_tf_t0);
    
%     [F] = Int_CR3Bn(0,X_new(end,1:6)',prms);

%     mat = [stm_tf_t0(4,3), stm_tf_t0(4,5); stm_tf_t0(6,3), stm_tf_t0(6,5)] - (1/F(5)).*[F(4); F(6)]*[stm_tf_t0(2,3), stm_tf_t0(2,5)];
    mat = [stm_tf_t0(4,3), stm_tf_t0(4,5); stm_tf_t0(6,3), stm_tf_t0(6,5)];

    invMat = inv(mat);
    
    %%% Map dXf to dX0
%     dX0 = stm_t0_tf * dXf;
%     dX0 = invMat * [dXf(4); dXf(6)].*1.9;
    dX0 = mat \ [dXf(4); dXf(6)];
    
    %%% Apply correction to X0
    X0_new = [X0_new(1:6) - [0; 0; dX0(1); 0; dX0(2); 0]; stm0_vec];
    
    
    
    %%% Integate new correction
    [time_new, X_new] = ode113(@Int_CR3BnSTM, time0_n, X0_new, options_event, prms);
    
    
    
    %%% Store results
    errors = [errors; error];
    fprintf('%04d: %1.2e\n',iteration, error)
%     error
%     iteration
%     error
    if iteration > iterLimit
        break
    end
end


% ------------------------------------
%%% Integrate full final PO
% ------------------------------------
timeF_n = linspace(0,time_new(end)*2,n_dt);
[time_PO, X_PO] = ode113(@Int_CR3BnSTM, timeF_n, X0_new, options, prms);

% ========================================================================
%%% Plotting
% ========================================================================
%%% Integrating
[time2_n, X_BCR_n] = ode113(@Int_CR3BnSTM, time_PO, X0, options, prms);

figure; hold all
% plot3(X_BCR_n(:,1),X_BCR_n(:,2),X_BCR_n(:,3),'b','linewidth',2)
plot3(X_BCR_n(:,1),X_BCR_n(:,2),X_BCR_n(:,3),'r','linewidth',2)
% plotBodyTexture3(primary.R/rNorm,[-secondary.MR, 0, 0],primary.img)
plotBodyTexture3(secondary.R_n, [1-secondary.MR, 0, 0],secondary.img)
PlotBoi3('$X_n$','$Y_n$','$Z_n$',18,'LaTex')
axis equal

% figure; hold all
% % plot3(X_BCR_n(:,1),X_BCR_n(:,2),X_BCR_n(:,3),'b','linewidth',2)
% plot3(X_new(:,1),X_new(:,2),X_new(:,3),'r','linewidth',1.5)
% % plotBodyTexture3(primary.R/rNorm,[-secondary.MR, 0, 0],primary.img)
% plotBodyTexture3(secondary.R_n, [1-secondary.MR, 0, 0],secondary.img)
% PlotBoi3('$X_n$','$Y_n$','$Z_n$',18,'LaTex')
% axis equal

figure; hold all
% plot3(X_BCR_n(:,1),X_BCR_n(:,2),X_BCR_n(:,3),'b','linewidth',2)
plot3(X_PO(:,1),X_PO(:,2),X_PO(:,3),'b','linewidth',1.5)
% plotBodyTexture3(primary.R/rNorm,[-secondary.MR, 0, 0],primary.img)
plotBodyTexture3(secondary.R_n, [1-secondary.MR, 0, 0],secondary.img)
PlotBoi3('$X_n$','$Y_n$','$Z_n$',18,'LaTex')
axis equal



toc

function [value, isterminal, direction] = event_xAxisNegYCrossing(t, X, prms)
%%% Description
%       Event function designed to stop trajectory when x-axis is crossed
%       from positive-Y to negative-Y
%       
% --------------------------------------------------------------
%%% Inputs
%       t - [1x1] time
%       X - [nx1] current state
%       prms - [structure] structure of necessary quantities for
%               integration (u, R2_n, L1x, L2x)
% --------------------------------------------------------------
%%% Outputs
%       value - [1xp] relationship to be zero 
%       isterminal - [1xp] should the integration be stopped? (1 or 0)
%       direction - [1xp] sign-direction to watch for on 'value' (1, 0, -1)
% --------------------------------------------------------------
%%% Author
%       Luke Bury, luke.bury@colorado.edu
% ===============================================================
% ------------------------------------
%%% Setting event options
% ------------------------------------
value = X(2);
isterminal = 1;
direction = 0;


end



















