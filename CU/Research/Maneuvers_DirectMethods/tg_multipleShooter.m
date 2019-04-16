clear
clc
close all
mbinPath = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/mbin';
addpath(genpath(mbinPath))

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
%%% Preferences
% ------------------------------------
perturbation_X0 = [100; -10; 0; 0.001; -0.01; 0].*100000;

n_iterations = 15;

% ------------------------------------
%%% System
% ------------------------------------
Earth = bodies.earth;

% ------------------------------------
%%% Truth orbit
% ------------------------------------
%%% Orbital elements
a_truth = Earth.R + 2000;
e_truth = 0.1;
i_truth = 0;
w_truth = 270*pi/180;
raan_truth = 0;
ta_truth = 0;

% ------------------------------------
%%% Integrator options
% ------------------------------------
%%% Time
t_i = 0; % sec
t_f = 0.8*3600; 
n_dt = 1000;
time0 = linspace(t_i,t_f,n_dt);

%%% Choosing ode tolerance
tol = 1e-13;

%%% Setting integrator options
options = odeset('RelTol',tol,'AbsTol',tol);

% ------------------------------------
%%% Integrating truth to obtain a target
% ------------------------------------
%%% Truth state at t0
[r0_truth, v0_truth] = OE2ECI(a_truth, e_truth, i_truth, raan_truth, w_truth, ta_truth, Earth.u);
X0_truth = [r0_truth; v0_truth];

%%% Propagating truth state to get target state
[time_full, X_truth] = ode113(@Int_2BI, time0, X0_truth, options, Earth.u);

%%% Setting target
X_target = X_truth(end,1:6)';

% ========================================================================
%%% Perturbing truth to get a reference, getting nodes, propagating nodes
% ========================================================================
% ------------------------------------
%%% Perturbing truth and integrating to get first dXf
% ------------------------------------
%%% Perturbing X0_truth
X0_perturbed = X0_truth + perturbation_X0;

%%% First state transition matrix
stm0_vec = reshape(eye(6),36,1);
X0_perturbed = [X0_perturbed; stm0_vec];

%%% Propagating truth state to get target state
[time_full, X_perturbed] = ode113(@Int_2BI_STM, time0, X0_perturbed, options, Earth.u);

% ------------------------------------
%%% Getting nodes
% ------------------------------------
%%% Number of nodes
N = 3;

%%% Getting time points of nodes
ms_tf_vec = [0, time_full(end)/N, time_full(end)*2/N, time_full(end)*3/N];

%%% Propagating truth state to get target state
[ms_tf_vec_out, X0s_ms] = ode113(@Int_2BI_STM, ms_tf_vec, X0_perturbed, options, Earth.u);

X0_N0 = [X0s_ms(1,1:6)'; stm0_vec];
X0_N1 = [X0s_ms(2,1:6)'; stm0_vec];
X0_N2 = [X0s_ms(3,1:6)'; stm0_vec];


% figure; hold all
% plot3(X_perturbed(:,1),X_perturbed(:,2),X_perturbed(:,3),'b','linewidth',2)
% plot3(X0_N0(1),X0_N0(2),X0_N0(3),'k.','markersize',15)
% plot3(X1_N0(1),X1_N0(2),X1_N0(3),'k.','markersize',15)
% plot3(X2_N0(1),X2_N0(2),X2_N0(3),'k.','markersize',15)


% ------------------------------------
%%% Propagating nodes
% ------------------------------------
tIn_N0 = linspace(ms_tf_vec(1),ms_tf_vec(2),n_dt);
tIn_N1 = linspace(ms_tf_vec(2),ms_tf_vec(3),n_dt);
tIn_N2 = linspace(ms_tf_vec(3),ms_tf_vec(4),n_dt);

[tOut_N0, X_N0] = ode113(@Int_2BI_STM, tIn_N0, X0_N0, options, Earth.u);

[tOut_N1, X_N1] = ode113(@Int_2BI_STM, tIn_N1, X0_N1, options, Earth.u);

[tOut_N2, X_N2] = ode113(@Int_2BI_STM, tIn_N2, X0_N2, options, Earth.u);

% figure; hold all
% plot3(X_N0(:,1),X_N0(:,2),X_N0(:,3),'b','linewidth',2)
% plot3(X_N1(:,1),X_N1(:,2),X_N1(:,3),'b','linewidth',2)
% plot3(X_N2(:,1),X_N2(:,2),X_N2(:,3),'b','linewidth',2)
% axis equal
% ========================================================================
%%% Finding deltas and begin the multiple shooting
% ========================================================================
X0_N0_new = X0_N0;
X0_N1_new = X0_N1;
X0_N2_new = X0_N2;

X_N0_new = X_N0;
X_N1_new = X_N1;
X_N2_new = X_N2;


figure('position',[-600 384 560 420]); hold all
plot3(X_truth(:,1),X_truth(:,2),X_truth(:,3),'k','linewidth',2)
plot3(X_truth(1,1),X_truth(1,2),X_truth(1,3),'ok','markersize',10)
plot3(X_truth(end,1),X_truth(end,2),X_truth(end,3),'xk','markersize',10)

plot3(X_N0(:,1),X_N0(:,2),X_N0(:,3),'r','linewidth',2)
plot3(X_N0(1,1),X_N0(1,2),X_N0(1,3),'or','markersize',10)
plot3(X_N0(end,1),X_N0(end,2),X_N0(end,3),'xr','markersize',10)

plot3(X_N1(:,1),X_N1(:,2),X_N1(:,3),'g','linewidth',2)
plot3(X_N1(1,1),X_N1(1,2),X_N1(1,3),'og','markersize',10)
plot3(X_N1(end,1),X_N1(end,2),X_N1(end,3),'xg','markersize',10)

plot3(X_N2(:,1),X_N2(:,2),X_N2(:,3),'b','linewidth',2)
plot3(X_N2(1,1),X_N2(1,2),X_N2(1,3),'ob','markersize',10)
plot3(X_N2(end,1),X_N2(end,2),X_N2(end,3),'xb','markersize',10)

plotBodyTexture3(bodies.earth.R,[0,0,0],bodies.earth.img)
axis equal
PlotBoi3('$X$, km','$Y$, km','$Z$, km',18,'LaTex')
view(0,90)

for kk = 1:n_iterations
    % ------------------------------------
    %%% Determine corrections
    % ------------------------------------
    %%% Gather STMs
    stmN0_tf_t0 = reshape(X_N0_new(end,7:42),6,6);
    stmN1_tf_t0 = reshape(X_N1_new(end,7:42),6,6);
    stmN2_tf_t0 = reshape(X_N2_new(end,7:42),6,6);
    
    %%% Invert so they map from tf to t0
    stmN0_t0_tf = inv(stmN0_tf_t0);
    stmN1_t0_tf = inv(stmN1_tf_t0);
    stmN2_t0_tf = inv(stmN2_tf_t0);
    
    %%% Make STM correction matrix
    stm_correctionMatrix = [stmN0_t0_tf, eye(size(stmN0_t0_tf,1)), zeros(size(stmN0_t0_tf,1));...
                            zeros(size(stmN1_t0_tf,1)), stmN1_t0_tf, eye(size(stmN1_t0_tf,1));...
                            zeros(size(stmN1_t0_tf,1)), zeros(size(stmN1_t0_tf,1)), stmN2_t0_tf];
    
    %%% Make dXf correction matrix  
    dXf_N0_new = (X_N0_new(end,1:6) - X_N1_new(1,1:6))';
    dXf_N1_new = (X_N1_new(end,1:6) - X_N2_new(1,1:6))';
    dXf_N2_new = (X_N2_new(end,1:6) - X_target')';
    
    dXf_correctionMatrix = [dXf_N0_new;...
                            dXf_N1_new;...
                            dXf_N2_new];
    
    %%% Multiply to get stack of dX0 corrections
    dX0_correctionMatrix = stm_correctionMatrix * dXf_correctionMatrix; % [n*N x 1]
    
    %%% Grab corrections
    dX0_N0 = dX0_correctionMatrix(1:6);
    dX0_N1 = dX0_correctionMatrix(7:12);
    dX0_N2 = dX0_correctionMatrix(13:18);
    
    %%% Apply corrections
    X0_N0_new = [X0_N0_new(1:6) - dX0_N0; stm0_vec]; % full state
    X0_N1_new = [X0_N1_new(1:6) - dX0_N1; stm0_vec]; % full state
    X0_N2_new = [X0_N2_new(1:6) - dX0_N2; stm0_vec]; % full state
    
    % ------------------------------------
    %%% Integrate new corrections
    % ------------------------------------
    %%% Integrate new states
    [tOut_N0, X_N0_new] = ode113(@Int_2BI_STM, tIn_N0, X0_N0_new, options, Earth.u);
    [tOut_N1, X_N1_new] = ode113(@Int_2BI_STM, tIn_N1, X0_N1_new, options, Earth.u);
    [tOut_N2, X_N2_new] = ode113(@Int_2BI_STM, tIn_N2, X0_N2_new, options, Earth.u);
    
    %%% Print results from current iteration
    fprintf('-------------------------------------------- %1d\n',kk)
    fprintf('Norm of dXf_N0:  \t%1.4e\n',norm(dXf_N0_new))
    fprintf('Norm of dXf_N1:  \t%1.4e\n',norm(dXf_N1_new))
    fprintf('Norm of dXf_N2:  \t%1.4e\n',norm(dXf_N2_new))
    
    
    
    figure(1); hold all
    plot3(X_N0_new(:,1),X_N0_new(:,2),X_N0_new(:,3),'r','linewidth',2)
    plot3(X_N0_new(1,1),X_N0_new(1,2),X_N0_new(1,3),'or','markersize',10)
    plot3(X_N0_new(end,1),X_N0_new(end,2),X_N0_new(end,3),'xr','markersize',10)
    
    plot3(X_N1_new(:,1),X_N1_new(:,2),X_N1_new(:,3),'g','linewidth',2)
    plot3(X_N1_new(1,1),X_N1_new(1,2),X_N1_new(1,3),'og','markersize',10)
    plot3(X_N1_new(end,1),X_N1_new(end,2),X_N1_new(end,3),'xg','markersize',10)
    
    plot3(X_N2_new(:,1),X_N2_new(:,2),X_N2_new(:,3),'b','linewidth',2)
    plot3(X_N2_new(1,1),X_N2_new(1,2),X_N2_new(1,3),'ob','markersize',10)
    plot3(X_N2_new(end,1),X_N2_new(end,2),X_N2_new(end,3),'xb','markersize',10)
    
    
    
    drawnow
    
    
end























