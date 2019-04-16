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
%%% System
% ------------------------------------
primary = bodies.jupiter;
secondary = bodies.europa;

L123 = EquilibriumPoints(secondary.MR,1:3); % [3x3] of L1, L2, L3 normalized BCR coordinates

%%% Normalizing constants
rNorm = secondary.a;         % n <-> km
tNorm = 1/secondary.meanMot; % n <-> sec
vNorm = rNorm / tNorm;       % n <-> km/sec

% ------------------------------------
%%% Preferences
% ------------------------------------
perturbation_r0 = [0; -1; 1]./rNorm;
perturbation_v0 = [-0.001; 0; 0]./vNorm;
perturbation_X0 = [perturbation_r0; perturbation_v0];

n_iterations = 50;

% ------------------------------------
%%% Truth orbit
% ------------------------------------
%%% Set an initial truth state - this initial condition set impacts Europa
%%% at (LatLon [11.713313624820113, -17.884817025059917]).
X0_truth = [1.0204617015266166,-0.0019850725823437,0.0000000000000000,-0.0010976413130697,-0.0030157447223195,0.0098771743854318]';

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
%%% Integrating truth to obtain a target
% ------------------------------------
%%% Integrate truth
[time_truth, X_BCR_truth] = ode113(@Int_CR3Bn, time0_n, X0_truth, options, prms);

%%% Finding truth impact site of trajectory
[targetLat_truth, targetLon_truth] = BCR2latlon(X_BCR_truth(end,1:3), 'secondary', secondary.MR);

%%% Hard-coding target state
X_BCR_Target = X_BCR_truth(end,1:6)';

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
[time_full, X_BCR_perturbed] = ode113(@Int_CR3BnSTM, time0_n, X0_perturbed, options, prms);
% [time_full, X_BCR_perturbed] = ode113(@Int_CR3BnSTM, time_truth, X0_perturbed, options, prms);

% ------------------------------------
%%% Getting nodes
% ------------------------------------
%%% Number of nodes
N = 6;

%%% Getting time points of nodes
% t0_nodes = [0, time_full(end)/N, time_full(end)*2/N, time_full(end)*3/N, time_full(end)*4/N, time_full(end)*5/N, time_full(end)*6/N];
% t0_nodes = [0, time_truth(end)/N, time_truth(end)*2/N, time_truth(end)*3/N, time_truth(end)*4/N, time_truth(end)*5/N, time_truth(end)*6/N];
t0_nodes = [0, time_truth(end)/N, time_truth(end)*2/N, time_truth(end)*3/N, time_truth(end)*4/N, time_truth(end)*(9/10), time_truth(end)*6/N];

%%% Propagating perturbed state to get states at nodes
[t_nodes, X0s_nodes] = ode113(@Int_CR3BnSTM, t0_nodes, X0_perturbed, options, prms);

X0_N1 = [X0s_nodes(1,1:6)'; stm0_vec];
X0_N2 = [X0s_nodes(2,1:6)'; stm0_vec];
X0_N3 = [X0s_nodes(3,1:6)'; stm0_vec];
X0_N4 = [X0s_nodes(4,1:6)'; stm0_vec];
X0_N5 = [X0s_nodes(5,1:6)'; stm0_vec];
X0_N6 = [X0s_nodes(6,1:6)'; stm0_vec];

figure; hold all
plot3(X_BCR_perturbed(:,1),X_BCR_perturbed(:,2),X_BCR_perturbed(:,3),'b','linewidth',2)
plot3(X0_N1(1),X0_N1(2),X0_N1(3),'k.','markersize',15)
plot3(X0_N2(1),X0_N2(2),X0_N2(3),'k.','markersize',15)
plot3(X0_N3(1),X0_N3(2),X0_N3(3),'k.','markersize',15)
plot3(X0_N4(1),X0_N4(2),X0_N4(3),'k.','markersize',15)
plot3(X0_N5(1),X0_N5(2),X0_N5(3),'k.','markersize',15)
plot3(X0_N6(1),X0_N6(2),X0_N6(3),'k.','markersize',15)
plotBodyTexture3(secondary.R_n, [1-secondary.MR, 0, 0],secondary.img)
axis equal
PlotBoi3('$X_n$','$Y_n$','$Z_n$',18,'LaTex')


% ------------------------------------
%%% Propagating nodes
% ------------------------------------
tIn_N1 = linspace(t0_nodes(1),t0_nodes(2),n_dt);
tIn_N2 = linspace(t0_nodes(2),t0_nodes(3),n_dt);
tIn_N3 = linspace(t0_nodes(3),t0_nodes(4),n_dt);
tIn_N4 = linspace(t0_nodes(4),t0_nodes(5),n_dt);
tIn_N5 = linspace(t0_nodes(5),t0_nodes(6),n_dt);
tIn_N6 = linspace(t0_nodes(6),t0_nodes(7),n_dt);

[tOut_N1, X_N1] = ode113(@Int_CR3BnSTM, tIn_N1, X0_N1, options, prms);
[tOut_N2, X_N2] = ode113(@Int_CR3BnSTM, tIn_N2, X0_N2, options, prms);
[tOut_N3, X_N3] = ode113(@Int_CR3BnSTM, tIn_N3, X0_N3, options, prms);
[tOut_N4, X_N4] = ode113(@Int_CR3BnSTM, tIn_N4, X0_N4, options, prms);
[tOut_N5, X_N5] = ode113(@Int_CR3BnSTM, tIn_N5, X0_N5, options, prms);
[tOut_N6, X_N6] = ode113(@Int_CR3BnSTM, tIn_N6, X0_N6, options, prms);

figure; hold all
plot3(X_N1(:,1),X_N1(:,2),X_N1(:,3),'b','linewidth',2)
plot3(X_N2(:,1),X_N2(:,2),X_N2(:,3),'b','linewidth',2)
plot3(X_N3(:,1),X_N3(:,2),X_N3(:,3),'b','linewidth',2)
plot3(X_N4(:,1),X_N4(:,2),X_N4(:,3),'b','linewidth',2)
plot3(X_N5(:,1),X_N5(:,2),X_N5(:,3),'b','linewidth',2)
plot3(X_N6(:,1),X_N6(:,2),X_N6(:,3),'b','linewidth',2)
plotBodyTexture3(secondary.R_n, [1-secondary.MR, 0, 0],secondary.img)
PlotBoi3('$X_n$','$Y_n$','$Z_n$',18,'LaTex')
axis equal
% ========================================================================
%%% Finding deltas and begin the multiple shooting
% ========================================================================
X0_N1_new = X0_N1;
X0_N2_new = X0_N2;
X0_N3_new = X0_N3;
X0_N4_new = X0_N4;
X0_N5_new = X0_N5;
X0_N6_new = X0_N6;

X_N1_new = X_N1;
X_N2_new = X_N2;
X_N3_new = X_N3;
X_N4_new = X_N4;
X_N5_new = X_N5;
X_N6_new = X_N6;


figure('position',[-600 384 560 420]); hold all
plot3(X_BCR_truth(:,1),X_BCR_truth(:,2),X_BCR_truth(:,3),'k','linewidth',2)
plot3(X_BCR_truth(1,1),X_BCR_truth(1,2),X_BCR_truth(1,3),'ok','markersize',10)
plot3(X_BCR_truth(end,1),X_BCR_truth(end,2),X_BCR_truth(end,3),'xk','markersize',10)

plot3(X_N1(:,1),X_N1(:,2),X_N1(:,3),'r','linewidth',2)
plot3(X_N1(1,1),X_N1(1,2),X_N1(1,3),'or','markersize',10)
plot3(X_N1(end,1),X_N1(end,2),X_N1(end,3),'xr','markersize',10)

plot3(X_N2(:,1),X_N2(:,2),X_N2(:,3),'g','linewidth',2)
plot3(X_N2(1,1),X_N2(1,2),X_N2(1,3),'og','markersize',10)
plot3(X_N2(end,1),X_N2(end,2),X_N2(end,3),'xg','markersize',10)

plot3(X_N3(:,1),X_N3(:,2),X_N3(:,3),'b','linewidth',2)
plot3(X_N3(1,1),X_N3(1,2),X_N3(1,3),'ob','markersize',10)
plot3(X_N3(end,1),X_N3(end,2),X_N3(end,3),'xb','markersize',10)

plot3(X_N4(:,1),X_N4(:,2),X_N4(:,3),'m','linewidth',2)
plot3(X_N4(1,1),X_N4(1,2),X_N4(1,3),'ob','markersize',10)
plot3(X_N4(end,1),X_N4(end,2),X_N4(end,3),'xb','markersize',10)

plot3(X_N5(:,1),X_N5(:,2),X_N5(:,3),'c','linewidth',2)
plot3(X_N5(1,1),X_N5(1,2),X_N5(1,3),'ob','markersize',10)
plot3(X_N5(end,1),X_N5(end,2),X_N5(end,3),'xb','markersize',10)

plot3(X_N6(:,1),X_N6(:,2),X_N6(:,3),'p','linewidth',2)
plot3(X_N6(1,1),X_N6(1,2),X_N6(1,3),'ob','markersize',10)
plot3(X_N6(end,1),X_N6(end,2),X_N6(end,3),'xb','markersize',10)

plotBodyTexture3(secondary.R_n, [1-secondary.MR, 0, 0],secondary.img)
axis equal
PlotBoi3('$X_n$','$Y_n$','$Z_n$',18,'LaTex')

for kk = 1:n_iterations
    % ------------------------------------
    %%% Determine corrections
    % ------------------------------------
    %%% Gather STMs
    stmN1_tf_t0 = reshape(X_N1_new(end,7:42),6,6);
    stmN2_tf_t0 = reshape(X_N2_new(end,7:42),6,6);
    stmN3_tf_t0 = reshape(X_N3_new(end,7:42),6,6);
    stmN4_tf_t0 = reshape(X_N4_new(end,7:42),6,6);
    stmN5_tf_t0 = reshape(X_N5_new(end,7:42),6,6);
    stmN6_tf_t0 = reshape(X_N6_new(end,7:42),6,6);
    
    %%% Invert so they map from tf to t0
    stmN1_t0_tf = inv(stmN1_tf_t0);
    stmN2_t0_tf = inv(stmN2_tf_t0);
    stmN3_t0_tf = inv(stmN3_tf_t0);
    stmN4_t0_tf = inv(stmN4_tf_t0);
    stmN5_t0_tf = inv(stmN5_tf_t0);
    stmN6_t0_tf = inv(stmN6_tf_t0);
    
    
    %%% Make STM correction matrix
%     stm_correctionMatrix = [stmN1_t0_tf,                  eye(size(stmN1_t0_tf,1)),     zeros(size(stmN1_t0_tf,1));...
%                             zeros(size(stmN2_t0_tf,1)),   stmN2_t0_tf,                  eye(size(stmN2_t0_tf,1));...
%                             zeros(size(stmN2_t0_tf,1)),   zeros(size(stmN2_t0_tf,1)),   stmN3_t0_tf];
    
    zeros_nxn = zeros(size(stmN2_t0_tf,1));
    eye_nxn = eye(size(stmN1_t0_tf,1));
    
%     stm_correctionMatrix = [stmN1_t0_tf, eye_nxn, zeros_nxn;...
%                             zeros_nxn, stmN2_t0_tf, eye_nxn;...
%                             zeros_nxn, zeros_nxn, stmN3_t0_tf];
                        
    stm_correctionMatrix = [stmN1_t0_tf, eye_nxn,     zeros_nxn,   zeros_nxn,   zeros_nxn,   zeros_nxn;...
                            zeros_nxn,   stmN2_t0_tf, eye_nxn,     zeros_nxn,   zeros_nxn,   zeros_nxn;...
                            zeros_nxn,   zeros_nxn,   stmN3_t0_tf, eye_nxn,     zeros_nxn,   zeros_nxn;...
                            zeros_nxn,   zeros_nxn,   zeros_nxn,   stmN4_t0_tf, eye_nxn,     zeros_nxn;...
                            zeros_nxn,   zeros_nxn,   zeros_nxn,   zeros_nxn,   stmN5_t0_tf, eye_nxn;...
                            zeros_nxn,   zeros_nxn,   zeros_nxn,   zeros_nxn,   zeros_nxn,   stmN6_t0_tf];
          
    %%% Make dXf correction matrix  
    dXf_N1_new = (X_N1_new(end,1:6) - X_N2_new(1,1:6))';
    dXf_N2_new = (X_N2_new(end,1:6) - X_N3_new(1,1:6))';
    dXf_N3_new = (X_N3_new(end,1:6) - X_N4_new(1,1:6))';
    dXf_N4_new = (X_N4_new(end,1:6) - X_N5_new(1,1:6))';
    dXf_N5_new = (X_N5_new(end,1:6) - X_N6_new(1,1:6))';
    dXf_N6_new = (X_N6_new(end,1:6) - X_BCR_Target')';
    
    dXf_correctionMatrix = [dXf_N1_new;...
                            dXf_N2_new;...
                            dXf_N3_new;...
                            dXf_N4_new;...
                            dXf_N5_new;...
                            dXf_N6_new];
    
    %%% Multiply to get stack of dX0 corrections
    dX0_correctionMatrix = stm_correctionMatrix * dXf_correctionMatrix; % [n*N x 1]
    
    %%% Grab corrections
    dX0_N1 = dX0_correctionMatrix(1:6);
    dX0_N2 = dX0_correctionMatrix(7:12);
    dX0_N3 = dX0_correctionMatrix(13:18);
    dX0_N4 = dX0_correctionMatrix(19:24);
    dX0_N5 = dX0_correctionMatrix(25:30);
    dX0_N6 = dX0_correctionMatrix(31:36);
    
    %%% Apply corrections
    X0_N1_new = [X0_N1_new(1:6) - dX0_N1; stm0_vec]; % full state
    X0_N2_new = [X0_N2_new(1:6) - dX0_N2; stm0_vec]; % full state
    X0_N3_new = [X0_N3_new(1:6) - dX0_N3; stm0_vec]; % full state
    X0_N4_new = [X0_N4_new(1:6) - dX0_N4; stm0_vec]; % full state
    X0_N5_new = [X0_N5_new(1:6) - dX0_N5; stm0_vec]; % full state
    X0_N6_new = [X0_N6_new(1:6) - dX0_N6; stm0_vec]; % full state
    
    % ------------------------------------
    %%% Integrate new corrections
    % ------------------------------------
    %%% Integrate new states
    [tOut_N1, X_N1_new] = ode113(@Int_CR3BnSTM, tIn_N1, X0_N1_new, options, prms);
    [tOut_N2, X_N2_new] = ode113(@Int_CR3BnSTM, tIn_N2, X0_N2_new, options, prms);
    [tOut_N3, X_N3_new] = ode113(@Int_CR3BnSTM, tIn_N3, X0_N3_new, options, prms);
    [tOut_N4, X_N4_new] = ode113(@Int_CR3BnSTM, tIn_N4, X0_N4_new, options, prms);
    [tOut_N5, X_N5_new] = ode113(@Int_CR3BnSTM, tIn_N5, X0_N5_new, options, prms);
    [tOut_N6, X_N6_new] = ode113(@Int_CR3BnSTM, tIn_N6, X0_N6_new, options, prms);
    
    %%% Print results from current iteration
    fprintf('-------------------------------------------- %1d\n',kk)
    fprintf('Norm of dXf_N1:  \t%1.4e\n',norm(dXf_N1_new))
    fprintf('Norm of dXf_N2:  \t%1.4e\n',norm(dXf_N2_new))
    fprintf('Norm of dXf_N3:  \t%1.4e\n',norm(dXf_N3_new))
    fprintf('Norm of dXf_N4:  \t%1.4e\n',norm(dXf_N4_new))
    fprintf('Norm of dXf_N5:  \t%1.4e\n',norm(dXf_N5_new))
    fprintf('Norm of dXf_N6:  \t%1.4e\n',norm(dXf_N6_new))
    
    
    
    plot3(X_N1_new(:,1),X_N1_new(:,2),X_N1_new(:,3),'r','linewidth',2)
    plot3(X_N1_new(1,1),X_N1_new(1,2),X_N1_new(1,3),'or','markersize',10)
    plot3(X_N1_new(end,1),X_N1_new(end,2),X_N1_new(end,3),'xr','markersize',10)
    
    plot3(X_N2_new(:,1),X_N2_new(:,2),X_N2_new(:,3),'g','linewidth',2)
    plot3(X_N2_new(1,1),X_N2_new(1,2),X_N2_new(1,3),'og','markersize',10)
    plot3(X_N2_new(end,1),X_N2_new(end,2),X_N2_new(end,3),'xg','markersize',10)
    
    plot3(X_N3_new(:,1),X_N3_new(:,2),X_N3_new(:,3),'b','linewidth',2)
    plot3(X_N3_new(1,1),X_N3_new(1,2),X_N3_new(1,3),'ob','markersize',10)
    plot3(X_N3_new(end,1),X_N3_new(end,2),X_N3_new(end,3),'xb','markersize',10)
    
    plot3(X_N4_new(:,1),X_N4_new(:,2),X_N4_new(:,3),'m','linewidth',2)
    plot3(X_N4_new(1,1),X_N4_new(1,2),X_N4_new(1,3),'om','markersize',10)
    plot3(X_N4_new(end,1),X_N4_new(end,2),X_N4_new(end,3),'xm','markersize',10)
    
    plot3(X_N5_new(:,1),X_N5_new(:,2),X_N5_new(:,3),'c','linewidth',2)
    plot3(X_N5_new(1,1),X_N5_new(1,2),X_N5_new(1,3),'oc','markersize',10)
    plot3(X_N5_new(end,1),X_N5_new(end,2),X_N5_new(end,3),'xc','markersize',10)
    
    plot3(X_N6_new(:,1),X_N6_new(:,2),X_N6_new(:,3),'p','linewidth',2)
    plot3(X_N6_new(1,1),X_N6_new(1,2),X_N6_new(1,3),'op','markersize',10)
    plot3(X_N6_new(end,1),X_N6_new(end,2),X_N6_new(end,3),'xp','markersize',10)
    
    
    drawnow
    
    
end























