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
print_iterations = 1;

% ========================================================================
%%% Finding family of POs and plotting
% ========================================================================
% ------------------------------------
%%% System
% ------------------------------------
%%% Mass Ratio
prms.u = 0.01;

%%% Time
n_dt = 1000;
t0 = 0;
tf = 6*pi;
time0_n = linspace(t0,tf,n_dt);

%%% Equillibrium points
L123 = EquilibriumPoints(prms.u,1:3); % [3x3] of L1, L2, L3 normalized BCR coordinates

% ------------------------------------
%%% Setting up perturbations along x-axis and looping through
% ------------------------------------
%%% Different initial perturbations to loop through
r0_perts = linspace(0.001, 0.01, 10);

%%% Preallocating space for PO ICs
PO_X0s = zeros(7,length(r0_perts));

for PO_number = 1:length(r0_perts)


%%% IC
X0 = [L123(1,1)-r0_perts(PO_number); 0; 0; 0; 0.1; 0];
stm0 = eye(6);
stm0_vec = reshape(stm0,36,1);
X0 = [X0; stm0_vec];

% ------------------------------------
%%% PO Options
% ------------------------------------
%%% Tolerance on PO
tol_PO = 1e-12;

iterLimit = 100;

% ------------------------------------
%%% Integrator Options
% ------------------------------------
tol = 1e-12;
options = odeset('RelTol',tol,'AbsTol',tol);
options_event1 = odeset('event',@event_xAxisNegYCrossing,'RelTol',tol,'AbsTol',tol);


% ------------------------------------
%%% Integrating initial orbit
% ------------------------------------
%%% Integrating
[time_n, X_BCR_n] = ode113(@Int_CR3BnSTM, time0_n, X0, options_event1, prms);
989
return
if PO_number == 1
figure; hold all
% plot3(X_BCR_n(:,1),X_BCR_n(:,2),X_BCR_n(:,3),'b','linewidth',2)
plot3(X_BCR_n(:,1),X_BCR_n(:,2),X_BCR_n(:,3),'r','linewidth',2)
plot3(L123(1,1),0,0,'^','markeredgecolor',colors.std.black,'markerfacecolor',colors.std.ltgrn)
PlotBoi3('$X_n$','$Y_n$','$Z_n$',18,'LaTex')
axis equal
end


%% =======================================================================
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
    %%% Set correction at final time (dXf)
    dXf = X_new(end,1:6)' - [X_new(end,1:3)'; 0; X_new(end,5); 0];
    
    %%% Check error and count the iteration
    error = norm([dXf(4), dXf(6)]);
    iteration = iteration + 1;
    
    %%% Grab STM(tf,t0) and invert it to get STM(t0,tf)
    stm_tf_t0 = reshape(X_new(end,7:42),6,6);
    mat = [stm_tf_t0(4,3), stm_tf_t0(4,5); stm_tf_t0(6,3), stm_tf_t0(6,5)];

    invMat = inv(mat);
    
    %%% Map dXf to dX0
    dX0 = mat \ [dXf(4); dXf(6)];
    
    %%% Apply correction to X0
    X0_new = [X0_new(1:6) - [0; 0; dX0(1); 0; dX0(2); 0]; stm0_vec];
    
    %%% Integate new correction
    [time_new, X_new] = ode113(@Int_CR3BnSTM, time0_n, X0_new, options_event1, prms);
    
    %%% Store results
    errors = [errors; error];
    
    if print_iterations == 1
        fprintf('%04d: %1.2e\n',iteration, error)
    end

    if iteration > iterLimit
        break
    end
end

PO_X0s(:,PO_number) = [X0_new(1:6); time_new(end)*2];

end % PO_number = 1:5

figure; hold all
[~, X_temp] = ode113(@Int_CR3Bn, linspace(0,PO_X0s(7,1),1000), PO_X0s(1:6,1), options_event1, prms);
plot3(X_temp(:,1),X_temp(:,2),X_temp(:,3),'r','linewidth',2)
plot3(L123(1,1),0,0,'^','markeredgecolor',colors.std.black,'markerfacecolor',colors.std.ltgrn)
PlotBoi3('$X_n$','$Y_n$','$Z_n$',18,'LaTex')
axis equal
    
figure; hold all
for PO_number = 1:length(r0_perts)
    time_n = linspace(0,PO_X0s(7,PO_number),1000);
    [~, X_PO] = ode113(@Int_CR3Bn, time_n, PO_X0s(1:6,PO_number), options, prms);
    p1 = plot3(X_PO(:,1),X_PO(:,2),X_PO(:,3),'b','linewidth',2);
end
p2 = plot3(L123(1,1),0,0,'^','markeredgecolor',colors.std.black,'markerfacecolor',colors.std.ltgrn);
legend([p1, p2], 'Periodic Orbit Family','L_1')
PlotBoi3('$X_n$','$Y_n$','$Z_n$',18,'LaTex')
axis equal


%% =======================================================================
%%% Calculating and plotting stable and unstable manifolds for a PO
% ========================================================================
%%% Picking a PO and chopping it into 10 discrete points
PO_choice = 3;
n_points = 20;
time_PO_choice = linspace(0,PO_X0s(7,PO_choice),1000);
time_points_in = linspace(0,PO_X0s(7,PO_choice),n_points);
[~, X_PO_choice] = ode113(@Int_CR3Bn, time_PO_choice, PO_X0s(1:6,PO_choice), options, prms);
[time_points, X_PO_points] = ode113(@Int_CR3Bn, time_points_in, PO_X0s(1:6,PO_choice), options, prms);

%%% Adding perturbation to each point
PO_pert = 0.001;
X_PO_perturbedPoints_1 = X_PO_points + [PO_pert, PO_pert, 0, 0, 0, 0];
X_PO_perturbedPoints_2 = X_PO_points - [PO_pert, PO_pert, 0, 0, 0, 0];

%%% Propagate manifolds
time_unstable = linspace(0,0.553*pi,1000);
time_stable   = linspace(0.52*pi,0,1000);
figure; hold all
for kk = 1:n_points
    [~, X_PO_unstableManifold_1] = ode113(@Int_CR3Bn, time_unstable, X_PO_perturbedPoints_1(kk,:)', options, prms);
    [~, X_PO_stableManifold_1] = ode113(@Int_CR3Bn, time_stable, X_PO_perturbedPoints_1(kk,:)', options, prms);
    
    [~, X_PO_unstableManifold_2] = ode113(@Int_CR3Bn, time_unstable, X_PO_perturbedPoints_2(kk,:)', options, prms);
    [~, X_PO_stableManifold_2] = ode113(@Int_CR3Bn, time_stable, X_PO_perturbedPoints_2(kk,:)', options, prms);
    
    p2 = plot3(X_PO_unstableManifold_1(:,1),X_PO_unstableManifold_1(:,2),X_PO_unstableManifold_1(:,3),'r','linewidth',2);
    p3 = plot3(X_PO_stableManifold_1(:,1),X_PO_stableManifold_1(:,2),X_PO_stableManifold_1(:,3),'g','linewidth',2);
    
    plot3(X_PO_unstableManifold_2(:,1),X_PO_unstableManifold_2(:,2),X_PO_unstableManifold_2(:,3),'r','linewidth',2);
    plot3(X_PO_stableManifold_2(:,1),X_PO_stableManifold_2(:,2),X_PO_stableManifold_2(:,3),'g','linewidth',2);
end
p1 = plot3(X_PO_choice(:,1),X_PO_choice(:,2),X_PO_choice(:,3),'k','linewidth',2);
p4 = plot3(L123(1,1),0,0,'^','markeredgecolor',colors.std.black,'markerfacecolor',colors.std.ltgrn);
legend([p1, p2, p3, p4], 'Periodic Orbit','Unstable Manifolds','Stable Manifolds','L_1')
PlotBoi3('$X_n$','$Y_n$','$Z_n$',18,'LaTex')
axis equal


%% =======================================================================
%%% Functions
% ========================================================================
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



















