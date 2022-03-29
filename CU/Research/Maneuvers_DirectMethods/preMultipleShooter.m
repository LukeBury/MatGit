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
%%% The following are set up to be mutually exclusive (1 = on, 0 = off)
run_r0error = 0; % Perturb & correct only initial position
run_v0error = 0; % Perturb & correct only initial velocity
run_X0error = 1; % Perturb & correct full initial state

% ========================================================================
%%% Setup
% ========================================================================
% ------------------------------------
%%% Preferences
% ------------------------------------
perturbation_X0 = [100; -10; 0; 0.001; -0.01; 0];

n_iterations = 100;

dX0_scale = 1;
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
[time, X_truth] = ode113(@Int_2BI, time0, X0_truth, options, Earth.u);

%%% Setting target
X_target = X_truth(end,1:6)';

% ========================================================================
%%% Perturbing truth and correcting with single shooter
% ========================================================================
% ------------------------------------
%%% Perturbing truth and integrating to get first dXf
% ------------------------------------
%%% Perturbing X0_truth
X0_perturbed = X0_truth + perturbation_X0;

%%% First state transition matrix
stm0 = reshape(eye(6),36,1);
X0_perturbed = [X0_perturbed; stm0];

%%% Propagating truth state to get target state
[time, X_perturbed] = ode113(@Int_2BI_STM, time0, X0_perturbed, options, Earth.u);

% %%% Fist dXf
% dXf = X_perturbed(end,1:6)' - X_target;

%%% Initialize new variables to iterate over
X0_new = X0_perturbed;
X_new = X_perturbed;

%%% Preallocate
dXf_norms = zeros(n_iterations,1);
dX0_norms = zeros(n_iterations,1);

%%% Iterate
for kk = 1:n_iterations
    % ------------------------------------
    %%% Correct state
    % ------------------------------------
    %%% Set error in final state
%     if (run_r0error == 1) && (run_v0error == 0) && (run_X0error == 0)
%         dXf = X_BCR_n_new(end,1:3)' - rTarget;
% 
%     elseif (run_r0error == 0) && (run_v0error == 1) && (run_X0error == 0)
%         dXf = X_BCR_n_new(end,4:6)' - vTarget;
% 
    if (run_r0error == 0) && (run_v0error == 0) && (run_X0error == 1)
        dXf = X_new(end,1:6)' - X_target;
    end
    
    

    %%% Grab STM(tf, t0) and invert it to get STM(t0, tf)
    stm_tf_t0 = reshape(X_new(end,7:42),6,6);
%     if (run_r0error == 1) && (run_v0error == 0) && (run_X0error == 0)
%         stm_tf_t0 = stm_tf_t0(1:3,1:3); % position
% 
%     elseif (run_r0error == 0) && (run_v0error == 1) && (run_X0error == 0)
%         stm_tf_t0 = stm_tf_t0(4:6,4:6); % velocity
%     end
    stm_t0_tf = inv(stm_tf_t0);

    %%% Map dXf to t0
    dX0 = stm_t0_tf * dXf;

    %%% Correct X0 with dX0
%     if (run_r0error == 1) && (run_v0error == 0) && (run_X0error == 0)
%         X0_new = [X0_new(1:6) - [dX0; 0; 0; 0].*dX0_scale; stm0]; % position only
%     elseif (run_r0error == 0) && (run_v0error == 1) && (run_X0error == 0)
%         X0_new = [X0_new(1:6) - [0; 0; 0; dX0].*dX0_scale; stm0]; % velocity only
    if (run_r0error == 0) && (run_v0error == 0) && (run_X0error == 1)
        X0_new = [X0_new(1:6) - dX0.*dX0_scale; stm0]; % full state
    end

    % ------------------------------------
    %%% Integrate new correction
    % ------------------------------------
    %%% Integrate new state
    [time, X_new] = ode113(@Int_2BI_STM, time0, X0_new, options, Earth.u);
    
    % ------------------------------------
    %%% Store results and print
    % ------------------------------------
    dXf_norms(kk,:) = norm(dXf);
    dX0_norms(kk,:) = norm(dX0);
    %%% Print results from current iteration
    fprintf('-------------------------------------------- %1d\n',kk)
    fprintf('Norm of dX0:  \t%1.4e\n',norm(dX0))
    fprintf('Norm of dXf:  \t%1.4e\n',norm(dXf))

end

% ------------------------------------
%%% Plotting of results
% ------------------------------------
%%% Set color gradient for iterations
[ colorMatrix ] = colorScale([colors.std.mag; colors.std.cyan],n_iterations );

%%% Plot system
figure; hold all
p1 = plot3(X_truth(:,1),X_truth(:,2),X_truth(:,3),'b','linewidth',2);
p2 = plot3(X_perturbed(:,1),X_perturbed(:,2),X_perturbed(:,3),'--r','linewidth',2);
p3 = plot3(X_new(:,1),X_new(:,2),X_new(:,3),'color',colorMatrix(end,:),'linewidth',2);
plotBodyTexture3(Earth.R, [0,0,0], Earth.img)
PlotBoi3('$X$','$Y$','$Z$',18,'LaTex')
axis equal
view(0,90)
legend([p1 p2 p3],'Truth','Perturbed','Final Iteration','location','best')

%%% Plot iterations of misses (dX0 and dXf
figure
subplot(1,2,1); hold all
scatter(linspace(1,n_iterations,n_iterations), dX0_norms,100,colorMatrix,'filled','marker','o')
PlotBoi2('Iterations','$|\delta X_0|$',18,'LaTex')
setLogPlot('y')

subplot(1,2,2); hold all
scatter(linspace(1,n_iterations,n_iterations), dXf_norms,100,colorMatrix,'filled','marker','o')
PlotBoi2('Iterations','$|\delta X_f|$',18,'LaTex')
setLogPlot('y')








% % % 
% % % 
% % % %%% Initialize STM and turn into a vector
% % % stm0 = reshape(eye(6),36,1);
% % % 
% % % %%% Add to initial state
% % % X0 = [X0; stm0];
% % % 
% % % %%% Propagate state and STM
% % % [time, X_perturbed] = ode113(@Int_2BI_STM, time0, X0, options, Earth.u);
% % % 
% % % %%% Grab STM(tf, t0) 
% % % stm_tf_t0 = reshape(X_new(end,7:42),6,6);
% % % 
% % % %%% Invert to get STM(t0, tf)
% % % stm_t0_tf = inv(stm_tf_t0);







