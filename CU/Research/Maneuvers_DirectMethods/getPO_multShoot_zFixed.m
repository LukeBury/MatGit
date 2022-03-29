% ========================================================================
%%% Description
% ========================================================================
% Finding POs by also correcting time vector

% Created: 
% Author : Luke Bury, luke.bury@colorado.edu
% ========================================================================
%%% Initialization
% ========================================================================
clear
clc
close all
mbinPath = '~/CU_Google_Drive/Documents/MatGit/mbin';
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
run_CR3BP              = 1;
run_CR3BP_J2pJ4pJ6pJ2s = 0;


if (run_CR3BP == 1) && (run_CR3BP_J2pJ4pJ6pJ2s == 1)
    warning('These are mutually exclusive options')
    return
end

% ========================================================================
%%% Setup
% ========================================================================
% ------------------------------------------------- 
%%% Initial conditions
% -------------------------------------------------
% % --------------------------
% % X0 and Tp corresponding to L2 Lyapunov biffurcation at Europa
% % --------------------------
% primary = bodies.jupiter;   secondary = bodies.europa;
% X0_n = [1.017026611070739; 0.000000014441540; 0;...
%         0.000000002451836; 0.020022207280852; 0];
% T0_PO_n = 3.124263922352081;



% % --------------------------
% % X0 and Tp corresponding to L1 Lyapunov biffurcation at Moon
% % --------------------------
% primary = bodies.earth;     secondary = bodies.moon;
% X0_n = [0.853072381084468;
%  -0.000000000000000;
%  0.000000000000000;
%  -0.000000000000001;
%  -0.121880628358625;
%  0.000000000000000];
% T0_PO_n =  2.734302452595165;

% % --------------------------
% % X0 and Tp corresponding to L2 Lyapunov biffurcation at Moon
% % --------------------------
% primary = bodies.earth;     secondary = bodies.moon;
% X0_n = [1.120360542366735; 0.000036610925873; 0;...
%         0.000002187136222; 0.176125048984917; 0];
% T0_PO_n = 3.415558008018916;


% % --------------------------
% % X0 and Tp corresponding to L2 Vertical biffurcation at Europa
% % --------------------------
% primary = bodies.jupiter;   secondary = bodies.europa;
% X0_n = [1.015816379745207; 0.000000000411681; 0.000000001297984;...
%         0.000000002061456; -0.017939765578591; -0.056562126699769];
% T0_PO_n = 4.266711475227815;
% 


% % --------------------------
% % X0 and Tp corresponding to L2 Vertical biffurcation at ganymede
% % --------------------------
% primary = bodies.jupiter;   secondary = bodies.ganymede;
% X0_n = [1.033657240017886;
%  -0.000000000000000;
%  0.000000000000000;
%  0.000000000000012;
%  -0.027989575962354;
%  0.000000000000000];
% T0_PO_n = 3.148710407660965;

% % --------------------------
% % X0 and Tp corresponding to L2 Lyapunov biffurcation at Enceladus
% % --------------------------
% primary = bodies.saturn;    secondary = bodies.enceladus;
% X0_n = [1.003390444557048; -0.000805481191125; 0;...
%         -0.000224307203511; 0.003235968095636; 0];
% T0_PO_n = 3.086292283602011;
% 
% % --------------------------
% % X0 and Tp corresponding to L2 Vertical biffurcation at Enceladus
% % --------------------------
% primary = bodies.saturn;    secondary = bodies.enceladus;
% X0_n = [1.003104911775108; 0; -0.000000000003900;...
%         0.000000000011816; -0.003397299105638; 0.011089587166677];
% T0_PO_n = 4.248617774497027;




% % --------------------------
% % X0 and Tp corresponding to L4-untitled1 biffurcation at Enceladus
% % --------------------------
% primary = bodies.saturn;    secondary = bodies.enceladus;
% X0_n = [0.564947073975989; 0.302902430320629; -0.877900340143571;...
%         -0.454995222410374; -0.108438133722420; -0.238901417152178];
% T0_PO_n = 6.283185489835216;
% %^^^ may want a perturbation on the order of 1e-2









% % --------------------------
% % X0 and Tp corresponding to unlabeled1TrailingS biffurcation 1 at Enceladus
% % --------------------------
% primary = bodies.saturn;    secondary = bodies.enceladus;
% X0_n = [1.003867374501306; 0; -0.001615356503498;...
%         0.007420489615297; -0.010140084240107; 0.004290008196305];
% T0_PO_n = 4.449570709517759;

% % --------------------------
% % X0 and Tp corresponding to unlabeled1TrailingS biffurcation 2 at Enceladus
% % --------------------------
% primary = bodies.saturn;    secondary = bodies.enceladus;
% X0_n = [1.003328461055074; 0; -0.001179557781950;...
%         0.007734754444788; -0.009573295236007; 0.004683361178424];
% T0_PO_n = 4.055805113650098;

% % --------------------------
% % X0 and Tp corresponding to unlabeled1LeadingS biffurcation 1 at Enceladus
% % --------------------------
% primary = bodies.saturn;    secondary = bodies.enceladus;
% X0_n = [; 0; ;...
%         ; ; ];
% T0_PO_n = ;

% % --------------------------
% % X0 and Tp corresponding to unlabeled1LeadingS biffurcation 2 at Enceladus
% % --------------------------
% primary = bodies.saturn;    secondary = bodies.enceladus;
% X0_n = [; 0; ;...
%         ; ; ];
% T0_PO_n = ;


% % --------------------------
% % X0 and Tp corresponding to L2 Lyapunov biffurcation at Triton
% % --------------------------
primary = bodies.neptune;    secondary = bodies.triton;
X0_n = [1.046598413405509;
 -0.000000000000000;
 0.000000000000000;
 0.000000000001725;
 -0.036447831477736;
 0.000000000000000];
T0_PO_n = 3.164931597055242;








% % --------------------------
% % X0 and Tp corresponding to L2 Vertical biffurcation at Titan
% % --------------------------
% primary = bodies.saturn;    secondary = bodies.titan;
% X0_n = [1.033297105540779; 0; 0.000000000000013;...
%         -0.000000000000012; -0.037907926186969; 0.117177360844365];
% T0_PO_n = 4.226170603679089;



% % --------------------------
% % X0 and Tp corresponding to L2 Lyap at Cordelia
% % --------------------------
% primary = bodies.uranus;    secondary = bodies.cordelia;
% X0_n = [1.000538491321221; 0; 0;...
%         0; 0.000093855706042; 0];
% T0_PO_n = 3.034205183528419;

% % --------------------------
% % X0 and Tp corresponding to L1 Lyap at Ophelia
% % --------------------------
% primary = bodies.uranus;    secondary = bodies.ophelia;
% X0_n = [0.999411780886685; 0; 0;...
%         0; 0.000000860848623; 0];
% T0_PO_n = 3.031757208321406;

% % --------------------------
% % X0 and Tp corresponding to early L2 Halo at Europa w/ J2p,J4p,J6p,J2s
% % --------------------------
% primary = bodies.jupiter;    secondary = bodies.europa;
% X0_n = [1.017571399455083;
%          -0.000318686122917;
%          0.000000000000000;
%          -0.000080660391112;
%          0.017152892682995;
%          0.000000000000000];
% T0_PO_n = 3.115795317196816;

% --------------------------
% X0 and Tp corresponding to early L2 Halo at Enceladus w/ J2p,J4p,J6p,J2s
% --------------------------
% primary = bodies.saturn;    secondary = bodies.enceladus;
% % X0_n = [1.004561503146581;
% %          0.000807190907664;
% %          0.000000000000000;
% %          0.000699312550325;
% %          -0.002924761146237;
% %          0.000000000000000];
% % T0_PO_n =  3.223988625149796;



%%% Make Guess
X0_guess_n = X0_n - [0; 0; 1; 0; 0; 0].*1e-5;
% X0_guess_n = X0_n;

Tp_guess_n = T0_PO_n;
% -------------------------------------------------
%%% System
% -------------------------------------------------
%%% Normalizing constants
rNorm = secondary.a;         % n <-> km
tNorm = 1/secondary.meanMot; % n <-> sec
vNorm = rNorm / tNorm;       % n <-> km/sec

%%% Setting parameters
prms.u  = secondary.MR;
prms.R2 = secondary.R_n;
prms.R1 = primary.R / rNorm;

%%% Shortcut variables
mu   = secondary.MR;
R2_n = secondary.R_n;

if run_CR3BP_J2pJ4pJ6pJ2s == 1
    prms.J2p = primary.J2; prms.J4p = primary.J4; prms.J6p = primary.J6;
    prms.J2s = secondary.J2;
end

%%% Equillibrium Points
if run_CR3BP == 1
    rLPs_n = EquilibriumPoints(mu);
elseif run_CR3BP_J2pJ4pJ6pJ2s == 1
    rLPs_n = collinearEquilibriumPoints_ZH(prms);
end

% ------------------------------------------------- 
%%% Integration Options
% -------------------------------------------------
tol = 1e-12;
options = odeset('RelTol',tol,'AbsTol',tol);
t0 = 0;


% ------------------------------------------------- 
%%% Periodic Orbit options
% -------------------------------------------------
%%% Error tolerance for constraint vector
error_tol = 1e-12;

%%% Number of nodes for multiple shooter
n_Nodes = 5; 

%%% Maximum number of iterations for multiple shooter
iterMax = 500;

%%% Step size for multiple shooter (usually 1, sometimes 2 or maybe even 3)
% stepSize = 0.01;
stepSize = 1;

% ------------------------------------------------- 
%%% Integrate and plot guess
% -------------------------------------------------
%%% Setup
stm0 = eye(6);
stm0_vec = reshape(stm0,36,1);

%%% Integrate guess
if run_CR3BP == 1
    [~, X_initial] = ode113(@Int_CR3Bn, linspace(0, Tp_guess_n,2), X0_guess_n, options, prms);
elseif run_CR3BP_J2pJ4pJ6pJ2s == 1
    [~, X_initial] = ode113(@Int_CR3Bn_J2pJ4pJ6pJ2s, linspace(0, Tp_guess_n,2), X0_guess_n, options, prms);
end

%%% Plot Guess
figure; hold all
plot3(X_initial(:,1),X_initial(:,2),X_initial(:,3),'r','linewidth',1.5)
PlotBoi3('$X_n$','$Y_n$','$Z_n$',20,'LaTex')
view(0,90)

% ========================================================================
%%% Create nodes and send the guess to a multiple shooter
% ========================================================================
% ------------------------------------------------- 
%%% Integrate guess and divide into nodes for multiple shooter
% -------------------------------------------------
%%% Integrate again to discretize into nodes evenly spaced in time
if run_CR3BP == 1
    [T_nodes, X_nodes] = ode113(@Int_CR3Bn, linspace(0, Tp_guess_n,n_Nodes+1), X0_guess_n, options, prms);
elseif run_CR3BP_J2pJ4pJ6pJ2s == 1
    [T_nodes, X_nodes] = ode113(@Int_CR3Bn_J2pJ4pJ6pJ2s, linspace(0, Tp_guess_n,n_Nodes+1), X0_guess_n, options, prms);
end

%%% Setting initial nodes (getting rid of repeat end state)
X_nodes = X_nodes(1:n_Nodes,1:6);

%%% Create first free-variable column vector
F_new   = X_nodes';
F_new   = F_new(:);
F_new   = [F_new; Tp_guess_n];

F_new = [F_new(1:2); F_new(4:end)];

% -------------------------------------------------
%%% Pre-While-Loop work
% -------------------------------------------------
%%% Set while-loop/convergence variables
iteration_counter = 0;
constraint_error  = 100;


% -------------------------------------------------
% Enter while loop - attempt to converge on next PO in family
% -------------------------------------------------
while (constraint_error > error_tol) && (iteration_counter < iterMax)

    %%% Count iteration
    iteration_counter = iteration_counter + 1;

    % -------------------------------------------------
    %%% Using multiple shooting - loop through nodes, integrate, and 
    %%% populate constraint vector and create DF matrix
    % -------------------------------------------------  
    %%% Integrate again to discretize into nodes evenly spaced in time
    if run_CR3BP == 1
        [DF_mat, constraints] = multShooter_stateContinuity_zFixed(n_Nodes, X0_guess_n(3), F_new, @Int_CR3BnSTM, options, prms);
    elseif run_CR3BP_J2pJ4pJ6pJ2s == 1
%         [DF_mat, constraints] = multShooter_stateContinuity(n_Nodes, F_new, @Int_CR3BnSTM_J2pJ4pJ6pJ2s, options, prms);
        [DF_mat, constraints] = multShooter_stateContinuity_zFixed(n_Nodes, X0_guess_n(3), F_new, @Int_CR3BnSTM_J2pJ4pJ6pJ2s, options, prms);
        
    end
    
    % -------------------------------------------------
    %%% Determine current error and create next F if necessary
    % -------------------------------------------------
    %%% Compute error
    constraint_error = norm(constraints);

    %%% Compute new free-variable vector if error is not converged
    if (constraint_error > error_tol)
        warning('off','MATLAB:nearlySingularMatrix')
        F_new = F_new - stepSize * DF_mat'*((DF_mat*(DF_mat'))\constraints);
        warning('on','MATLAB:nearlySingularMatrix')
    end
    
end % while constraint_error > error_tol

if constraint_error > error_tol
    warning('Solution not converged on')
end



% ========================================================================
%%% Integrate solution, plot and print result
% ========================================================================
%%% Integrate Solution
if run_CR3BP == 1
    [T_sol_n, X_sol_n] = ode113(@Int_CR3Bn, [0, F_new(end)], [F_new(1:2); X0_guess_n(3); F_new(3:5)], options, prms);
elseif run_CR3BP_J2pJ4pJ6pJ2s == 1
    [T_sol_n, X_sol_n] = ode113(@Int_CR3Bn_J2pJ4pJ6pJ2s, [0, F_new(end)], [F_new(1:2); X0_guess_n(3); F_new(3:5)], options, prms);
end

%%% Plot Solution
figure; hold all
plot3(X_sol_n(:,1),X_sol_n(:,2),X_sol_n(:,3),'r','linewidth',1.5)
PlotBoi3('$X_n$','$Y_n$','$Z_n$',20,'LaTex')
% axis equal

plot3(rLPs_n(2,1), rLPs_n(2,2), rLPs_n(2,3), '^', 'markersize', 8, 'markerfacecolor', colors.std.ltgrn, 'markeredgecolor', colors.std.blue)
% plot3(rLPs_n(1,1), rLPs_n(1,2), rLPs_n(1,3), '^', 'markersize', 8, 'markerfacecolor', colors.std.ltgrn, 'markeredgecolor', colors.std.blue)

% plotSecondary(secondary)

view(0,0)

%%% Print Solution
fprintf('Iterations: %1.0d\n',iteration_counter)
fprintf('Error:      %1.1e\n\n',constraint_error)
prettyColVec([F_new(1:2); X0_guess_n(3); F_new(3:5); F_new(end)])





% ========================================================================
%%% Formatting Structures
% ========================================================================
% -------------------------------------------------
%%% 
% -------------------------------------------------

% --------------------------
% 
% --------------------------

% ========================================================================
%%% Closeout
% ========================================================================
toc(ticWhole)

















