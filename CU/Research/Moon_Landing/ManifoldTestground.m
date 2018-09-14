clear
clc
% close all
addpath(genpath('/Users/lukebury/Documents/MATLAB/mbin'))

% ========================================================================
%%% Run switches
% ========================================================================
run_symbollicWork = 0;
% ========================================================================
%%% Importing Data
% ========================================================================
%%% General data on solar system bodies
bodies = getBodyData();

%%% Color options/schemes
colors = get_colors();

% ========================================================================
%%% Setup
% ========================================================================
% -----------------------------------------------------------
%%% Initial Condition Sets
% -----------------------------------------------------------
% ----------------------
%%% Earth-Moon
% ----------------------
% primary = bodies.earth;
% secondary = bodies.moon;

% %%% L1 - Planar (earth-moon)
% T0 = 2.691624497299714;
% X0 = [0.836571900577714;...
%       0.000000055924169;...
%                       0;...
%       0.000000036748495;...
%       0.002968035126446;...
%                       0];
% T_man = 2*T0; % Manifold propagation time

% %%% L2 - Clamshell (earth-moon)
% T0 = 3.524355367498158;
% X0 = [1.154463466857942;...
%       0.000000000006268;...
%       0.000000000055305;...
%       0.000000000012753;...
%      -0.002743204741188;...
%      -0.068069904553603];
% T_man = 2*T0; % Manifold propagation time

% %%% L3 - Planar (earth-moon)
% T0 = 6.218401292006472;
% X0 = [-1.006379590524185;...
%       -0.000000000012063;...
%                        0;...
%       -0.000000000009193;...
%        0.002662662555309;...
%                        0];
% T_man = 2*T0; % Manifold propagation time

% ----------------------
%%% Jupiter-Europa
% ----------------------
primary = bodies.jupiter;
secondary = bodies.europa;

% %%% L1 - Planar (jupiter-europa) (POs-30, from n = 100)
% T0 = 2.990191959077223;
% X0 = [0.979372439622443;...
%       0.000000025997677;...
%                       0;...
%       0.000000018078227;...
%       0.002723419684919;...
%                       0];
% T_man = 2*T0; % Manifold propagation time

% %%% L1 - Clamshell-500 (jupiter-europa) (from n = 1000)
% T0 = 3.138956088377439;
% X0 = [0.980217769439610;...
%      -0.000000000037079;...
%      -0.000000000645616;...
%      -0.000000000238659;...
%       0.000932691862210;...
%       0.016663538339507];
% T_man = 2*T0; % Manifold propagation time

% %%% L1 - Clamshell-10000 (jupiter-europa)
% T0 = 3.188476570848546;
% X0 = [0.980656399071857;...
%      -0.000000000181723;...
%      -0.000000002212019;...
%      -0.000000001176243;...
%       0.001912360187174;...
%       0.023431862692639];
% T_man = 2*T0; % Manifold propagation time

% %%% L2 - Planar (jupiter-europa)
% T0 = 3.077407327172548;
% X0 = [1.020027801952722;...
%      -0.000000000182340;...
%                       0;...
%      -0.000000000108287;...
%       0.002758777827967;...
%                       0];
% T_man = 2*T0; % Manifold propagation time

% %%% L2 - Clamshell-500 (jupiter-europa) (from n = 1000)
% T0 = 3.234441239111177;
% X0 = [1.019933901325254;...
%       0.000000000102434;...
%       0.000000001583409;...
%       0.000000000634031;...
%      -0.001116489630818;...
%      -0.017437037466781];
% T_man = 2*T0*1.5; % Manifold propagation time

%%% L2 - Clamshell-10000 (jupiter-europa)
T0 = 3.283922417064763;
X0 = [1.019430611934981;...
   0.000000000457065;...
   0.000000004888444;...
   0.000000002816334;...
  -0.002279148538495;...
  -0.024446534776545];
T_man = 2*T0*1.4; % Manifold propagation time

% %%% L3 - Planar (jupiter-europa)
% T0 = 6.283046325015956;
% X0 = [-0.998682368274998;...
%        0.000000004535936;...
%                        0;...
%        0.000000003458249;...
%       -0.002657273182115;...
%                       0];
% T_man = 15*T0; % Manifold propagation time

% ------------------------------------------------- 
%%% Manifold options
% -------------------------------------------------
%%% Manifold options
nMan = 100; % Number of manifolds to compute
% eVec_n = 1;
manPert = 5e-8;

stable_yin    = 0;
unstable_yin  = 1;
stable_yang   = 0;
unstable_yang = 0;

plot_PO = 0;

% ------------------------------------------------- 
%%% Normalizing constants and u (MR)
% -------------------------------------------------
%%% Normalizing constants
rNorm = secondary.a;         % n <-> km
tNorm = 1/secondary.meanMot; % n <-> sec
vNorm = rNorm / tNorm;       % n <-> km/sec

%%% Mass ratio of system
u = secondary.MR;

% ------------------------------------------------- 
%%% Integration Options
% -------------------------------------------------
tol = 1e-9;
options = odeset('RelTol',tol,'AbsTol',tol);

% ------------------------------------------------- 
%%% Finding Equilibrium Points
% -------------------------------------------------
Ls_n = EquilibriumPoints(secondary.MR);

% ========================================================================
%%% Periodic Orbit Integration
% ========================================================================
% ------------------------------------------------- 
%%% Setting up matrices for integration
% -------------------------------------------------
%%% Initial STM
stm0 = eye(6);

%%% Forming full initial conditions
X0 = [X0; reshape(stm0,36,1)];

%%% Setting time vector
time = [0, T0];

% ------------------------------------------------- 
%%% Integrating initial periodic orbit
% -------------------------------------------------
%%% Find Periodic Orbit with ode45
[T_PO, X_PO] = ode45(@int_CR3BnSTM, time, X0, options, u);

%%% Seperating monodromy matrix and getting eigenvalues and vectors
monodromy = reshape(X_PO(end,7:end),6,6);
[eVecs, eVals] = eig(monodromy);

%%% Defining eigenvalues and vectors associated with unstable and stable perturbations
eVal_S = min(real(diag(eVals)));
eVal_U = max(real(diag(eVals)));
eVec_S = eVecs(:,find(diag(eVals)==eVal_S));
eVec_U = eVecs(:,find(diag(eVals)==eVal_U));

% ========================================================================
%%% Finding manifolds
% ========================================================================
%%% Looping through manifolds
for kk = 1:nMan
    %%% Integrating PO for some fraction of its total period
    time = [T0*kk/nMan 0];
    [T_start, X_start] = ode45(@int_CR3BnSTM, time, X0, options, u);
    stm = reshape(X_start(end,7:end),6,6);
    
    % ----------------------------------------- 
    %%% Stable manifold to the "left" of PO
    % -----------------------------------------
    if stable_yin == 1
        %%% Perturb state in stable direction
        X_pert = X_start(end,1:6) + manPert*((stm*eVec_S)./sqrt((stm*eVec_S)'*(stm*eVec_S)))';
        
        %%% Store new ICs
        IC = [X_pert, reshape(stm0,1,36)];
        time = [T_man, 0];
        
        %%% Integrate stable manifold
        [T_int, X] = ode45(@int_CR3BnSTM, time, IC, options, u);
        
        %%% Plot stable manifold
        figure(1); hold all
        plot3(X(:,1),X(:,2),X(:,3),'color',colors.std.black);
        plot3(X(1,1),X(1,2),X(1,3),'o','markersize',5,'color',colors.std.black)
    end
    
    % ----------------------------------------- 
    %%% Unstable manifold to the "left" of PO
    % -----------------------------------------
    if unstable_yin == 1
        %%% Perturb state in unstable direction
        X_pert = X_start(end,1:6) + manPert*(stm*eVec_U/norm(stm*eVec_U))';
        
        %%% Store new ICs
        IC = [X_pert, reshape(stm0,1,36)];
        time = [0, T_man];

        %%% Integrate unstable manifold
        [T_int, X] = ode45(@int_CR3BnSTM, time, IC, options, u);
        
        %%% Plot unstable manifold
        figure(1); hold all
        plot3(X(:,1),X(:,2),X(:,3),'color',colors.std.black);
        plot3(X(1,1),X(1,2),X(1,3),'o','markersize',5,'color',colors.std.black)
    end
    
    % ----------------------------------------- 
    %%% Stable manifold to the "right" of PO
    % -----------------------------------------
    if stable_yang == 1
        %%% Perturb state in stable direction
        X_pert = X_start(end,1:6) - manPert*((stm*eVec_S)./sqrt((stm*eVec_S)'*(stm*eVec_S)))';
        
        %%% Store new ICs
        IC = [X_pert, reshape(stm0,1,36)];
        time = [T_man, 0];
        
        %%% Integrate stable manifold
        [T_int, X] = ode45(@int_CR3BnSTM, time, IC, options, u);
        
        %%% Plot stable manifold
        figure(1); hold all
        plot3(X(:,1),X(:,2),X(:,3),'color',colors.std.black);
        plot3(X(1,1),X(1,2),X(1,3),'o','markersize',5,'color',colors.std.black)
    end
    
    % ----------------------------------------- 
    %%% Unstable manifold to the "left" of PO
    % -----------------------------------------
    if unstable_yang == 1
        %%% Perturb state in unstable direction
        X_pert = X_start(end,1:6) - manPert*(stm*eVec_U/norm(stm*eVec_U))';
        
        %%% Store new ICs
        IC = [X_pert, reshape(stm0,1,36)];
        time = [0, T_man];

        %%% Integrate unstable manifold
        [T_int, X] = ode45(@int_CR3BnSTM, time, IC, options, u);
        
        %%% Plot unstable manifold
        figure(1); hold all
        plot3(X(:,1),X(:,2),X(:,3),'color',colors.std.black);
        plot3(X(1,1),X(1,2),X(1,3),'o','markersize',5,'color',colors.std.black)
    end
    
    
end

%%% Plotting periodic orbit
figure(1)
plot3(X_PO(:,1),X_PO(:,2),X_PO(:,3),'k','linewidth',2)

%%% Plotting bodies
plotBody2(secondary.R_n,[1-secondary.MR,0,0],secondary.color,colors.std.black,0.5,.5)
% plotBody2(primary.R/rNorm,[-secondary.MR,0,0],primary.color,colors.std.black,0.5)
% % plotBodyTexture3(secondary.R_n,[1-secondary.MR,0,0],secondary.img)
% % plotBodyTexture3(primary.R/rNorm,[-secondary.MR,0,0],primary.img)

%%% Plotting lagrange points
% plot3(Ls_n(1,1),Ls_n(1,2),Ls_n(1,3),'^','markerfacecolor',colors.std.red,'markeredgecolor',colors.std.black,'markersize',8)
% plot3(Ls_n(2,1),Ls_n(2,2),Ls_n(2,3),'^','markerfacecolor',colors.std.red,'markeredgecolor',colors.std.black,'markersize',8)

%%% Plot Options
PlotBoi3('$x_n$','$y_n$','$z_n$',16,'LaTex')

if plot_PO == 1
    figure(2); hold all
    plot3(X_PO(:,1),X_PO(:,2),X_PO(:,3),'k','linewidth',2)
    PlotBoi3('X','Y','Z',16)
end

[JC_L1]   = JacobiConstantCalculator(secondary.MR,Ls_n(1,:),[0,0,0])
[JC_L2]   = JacobiConstantCalculator(secondary.MR,Ls_n(2,:),[0,0,0])
[JC_Surf] = JacobiConstantCalculator(secondary.MR,[1-secondary.MR-secondary.R_n,0,0],[0,0,0])
[JC_PO]   = JacobiConstantCalculator(secondary.MR,X_PO(1,1:3),X_PO(1,4:6))







































