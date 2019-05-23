clear
clc
close all
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
% T0 = 2.691623089633314;
% X0 = [0.836571976076091;...
%       0.000000001773828;...
%                       0;...
%       0.000000001165175;...
%       0.002968035186365;...
%                       0];
% T_man = 2*T0; % Manifold propagation time

% %%% L2 - Planar (earth-moon)
% T0 = 3.373250668778623;
% X0 = [1.155132035976244;...
%       0.000000000021998;...
%                       0;...
%       0.000000000013362;...
%       0.002935065880136;...
%                       0];
% T_man = 2*T0; % Manifold propagation time

% %%% L3 - Planar (earth-moon)
% T0 = 6.218403969815443;
% X0 = [-1.006379740628474;...
%       -0.000000000011746;...
%                        0;...
%       -0.000000000009066;...
%        0.002662662551122;...
%                        0];
% T_man = 2*T0; % Manifold propagation time

% ----------------------
%%% Jupiter-Europa
% ----------------------
primary = bodies.jupiter;
secondary = bodies.europa;

% %%% L1 - Planar (jupiter-europa) ... Great demonstration of landing (yin)
% T0 = 2.985780105787740;
% X0 = [0.979401436352991;...
%       0.000000068848856;...
%                       0;...
%       0.000000047890083;...
%       0.002725303957529;...
%                       0];
% T_man = 2*T0; % Manifold propagation time

% %%% L1 - Clamshell-500 (jupiter-europa)
% T0 = 3.133881318523165;
% X0 = [0.980248383388469;...
%      -0.000000000068220;...
%      -0.000000001200126;...
%      -0.000000000445721;...
%       0.000936572582599;...
%       0.016707520453371];
% T_man = 2*T0; % Manifold propagation time

% %%% L1 - Clamshell-10000 (jupiter-europa)
% T0 = 3.183399243940943;
% X0 = [0.980689448569968;...
%      -0.000000000338156;...
%      -0.000000004121808;...
%      -0.000000002202635;...
%       0.001920589557276;...
%       0.023493334805381];
% T_man = 2*T0; % Manifold propagation time

% %%% L2 - Planar (jupiter-europa) ... Great demonstration of landing (yin)
% T0 = 3.081438015250640;
% X0 = [1.020052317117018;...
%      -0.000000001571777;...
%                       0;...
%      -0.000000000930173;...
%       0.002757023618697;...
%                       0];
% T_man = 4*T0; % Manifold propagation time

% %%% L2 - Clamshell-500 (jupiter-europa)
% T0 = 3.238693317042185;
% X0 = [1.019962098342030;...
%    0.000000000028209;...
%    0.000000000425073;...
%    0.000000000169420;...
%   -0.001111386718711;...
%   -0.017392960199180];
% T_man = 2*T0*1.5; % Manifold propagation time

%%% L2 - Clamshell-10000 (jupiter-europa)
T0 = 3.288177057601299;
X0 = [1.019461532206182;...
   0.000000000123262;...
   0.000000001310867;...
   0.000000000751399;...
  -0.002268518866163;...
  -0.024385906535592];
T_man = 2*T0*1.5; % Manifold propagation time

% %%% L3 - Planar (jupiter-europa) .... SOOO STABLE WTF
% T0 = 6.284553305974939;
% X0 = [-1.000088378334708;...
%        0.002657329394384;...
%                        0;...
%        0.001328060223484;...
%        0.000000220283315;...
%                        0];
% T_man = 2*T0*30; % Manifold propagation time

   

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
[T_PO, X_PO] = ode45(@int_CR3BnSTM_J2, time, X0, options, secondary.MR, primary.R/rNorm, secondary.R_n, primary.J2, 0);

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
    [T_start, X_start] = ode45(@int_CR3BnSTM_J2, time, X0, options, secondary.MR, primary.R/rNorm, secondary.R_n, primary.J2, 0);
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
        [T_int, X] = ode45(@int_CR3BnSTM_J2, time, IC, options, secondary.MR, primary.R/rNorm, secondary.R_n, primary.J2, 0);
        
        %%% Plot stable manifold
        figure(1); hold all
        p1 = plot3(X(:,1),X(:,2),X(:,3),'color',colors.std.grn);
        plot3(X(1,1),X(1,2),X(1,3),'o','markersize',5,'color',colors.std.grn)
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
        [T_int, X] = ode45(@int_CR3BnSTM_J2, time, IC, options, secondary.MR, primary.R/rNorm, secondary.R_n, primary.J2, 0);
        
        %%% Plot unstable manifold
        figure(1); hold all
        p2 = plot3(X(:,1),X(:,2),X(:,3),'r');
        plot3(X(1,1),X(1,2),X(1,3),'ro','markersize',5)
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
        [T_int, X] = ode45(@int_CR3BnSTM_J2, time, IC, options, secondary.MR, primary.R/rNorm, secondary.R_n, primary.J2, 0);
        
        %%% Plot stable manifold
        figure(1); hold all
        p1 = plot3(X(:,1),X(:,2),X(:,3),'color',colors.std.grn);
        plot3(X(1,1),X(1,2),X(1,3),'o','markersize',5,'color',colors.std.grn)
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
        [T_int, X] = ode45(@int_CR3BnSTM_J2, time, IC, options, secondary.MR, primary.R/rNorm, secondary.R_n, primary.J2, 0);
        
        %%% Plot unstable manifold
        figure(1); hold all
        p2 = plot3(X(:,1),X(:,2),X(:,3),'color',colors.std.red);
        plot3(X(1,1),X(1,2),X(1,3),'o','markersize',5,'color',colors.std.red)
    end
    
    
end

%%% Plotting periodic orbit
figure(1)
plot3(X_PO(:,1),X_PO(:,2),X_PO(:,3),'b','linewidth',2)

%%% Plotting bodies
% plotBody2(secondary.R_n,[1-secondary.MR,0,0],secondary.color,colors.std.black,0.5,.5)
plotBody3(secondary.R_n,[1-secondary.MR,0,0],secondary.color,0.5)
view(0,90)
% plotBody2(primary.R/rNorm,[-secondary.MR,0,0],primary.color,colors.std.black,0.5)
% plotBodyTexture3(secondary.R_n,[1-secondary.MR,0,0],secondary.img)
% plotBodyTexture3(primary.R/rNorm,[-secondary.MR,0,0],primary.img)

%%% Plotting lagrange points
plot3(Ls_n(1,1),Ls_n(1,2),Ls_n(1,3),'^','markerfacecolor',colors.std.red,'markeredgecolor',colors.std.black,'markersize',8)
plot3(Ls_n(2,1),Ls_n(2,2),Ls_n(2,3),'^','markerfacecolor',colors.std.red,'markeredgecolor',colors.std.black,'markersize',8)

%%% Plot Options
PlotBoi3('X','Y','Z',16)
[~,legHand] = legend([p1 p2],'Stable','Unstable');
hL=findobj(legHand,'type','line');
set(hL,'linewidth',2);

if plot_PO == 1
    figure(2); hold all
    plot3(X_PO(:,1),X_PO(:,2),X_PO(:,3),'b','linewidth',2)
    PlotBoi3('X','Y','Z',16)
end

[JC_L1]   = JacobiConstantCalculator_J2(secondary.MR,Ls_n(1,:),[0,0,0],primary.R/rNorm,secondary.R_n,primary.J2,0);
[JC_L2]   = JacobiConstantCalculator_J2(secondary.MR,Ls_n(2,:),[0,0,0],primary.R/rNorm,secondary.R_n,primary.J2,0);
[JC_Surf] = JacobiConstantCalculator_J2(secondary.MR,[1-secondary.MR-secondary.R_n,0,0],[0,0,0],primary.R/rNorm,secondary.R_n,primary.J2,0);
[JC_PO]   = JacobiConstantCalculator_J2(secondary.MR,X_PO(1,1:3),X_PO(1,4:6),primary.R/rNorm,secondary.R_n,primary.J2,0);







































