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
run_p1       = 0;
run_p2       = 0;
run_p3_p4_p5 = 1;

%% =======================================================================
%%% Problem 1
% ========================================================================
if run_p1 == 1
% ------------------------------------
%%% Setup
% ------------------------------------
%%% Mass ratio
prms.u = 0.1;

%%% Time vector
time0 = linspace(0,2*pi,1000);

%%% Arbitrary initial conditions
X1_0 = [0.5; 0; 0; 0; 0.1; 0.1];
X2_0 = [0.2; 0.3; 0.2; 0.2; 0; 0.1];
X3_0 = [1.1; 0.1; 0; 0.1; 0; -0.3];
X4_0 = [1.1; 0.1; 0; -0.1; 0; -0.3];
X5_0 = [0.5; 0.5; 0; 0; 0; 0];

%%% Colors
col1 = colors.std.black;
col2 = colors.std.ltred;
col3 = colors.std.cyan;
col4 = colors.std.blue;
col5 = colors.std.ltgrn;

% ------------------------------------
%%% Integration
% ------------------------------------
%%% Choosing ode tolerance
tol = 1e-11;

%%% Setting integrator options
options = odeset('RelTol',tol,'AbsTol',tol);

%%% Integrating states
[~, X1_BCRn] = ode113(@Int_CR3Bn, time0, X1_0, options, prms);
[~, X2_BCRn] = ode113(@Int_CR3Bn, time0, X2_0, options, prms);
[~, X3_BCRn] = ode113(@Int_CR3Bn, time0, X3_0, options, prms);
[~, X4_BCRn] = ode113(@Int_CR3Bn, time0, X4_0, options, prms);
[~, X5_BCRn] = ode113(@Int_CR3Bn, time0, X5_0, options, prms);

% ------------------------------------
%%% Gathering Jacobi constants
% ------------------------------------
JCs_X1 = JacobiConstantCalculator(prms.u,X1_BCRn(:,1:3),X1_BCRn(:,4:6));
JCs_X2 = JacobiConstantCalculator(prms.u,X2_BCRn(:,1:3),X2_BCRn(:,4:6));
JCs_X3 = JacobiConstantCalculator(prms.u,X3_BCRn(:,1:3),X3_BCRn(:,4:6));
JCs_X4 = JacobiConstantCalculator(prms.u,X4_BCRn(:,1:3),X4_BCRn(:,4:6));
JCs_X5 = JacobiConstantCalculator(prms.u,X5_BCRn(:,1:3),X5_BCRn(:,4:6));

% ------------------------------------
%%% Plotting Results
% ------------------------------------
figure; hold all
p1 = plot3(X1_BCRn(:,1),X1_BCRn(:,2),X1_BCRn(:,3),'linewidth',1.5,'color',col1);
p2 = plot3(X2_BCRn(:,1),X2_BCRn(:,2),X2_BCRn(:,3),'linewidth',1.5,'color',col2);
p3 = plot3(X3_BCRn(:,1),X3_BCRn(:,2),X3_BCRn(:,3),'linewidth',1.5,'color',col3);
p4 = plot3(X4_BCRn(:,1),X4_BCRn(:,2),X4_BCRn(:,3),'linewidth',1.5,'color',col4);
p5 = plot3(X5_BCRn(:,1),X5_BCRn(:,2),X5_BCRn(:,3),'linewidth',1.5,'color',col5);

plotBody2(0.03,[-prms.u, 0, 0],colors.std.red,colors.std.black,1,1)
plotBody2(0.01,[1-prms.u, 0, 0],colors.std.red,colors.std.black,1,1)

legend([p1, p2, p3, p4, p5],'X1','X2','X3','X4','X5')
PlotBoi3('$X_n$','$Y_n$','$Z_n$',18,'LaTex')
axis equal



figure; hold all 
p1 = plot(time0, JCs_X1,'linewidth',2,'color',col1);
p2 = plot(time0, JCs_X2,'linewidth',2,'color',col2);
p3 = plot(time0, JCs_X3,'linewidth',2,'color',col3);
p4 = plot(time0, JCs_X4,'linewidth',2,'color',col4);
p5 = plot(time0, JCs_X5,'linewidth',2,'color',col5);

legend([p1, p2, p3, p4, p5],'X1','X2','X3','X4','X5')
PlotBoi2('Time$_n$','Jacobi Constants',18,'LaTex')


figure
subplot(5,1,1)
plot(time0, percentchange(JCs_X1),'linewidth',2,'color',col1)
PlotBoi2('','$JC, X_1$',18,'LaTex')
subplot(5,1,2)
plot(time0, percentchange(JCs_X2),'linewidth',2,'color',col2)
PlotBoi2('','$JC, X_2$',18,'LaTex')
subplot(5,1,3)
plot(time0, percentchange(JCs_X3),'linewidth',2,'color',col3)
PlotBoi2('','$JC, X_3$',18,'LaTex')
subplot(5,1,4)
plot(time0, percentchange(JCs_X4),'linewidth',2,'color',col4)
PlotBoi2('','$JC, X_4$',18,'LaTex')
subplot(5,1,5)
plot(time0, percentchange(JCs_X5),'linewidth',2,'color',col5)
PlotBoi2('Time$_n$','$JC, X_5$',18,'LaTex')
end % run_p1 == 1


%% =======================================================================
%%% Problem 2
% ========================================================================
if run_p2 == 1
% ------------------------------------
%%% Setup
% ------------------------------------
%%% Number of mass ratios to test
nMRs = 100;

%%% Set up logspace of MRs between 1e-6 and 0.5
MRs = logspace(-6,0,nMRs);
MRs(end) = MRs(end)*0.5;

%%% Loop through MRs and store x-locations of equilibrium points
L123s = zeros(10,3);
for kk = 1:nMRs
L123 = EquilibriumPoints(MRs(kk),1:3); % [3x3] of L1, L2, L3 normalized BCR coordinates
L123s(kk,:) = L123(:,1)';
end

%%% Plot x-locations of equilibrium points vs mass ratio
figure; hold all
p1 = plot(L123s(:,1), MRs,'b.','markersize',20);
p2 = plot(L123s(:,2), MRs,'r.','markersize',20);
p3 = plot(L123s(:,3), MRs,'g.','markersize',20);
setLogPlot('y')
legend([p1, p2, p3],'L1','L2','L3','location','southeast')
PlotBoi2('X-Location of Equillibria','Mass Ratio, $\mu$',18,'LaTex')


end % run_p2 == 1

%% =======================================================================
%%% Problems 3, 4, & 5
% ========================================================================
%%% Problem 3
% ========================================================================
if run_p3_p4_p5 == 1
% ------------------------------------
%%% Positions of bodies and L4 / L5
% ------------------------------------
givenValue = 0.5 * (1 - sqrt(23/27));

u1 = givenValue * 0.5;
u2 = givenValue * 2;

L4_std = R3([1,0,0], 60*pi/180);
L5_std = R3([1,0,0], -60*pi/180);

x1_u1 = [-u1, 0, 0];
x2_u1 = [1-u1, 0, 0];

x1_u2 = [-u2, 0, 0];
x2_u2 = [1-u2, 0, 0];

L4_u1 = L4_std + x1_u1;
L5_u1 = L5_std + x1_u1;

L4_u2 = L4_std + x1_u2;
L5_u2 = L5_std + x1_u2;

L123_u1 = EquilibriumPoints(u1,1:3); % [3x3] of L1, L2, L3 normalized BCR coordinates
L123_u2 = EquilibriumPoints(u2,1:3); % [3x3] of L1, L2, L3 normalized BCR coordinates
% ------------------------------------
%%% X0s in vicinity of L4s/L5s
% ------------------------------------
X0s_L4_u1 = [L4_u1' + (rand(3,100)-0.5)*0.01; zeros(3,100)];
X0s_L5_u1 = [L5_u1' + (rand(3,100)-0.5)*0.01; zeros(3,100)];

X0s_L4_u2 = [L4_u2' + (rand(3,100)-0.5)*0.0001; zeros(3,100)];
X0s_L5_u2 = [L5_u2' + (rand(3,100)-0.5)*0.0001; zeros(3,100)];

% ------------------------------------
%%% Propagating states
% ------------------------------------
%%% Time vector
nTime = 300;
time02 = linspace(0,6*pi,nTime);

%%% Choosing ode tolerance
tol = 1e-12;

%%% Setting integrator options
options = odeset('RelTol',tol,'AbsTol',tol);

%%% Setting parameters
prms1.u = u1;
prms2.u = u2;

% ------------------------------------
%%% Integrating and plotting u1 case
% ------------------------------------
figure; hold all
for kk = 1:100
    [~, X1_L4_u1_temp] = ode113(@Int_CR3Bn, time02, X0s_L4_u1(:,kk), options, prms1);
    [~, X1_L5_u1_temp] = ode113(@Int_CR3Bn, time02, X0s_L5_u1(:,kk), options, prms1);
    plot3(X1_L4_u1_temp(:,1),X1_L4_u1_temp(:,2),X1_L4_u1_temp(:,3),'b')
    plot3(X1_L5_u1_temp(:,1),X1_L5_u1_temp(:,2),X1_L5_u1_temp(:,3),'b')
end

plot3(L4_u1(1),L4_u1(2),L4_u1(3),'^','markeredgecolor',colors.std.black,'markerfacecolor',colors.std.ltgrn)
plot3(L5_u1(1),L5_u1(2),L5_u1(3),'^','markeredgecolor',colors.std.black,'markerfacecolor',colors.std.ltgrn)
plot3(L123_u1(1,1),L123_u1(1,2),L123_u1(1,3),'^','markeredgecolor',colors.std.black,'markerfacecolor',colors.std.ltgrn)
plot3(L123_u1(2,1),L123_u1(2,2),L123_u1(2,3),'^','markeredgecolor',colors.std.black,'markerfacecolor',colors.std.ltgrn)
plot3(L123_u1(3,1),L123_u1(3,2),L123_u1(3,3),'^','markeredgecolor',colors.std.black,'markerfacecolor',colors.std.ltgrn)
PlotBoi3('$X_n$','$Y_n$','$Z_n$',18,'LaTex')
axis equal
title(['\mu_1',sprintf(' = %1.4f',u1)])
% ------------------------------------
%%% Integrating and plotting u2 case
% ------------------------------------
figure; hold all
for kk = 1:100
    [~, X1_L4_u2_temp] = ode113(@Int_CR3Bn, time02, X0s_L4_u2(:,kk), options, prms2);
    [~, X1_L5_u2_temp] = ode113(@Int_CR3Bn, time02, X0s_L5_u2(:,kk), options, prms2);
    plot3(X1_L4_u2_temp(:,1),X1_L4_u2_temp(:,2),X1_L4_u2_temp(:,3),'b')
    plot3(X1_L5_u2_temp(:,1),X1_L5_u2_temp(:,2),X1_L5_u2_temp(:,3),'b')
end

plot3(L4_u2(1),L4_u2(2),L4_u2(3),'^','markeredgecolor',colors.std.black,'markerfacecolor',colors.std.ltgrn)
plot3(L5_u2(1),L5_u2(2),L5_u2(3),'^','markeredgecolor',colors.std.black,'markerfacecolor',colors.std.ltgrn)
plot3(L123_u2(1,1),L123_u2(1,2),L123_u2(1,3),'^','markeredgecolor',colors.std.black,'markerfacecolor',colors.std.ltgrn)
plot3(L123_u2(2,1),L123_u2(2,2),L123_u2(2,3),'^','markeredgecolor',colors.std.black,'markerfacecolor',colors.std.ltgrn)
plot3(L123_u2(3,1),L123_u2(3,2),L123_u2(3,3),'^','markeredgecolor',colors.std.black,'markerfacecolor',colors.std.ltgrn)
PlotBoi3('$X_n$','$Y_n$','$Z_n$',18,'LaTex')
axis equal
title(['\mu_2',sprintf(' = %1.4f',u2)])

% ------------------------------------
%%% Grabbing trajectories for mu1 and mu2 at L4 and L5 for problem 5
% ------------------------------------
X_L4_BCR_unstable = X1_L4_u2_temp;
X_L4_BCR_stable   = X1_L4_u1_temp;

X_L5_BCR_unstable = X1_L5_u2_temp;
X_L5_BCR_stable   = X1_L5_u1_temp;

%% =======================================================================
%%% Problem  4
% ========================================================================

% ------------------------------------
%%% Setup
% ------------------------------------
prms.u = 0.01;
L123 = EquilibriumPoints(prms.u,1:3); % [3x3] of L1, L2, L3 normalized BCR coordinates

%%% perturbation from equilibrium point
pert = 0.00001;

% ------------------------------------
%%% X0s on either side of L1 and L2 equilibrium points
% ------------------------------------
X01_L1 = [L123(1,1) + pert; zeros(5,1)];
X02_L1 = [L123(1,1) - pert; zeros(5,1)];

X01_L2 = [L123(2,1) + pert; zeros(5,1)];
X02_L2 = [L123(2,1) - pert; zeros(5,1)];

% ------------------------------------
%%% Integrating forward to find unstable manifolds
% ------------------------------------
%%% Time vector
% nTime = 300;
tPi = 6;
time0_unstable = linspace(0,tPi*pi,nTime);
time0_stable   = linspace(tPi*pi,0,nTime);

%%% Choosing ode tolerance
tol = 1e-12;

%%% Setting integrator options
options = odeset('RelTol',tol,'AbsTol',tol);

%%% Integrating unstable manifolds
[~, X1_L1_BCR_unstable] = ode113(@Int_CR3Bn, time0_unstable, X01_L1, options, prms);
[~, X2_L1_BCR_unstable] = ode113(@Int_CR3Bn, time0_unstable, X02_L1, options, prms);
[~, X1_L2_BCR_unstable] = ode113(@Int_CR3Bn, time0_unstable, X01_L2, options, prms);
[~, X2_L2_BCR_unstable] = ode113(@Int_CR3Bn, time0_unstable, X02_L2, options, prms);

%%% Integrating stable manifolds
[~, X1_L1_BCR_stable] = ode113(@Int_CR3Bn, time0_stable, X01_L1, options, prms);
[~, X2_L1_BCR_stable] = ode113(@Int_CR3Bn, time0_stable, X02_L1, options, prms);
[~, X1_L2_BCR_stable] = ode113(@Int_CR3Bn, time0_stable, X01_L2, options, prms);
[~, X2_L2_BCR_stable] = ode113(@Int_CR3Bn, time0_stable, X02_L2, options, prms);


% ------------------------------------
%%% plotting
% ------------------------------------
figure; hold all

p1 = plot3(X1_L1_BCR_unstable(:,1),X1_L1_BCR_unstable(:,2),X1_L1_BCR_unstable(:,3),'r','linewidth',2);
plot3(X2_L1_BCR_unstable(:,1),X2_L1_BCR_unstable(:,2),X2_L1_BCR_unstable(:,3),'r','linewidth',2)
plot3(X1_L2_BCR_unstable(:,1),X1_L2_BCR_unstable(:,2),X1_L2_BCR_unstable(:,3),'r','linewidth',2)
plot3(X2_L2_BCR_unstable(:,1),X2_L2_BCR_unstable(:,2),X2_L2_BCR_unstable(:,3),'r','linewidth',2)

p2 = plot3(X1_L1_BCR_stable(:,1),X1_L1_BCR_stable(:,2),X1_L1_BCR_stable(:,3),'g','linewidth',2);
plot3(X2_L1_BCR_stable(:,1),X2_L1_BCR_stable(:,2),X2_L1_BCR_stable(:,3),'g','linewidth',2)
plot3(X1_L2_BCR_stable(:,1),X1_L2_BCR_stable(:,2),X1_L2_BCR_stable(:,3),'g','linewidth',2)
plot3(X2_L2_BCR_stable(:,1),X2_L2_BCR_stable(:,2),X2_L2_BCR_stable(:,3),'g','linewidth',2)

p3 = plot3(L123(1,1),L123(1,2),L123(1,3),'^','markeredgecolor',colors.std.black,'markerfacecolor',colors.std.ltgrn,'markersize',8);
plot3(L123(2,1),L123(2,2),L123(2,3),'^','markeredgecolor',colors.std.black,'markerfacecolor',colors.std.ltgrn,'markersize',8)

PlotBoi3('$X_n$','$Y_n$','$Z_n$',18,'LaTex')
axis equal

legend([p1 p2 p3],'Unstable','Stable','L_1 & L_2')
title([sprintf('t = %1d',max(time0_unstable)/pi),'\pi'])

%% =======================================================================
%%% Problem 5
% ========================================================================
% ------------------------------------
%%% Build rotation matrix
% ------------------------------------
%%% Setting up rotation matrix
T = @(t) [cos(t), -sin(t), 0;...
          sin(t), cos(t),  0;...
          0,      0,       1];
      
      
% ------------------------------------
%%% Rotate frames into inertial
% ------------------------------------
%%% Preallocating
X1_L1_BCI_unstable = zeros(size(X1_L1_BCR_unstable));
X2_L1_BCI_unstable = zeros(size(X2_L1_BCR_unstable));
X1_L2_BCI_unstable = zeros(size(X1_L2_BCR_unstable));
X2_L2_BCI_unstable = zeros(size(X2_L2_BCR_unstable));

X1_L1_BCI_stable = zeros(size(X1_L1_BCR_stable));
X2_L1_BCI_stable = zeros(size(X2_L1_BCR_stable));
X1_L2_BCI_stable = zeros(size(X1_L2_BCR_stable));
X2_L2_BCI_stable = zeros(size(X2_L2_BCR_stable));


X_L4_BCI_unstable = zeros(size(X_L4_BCR_unstable));
X_L4_BCI_stable   = zeros(size(X_L4_BCR_stable));
X_L5_BCI_unstable = zeros(size(X_L5_BCR_unstable));
X_L5_BCI_stable   = zeros(size(X_L5_BCR_stable));

%%% Looping through times to rotate
for kk = 1:nTime
    %%% Grab current times
    t_unstable_kk = time0_unstable(kk);
    t_stable_kk   = time0_stable(kk);
    
    %%% Build current rotation matrices
    T_unstable = T(t_unstable_kk);
    T_stable   = T(t_stable_kk);
    
    %%% Rotate states into inertial frame
    X1_L1_BCI_unstable(kk,:) = [T_unstable * X1_L1_BCR_unstable(kk,1:3)';...
                                T_unstable * (X1_L1_BCR_unstable(kk,4:6) + [-X1_L1_BCR_unstable(kk,2), X1_L1_BCR_unstable(kk,1), 0])'];
    X2_L1_BCI_unstable(kk,:) = [T_unstable * X2_L1_BCR_unstable(kk,1:3)';...
                                T_unstable * (X2_L1_BCR_unstable(kk,4:6) + [-X2_L1_BCR_unstable(kk,2), X2_L1_BCR_unstable(kk,1), 0])'];
    X1_L2_BCI_unstable(kk,:) = [T_unstable * X1_L2_BCR_unstable(kk,1:3)';...
                                T_unstable * (X1_L2_BCR_unstable(kk,4:6) + [-X1_L2_BCR_unstable(kk,2), X1_L2_BCR_unstable(kk,1), 0])'];
    X2_L2_BCI_unstable(kk,:) = [T_unstable * X2_L2_BCR_unstable(kk,1:3)';...
                                T_unstable * (X2_L2_BCR_unstable(kk,4:6) + [-X2_L2_BCR_unstable(kk,2), X2_L2_BCR_unstable(kk,1), 0])'];
                            
    X1_L1_BCI_stable(kk,:) = [T_stable * X1_L1_BCR_stable(kk,1:3)';...
                                T_stable * (X1_L1_BCR_stable(kk,4:6) + [-X1_L1_BCR_stable(kk,2), X1_L1_BCR_stable(kk,1), 0])'];
    X2_L1_BCI_stable(kk,:) = [T_stable * X2_L1_BCR_stable(kk,1:3)';...
                                T_stable * (X2_L1_BCR_stable(kk,4:6) + [-X2_L1_BCR_stable(kk,2), X2_L1_BCR_stable(kk,1), 0])'];
    X1_L2_BCI_stable(kk,:) = [T_stable * X1_L2_BCR_stable(kk,1:3)';...
                                T_stable * (X1_L2_BCR_stable(kk,4:6) + [-X1_L2_BCR_stable(kk,2), X1_L2_BCR_stable(kk,1), 0])'];
    X2_L2_BCI_stable(kk,:) = [T_stable * X2_L2_BCR_stable(kk,1:3)';...
                                T_stable * (X2_L2_BCR_stable(kk,4:6) + [-X2_L2_BCR_stable(kk,2), X2_L2_BCR_stable(kk,1), 0])'];
    
    X_L4_BCI_unstable(kk,:) = [T_stable * X_L4_BCR_unstable(kk,1:3)';...
                                T_stable * (X_L4_BCR_unstable(kk,4:6) + [-X_L4_BCR_unstable(kk,2), X_L4_BCR_unstable(kk,1), 0])'];
    X_L4_BCI_stable(kk,:) = [T_stable * X_L4_BCR_stable(kk,1:3)';...
                                T_stable * (X_L4_BCR_stable(kk,4:6) + [-X_L4_BCR_stable(kk,2), X_L4_BCR_stable(kk,1), 0])'];
    X_L5_BCI_unstable(kk,:) = [T_stable * X_L5_BCR_unstable(kk,1:3)';...
                                T_stable * (X_L5_BCR_unstable(kk,4:6) + [-X_L5_BCR_unstable(kk,2), X_L5_BCR_unstable(kk,1), 0])'];
    X_L5_BCI_stable(kk,:) = [T_stable * X_L5_BCR_stable(kk,1:3)';...
                                T_stable * (X_L5_BCR_stable(kk,4:6) + [-X_L5_BCR_stable(kk,2), X_L5_BCR_stable(kk,1), 0])'];
end

% ------------------------------------
%%% Plotting Inertial States
% ------------------------------------
figure; hold all

p1 = plot3(X1_L1_BCI_unstable(:,1),X1_L1_BCI_unstable(:,2),X1_L1_BCI_unstable(:,3),'r','linewidth',2);
plot3(X2_L1_BCI_unstable(:,1),X2_L1_BCI_unstable(:,2),X2_L1_BCI_unstable(:,3),'r','linewidth',2)
plot3(X1_L2_BCI_unstable(:,1),X1_L2_BCI_unstable(:,2),X1_L2_BCI_unstable(:,3),'r','linewidth',2)
plot3(X2_L2_BCI_unstable(:,1),X2_L2_BCI_unstable(:,2),X2_L2_BCI_unstable(:,3),'r','linewidth',2)

p2 = plot3(X1_L1_BCI_stable(:,1),X1_L1_BCI_stable(:,2),X1_L1_BCI_stable(:,3),'g','linewidth',2);
plot3(X2_L1_BCI_stable(:,1),X2_L1_BCI_stable(:,2),X2_L1_BCI_stable(:,3),'g','linewidth',2)
plot3(X1_L2_BCI_stable(:,1),X1_L2_BCI_stable(:,2),X1_L2_BCI_stable(:,3),'g','linewidth',2)
plot3(X2_L2_BCI_stable(:,1),X2_L2_BCI_stable(:,2),X2_L2_BCI_stable(:,3),'g','linewidth',2)

PlotBoi3('$X_n$','$Y_n$','$Z_n$',18,'LaTex')
axis equal

legend([p1 p2],'Unstable','Stable')
title([sprintf('t = %1d',max(time0_unstable)/pi),'\pi'])

% ------------------------------------
%%% Calculating orbital elements
% ------------------------------------
%%% Preallocating
aeiwt_X1_L1_unstable = zeros(size(X1_L1_BCI_unstable,1),4);
aeiwt_X2_L1_unstable = zeros(size(X2_L1_BCI_unstable,1),4);
aeiwt_X1_L2_unstable = zeros(size(X1_L2_BCI_unstable,1),4);
aeiwt_X2_L2_unstable = zeros(size(X2_L2_BCI_unstable,1),4);

aeiwt_X1_L1_stable = zeros(size(X1_L1_BCI_stable,1),4);
aeiwt_X2_L1_stable = zeros(size(X2_L1_BCI_stable,1),4);
aeiwt_X1_L2_stable = zeros(size(X1_L2_BCI_stable,1),4);
aeiwt_X2_L2_stable = zeros(size(X2_L2_BCI_stable,1),4);

aeiwt_X_L4_unstable = zeros(size(X_L4_BCI_unstable,1),4);
aeiwt_X_L4_stable   = zeros(size(X_L4_BCI_stable,1),4);
aeiwt_X_L5_unstable = zeros(size(X_L5_BCI_unstable,1),4);
aeiwt_X_L5_stable   = zeros(size(X_L5_BCI_stable,1),4);

for kk = 1:nTime
    %%% Calculating orbital elements
    [a_X1L1un,e_X1L1un,i_X1L1un,raan_X1L1un,w_X1L1un,~] = RV2COE(X1_L1_BCI_unstable(kk,1:3),X1_L1_BCI_unstable(kk,4:6),1);
    [a_X2L1un,e_X2L1un,i_X2L1un,raan_X2L1un,w_X2L1un,~] = RV2COE(X2_L1_BCI_unstable(kk,1:3),X2_L1_BCI_unstable(kk,4:6),1);
    [a_X1L2un,e_X1L2un,i_X1L2un,raan_X1L2un,w_X1L2un,~] = RV2COE(X1_L2_BCI_unstable(kk,1:3),X1_L2_BCI_unstable(kk,4:6),1);
    [a_X2L2un,e_X2L2un,i_X2L2un,raan_X2L2un,w_X2L2un,~] = RV2COE(X2_L2_BCI_unstable(kk,1:3),X2_L2_BCI_unstable(kk,4:6),1);
    
    [a_X1L1st,e_X1L1st,i_X1L1st,raan_X1L1st,w_X1L1st,~] = RV2COE(X1_L1_BCI_stable(kk,1:3),X1_L1_BCI_stable(kk,4:6),1);
    [a_X2L1st,e_X2L1st,i_X2L1st,raan_X2L1st,w_X2L1st,~] = RV2COE(X2_L1_BCI_stable(kk,1:3),X2_L1_BCI_stable(kk,4:6),1);
    [a_X1L2st,e_X1L2st,i_X1L2st,raan_X1L2st,w_X1L2st,~] = RV2COE(X1_L2_BCI_stable(kk,1:3),X1_L2_BCI_stable(kk,4:6),1);
    [a_X2L2st,e_X2L2st,i_X2L2st,raan_X2L2st,w_X2L2st,~] = RV2COE(X2_L2_BCI_stable(kk,1:3),X2_L2_BCI_stable(kk,4:6),1);
    
    [a_L4un,e_L4un,i_L4un,raan_L4un,w_L4un,~] = RV2COE(X_L4_BCI_unstable(kk,1:3),X_L4_BCI_unstable(kk,4:6),1);
    [a_L4st,e_L4st,i_L4st,raan_L4st,w_L4st,~] = RV2COE(X_L4_BCI_stable(kk,1:3),X_L4_BCI_stable(kk,4:6),1);
    [a_L5un,e_L5un,i_L5un,raan_L5un,w_L5un,~] = RV2COE(X_L5_BCI_unstable(kk,1:3),X_L5_BCI_unstable(kk,4:6),1);
    [a_L5st,e_L5st,i_L5st,raan_L5st,w_L5st,~] = RV2COE(X_L5_BCI_stable(kk,1:3),X_L5_BCI_stable(kk,4:6),1);
    
    %%% Storing orbital elements
    aeiwt_X1_L1_unstable(kk,:) = [a_X1L1un,e_X1L1un,i_X1L1un,w_X1L1un + raan_X1L1un];
    aeiwt_X2_L1_unstable(kk,:) = [a_X2L1un,e_X2L1un,i_X2L1un,w_X2L1un + raan_X2L1un];
    aeiwt_X1_L2_unstable(kk,:) = [a_X1L2un,e_X1L2un,i_X1L2un,w_X1L2un + raan_X1L2un];
    aeiwt_X2_L2_unstable(kk,:) = [a_X2L2un,e_X2L2un,i_X2L2un,w_X2L2un + raan_X2L2un];

    aeiwt_X1_L1_stable(kk,:) = [a_X1L1st,e_X1L1st,i_X1L1st,w_X1L1st + raan_X1L1st];
    aeiwt_X2_L1_stable(kk,:) = [a_X2L1st,e_X2L1st,i_X2L1st,w_X2L1st + raan_X2L1st];
    aeiwt_X1_L2_stable(kk,:) = [a_X1L2st,e_X1L2st,i_X1L2st,w_X1L2st + raan_X1L2st];
    aeiwt_X2_L2_stable(kk,:) = [a_X2L2st,e_X2L2st,i_X2L2st,w_X2L2st + raan_X2L2st];
    
    aeiwt_X_L4_unstable(kk,:) = [a_L4un,e_L4un,i_L4un,w_L4un + raan_L4un];
    aeiwt_X_L4_stable(kk,:)   = [a_L4st,e_L4st,i_L4st,w_L4st + raan_L4st];
    aeiwt_X_L5_unstable(kk,:) = [a_L5un,e_L5un,i_L5un,w_L5un + raan_L5un];
    aeiwt_X_L5_stable(kk,:)   = [a_L5st,e_L5st,i_L5st,w_L5st + raan_L5st];

end


% ------------------------------------
%%% Plotting Orbital Elements
% ------------------------------------
close all
figure; hold all
plot_aeiwt(time0_unstable,aeiwt_X2_L1_unstable,'L1, Unstable Manifold')

figure; hold all
plot_aeiwt(time0_stable,aeiwt_X2_L1_stable, 'L1, Stable Manifold')

figure; hold all
plot_aeiwt(time0_unstable,aeiwt_X2_L2_unstable, 'L2, Unstable Manifold')

figure; hold all
plot_aeiwt(time0_stable,aeiwt_X2_L2_stable, 'L2, Stable Manifold')

% figure; hold all
% plot_aeiwt(time02,aeiwt_X_L4_unstable, 'L4, Unstable \mu')
% 
% figure; hold all
% plot_aeiwt(time02,aeiwt_X_L4_stable, 'L4, Stable \mu')
% 
% figure; hold all
% plot_aeiwt(time02,aeiwt_X_L5_unstable, 'L5, Unstable \mu')
% 
% figure; hold all
% plot_aeiwt(time02,aeiwt_X_L5_stable, 'L5, Stable \mu')

end % run_p3_p4_p5 == 1












function plot_aeiwt(time_vec, aeiwt_mat,title_string)
subplot(4,1,1); hold all
title(title_string)
plot(time_vec, aeiwt_mat(:,1),'linewidth',2)
PlotBoi2('','$a$',18,'LaTex')

subplot(4,1,2); hold all
plot(time_vec, aeiwt_mat(:,2),'linewidth',2)
PlotBoi2('','$e$',18,'LaTex')

subplot(4,1,3); hold all
plot(time_vec, aeiwt_mat(:,3),'linewidth',2)
PlotBoi2('','$i, ^\circ$',18,'LaTex')

subplot(4,1,4); hold all
plot(time_vec, aeiwt_mat(:,4),'linewidth',2)
PlotBoi2('$Time_n$','$\widetilde{\omega}$, $^\circ$',18,'LaTex')
end









