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
u = 1;

% ------------------------------------
%%% Initial conditions
% ------------------------------------
r0 = [1; 0 ; 0]; % 
v0 = [0; 1; 0]; % km/s
X0 = [r0; v0];

% ------------------------------------
%%% Perturbing Force
% ------------------------------------
aPerts = [0.001; 0.01; 0.1; 1];

% ------------------------------------
%%% Initial elements
% ------------------------------------
[a0,e0,i0,raan0,w0,ta0] = RV2COE(r0,v0,u);
Tp0 = 2*pi*sqrt((a0^3)/u);

% ------------------------------------
%%% Integrator prep
% ------------------------------------
%%% Time info
t0 = 0;
tf = 10*Tp0;
n_t = 10000;
time0 = linspace(t0, tf, n_t);

time_orbits = time0./Tp0;

%%% ODE Tolerance
tol = 1e-13;

%%% Setting integrator options
options = odeset('RelTol',tol,'AbsTol',tol);

% ========================================================================
%%% Integration and Analyzing
% ========================================================================
% ------------------------------------
%%% Integrating
% ------------------------------------
%%% Integrating each case
[~, X_1] = ode113(@Int_2BI_constPert, time0, X0, options, u, [aPerts(1); 0; 0]);
[~, X_2] = ode113(@Int_2BI_constPert, time0, X0, options, u, [aPerts(2); 0; 0]);
[~, X_3] = ode113(@Int_2BI_constPert, time0, X0, options, u, [aPerts(3); 0; 0]);
[~, X_4] = ode113(@Int_2BI_constPert, time0, X0, options, u, [aPerts(4); 0; 0]);

% ------------------------------------
%%% Calculating OEs
% ------------------------------------
%%% Preallocating
as    = zeros(n_t,4);
es    = zeros(n_t,4);
is    = zeros(n_t,4);
raans = zeros(n_t,4);
ws    = zeros(n_t,4);
tas   = zeros(n_t,4);

for kk = 1:n_t
    %%% Calculating OEs
    [a1,e1,i1,raan1,w1,ta1] = RV2COE(X_1(kk,1:3),X_1(kk,4:6),u);
    [a2,e2,i2,raan2,w2,ta2] = RV2COE(X_2(kk,1:3),X_2(kk,4:6),u);
    [a3,e3,i3,raan3,w3,ta3] = RV2COE(X_3(kk,1:3),X_3(kk,4:6),u);
    [a4,e4,i4,raan4,w4,ta4] = RV2COE(X_4(kk,1:3),X_4(kk,4:6),u);
    
    %%% Storing OEs
    as(kk,:)    = [   a1,    a2,    a3,    a4];
    es(kk,:)    = [   e1,    e2,    e3,    e4];
    is(kk,:)    = [   i1,    i2,    i3,    i4];
    raans(kk,:) = [raan1, raan2, raan3, raan4];
    ws(kk,:)    = [   w1,    w2,    w3,    w4];
    tas(kk,:)   = [  ta1,   ta2,   ta3,   ta4];
end




% ========================================================================
%%% Plotting
% ========================================================================
% ------------------------------------
%%% Plotting full systems
% ------------------------------------
figure; hold all
p1 = plot(X_1(:,1),X_1(:,2),'b','linewidth',1);
p2 = plot(0,0,'k.','markersize',30);
PlotBoi2('$X$','$Y$',18,'LaTex')
axis equal
legend([p1 p2],'Perturbed Orbit','Central Body')
title(sprintf('a_s = %1.3f',aPerts(1)))

figure; hold all
p1 = plot(X_2(:,1),X_2(:,2),'b','linewidth',1);
p2 = plot(0,0,'k.','markersize',30);
PlotBoi2('$X$','$Y$',18,'LaTex')
axis equal
legend([p1 p2],'Perturbed Orbit','Central Body')
title(sprintf('a_s = %1.3f',aPerts(2)))

figure; hold all
p1 = plot(X_3(:,1),X_3(:,2),'b','linewidth',1);
p2 = plot(0,0,'k.','markersize',30);
PlotBoi2('$X$','$Y$',18,'LaTex')
axis equal
legend([p1 p2],'Perturbed Orbit','Central Body')
title(sprintf('a_s = %1.3f',aPerts(3)))

figure; hold all
p1 = plot(X_4(:,1),X_4(:,2),'b','linewidth',1);
p2 = plot(0,0,'k.','markersize',30);
PlotBoi2('$X$','$Y$',18,'LaTex')
axis equal
legend([p1 p2],'Perturbed Orbit','Central Body')
title(sprintf('a_s = %1.3f',aPerts(4)))
% ------------------------------------
%%% Plotting OEs
% ------------------------------------
figure
subplot(2,1,1); hold all
title(sprintf('a_s = %1.3f',aPerts(1)))
plot(time_orbits,as(:,1),'b','linewidth',2)
PlotBoi2('','a',18,'LaTex')
subplot(2,1,2); hold all
plot(time_orbits,es(:,1),'b','linewidth',2)
PlotBoi2('Tp','e',18,'LaTex')
% subplot(3,1,3); hold all
% plot(time0,is(:,1),'b')
% PlotBoi2('','i',18,'LaTex')
% subplot(3,2,2); hold all
% plot(time0,raans(:,1),'b')
% subplot(3,2,4); hold all
% plot(time0,ws(:,1),'b')
% subplot(3,2,6); hold all
% plot(time0,tas(:,1),'b')


figure
subplot(2,1,1); hold all
title(sprintf('a_s = %1.3f',aPerts(2)))
plot(time_orbits,as(:,2),'b','linewidth',2)
PlotBoi2('','a',18,'LaTex')
subplot(2,1,2); hold all
plot(time_orbits,es(:,2),'b','linewidth',2)
PlotBoi2('Tp','e',18,'LaTex')


figure
subplot(2,1,1); hold all
title(sprintf('a_s = %1.3f',aPerts(3)))
plot(time_orbits,as(:,3),'b','linewidth',2)
PlotBoi2('','a',18,'LaTex')
subplot(2,1,2); hold all
plot(time_orbits,es(:,3),'b','linewidth',2)
PlotBoi2('Tp','e',18,'LaTex')


figure
subplot(2,1,1); hold all
title(sprintf('a_s = %1.3f',aPerts(4)))
plot(time_orbits,as(:,4),'b','linewidth',2)
PlotBoi2('','a',18,'LaTex')
subplot(2,1,2); hold all
plot(time_orbits,es(:,4),'b','linewidth',2)
PlotBoi2('Tp','e',18,'LaTex')











