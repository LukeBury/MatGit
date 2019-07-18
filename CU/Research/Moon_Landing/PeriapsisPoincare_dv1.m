clear
clc
close all
mbinPath = '~/CU_Google_Drive/Documents/MatGit/mbin';
addpath(genpath(mbinPath))
tic
% ========================================================================
%%% Run Switches
% ========================================================================
run_apsesAnalysis = 1;

% ========================================================================
%%% Importing Data
% ========================================================================
%%% General data on solar system bodies
bodies = getBodyData(mbinPath);

%%% Color options/schemes
colors = get_colors();

% ========================================================================
%%% Setup
% ========================================================================
% -------------------------------------------------
% Choose 3B system and Lagrange point
% -------------------------------------------------
%%% 3B system
primary   = bodies.jupiter;
secondary = bodies.europa;

%%% Collinear lagrange point of interest
Lpoint = 2;  % 1 or 2

%%% Normalizing constants
rNorm = secondary.a;         % n <-> km
tNorm = 1/secondary.meanMot; % n <-> sec
vNorm = rNorm / tNorm;       % n <-> km/sec

prms.u = secondary.MR;
prms.R2_n = secondary.R_n;
% -------------------------------------------------
% Equilibrium Points and Jacobi Constants
% -------------------------------------------------
%%% Acquire Collinear Lagrange points
L123 = EquilibriumPoints(secondary.MR,1:3); % [3x3] of L1, L2, L3 normalized BCR coordinates

%%% Jacobi constant of Lagrange point
[JC_Lp] = JacobiConstantCalculator(secondary.MR,L123(Lpoint,:),[0,0,0]);

% -------------------------------------------------
% Specifying set of initial conditions
% -------------------------------------------------
%%% Difference in JC value between Lagrange point and s/c initial state.
%%% Basically controls where we want our ZV curves
dJC_Lp = 0.000017; 

%%% s/c starting JC (JC_scDesired) is lower than JC_Lp because the dJC_Lp is accounted for
%%% in velocity
JC_scDesired = JC_Lp-dJC_Lp;

%%% Initial position
r0_n = [L123(Lpoint,1), 0.9*secondary.R_n, 0];

%%% With JC0 defined, starting velocity is a function of position. So first
%%% we must calculate the JC of the stationary starting position
JC_initialPos = JacobiConstantCalculator(secondary.MR,r0_n,[0,0,0]);

%%% Starting velocity is found from difference between s/c JC (JC_scDesired) and the
%%% JC of the stationary starting position (JC_initialPos)
dJC_forInitialVelocity = JC_initialPos - JC_scDesired;

%%% Find necessary velocity magnitude
if dJC_forInitialVelocity < 0
    warning('Spacecraft starting in a forbidden zone')
    return
elseif dJC_forInitialVelocity > 0
    v0i = sqrt(abs(dJC_forInitialVelocity));

    % Find direction for velocity
    vHat = R3([0, -1, 0], (.75)*2*pi);

    % Create initial velocity vector
    v0_n = vHat .* v0i;
end

%%% Store initial conditions
X0_n = [r0_n, v0_n];

% -------------------------------------------------
% Integration options
% -------------------------------------------------
%%% Setting time vector
t_i = 0; % sec
t_f = 2*pi;
dt = t_f/10000;
time0_n = t_i:dt:t_f;

%%% Choosing ode45 tolerance
tol = 1e-10;

%%% Setting integrator options
options = odeset('RelTol',tol,'AbsTol',tol);
options_Apsis = odeset('Events',@event_Apsis,'RelTol',tol,'AbsTol',tol);
options_Impact = odeset('Events',@event_Impact_CR3Bn,'RelTol',tol,'AbsTol',tol);

% -------------------------------------------------
% dV Specifications
% -------------------------------------------------
% t_dV = 0.395*pi;
% % t_dV = 0.1*pi;
dVMag_metersPerSecond = 40; % m/s
dVMag_n = dVMag_metersPerSecond/1000/vNorm;
% % dV_headingAngle = (.5)*2*pi; % radians from direction of motion

% -------------------------------------------------
% Finding bounds on possible post-dV JCs and creating color scheme
% -------------------------------------------------
JC_min = JC_scDesired - dVMag_n^2;
JC_max = JC_scDesired + dVMag_n^2;

[ colorMatrix ] = colorScale([colors.std.mag; colors.std.cyan],6 );


% -------------------------------------------------
% Loop Specifications
% -------------------------------------------------
for t_dV = [0.05:0.05:0.41].*pi
for dV_headingAngle = 0:(2*pi/8):(2*pi)
% ========================================================================
%%% Looping through integrations and storing periapses
% ========================================================================
% -------------------------------------------------
% Preallocation
% -------------------------------------------------
X_apses_SCR_n = [];
time_apses = [];

% -------------------------------------------------
% Integration
% -------------------------------------------------
%%% Propagating Pre-Burn States
[time_n_predV, X_BCR_n_predV] = ode45(@Int_CR3Bn, [t_i:dt:t_dV], X0_n, options, prms);

%%% Adding tangent dV to states
% Finding vHat
vHat = X_BCR_n_predV(end,4:6)./norm(X_BCR_n_predV(end,4:6));
vHat = R3(vHat,dV_headingAngle);

% Determining dV
dV_n = vHat .* dVMag_n;

% Creating post-burn initial states
newVelocity = X_BCR_n_predV(end,4:6) + dV_n;
X0_postdV_n = [X_BCR_n_predV(end,1:3), newVelocity];

%%% Propagating the post-burn States with Apsis event
if run_apsesAnalysis == 1
    [time_n_Apsis, X_BCR_n_Apsis, time_eventApsis, X_eventApsis, index_eventApsis] = ode45(@Int_CR3Bn, [t_dV:dt:t_f], X0_postdV_n,...
                                                           options_Apsis, prms);
    %%% Storing periapses and apoapses
    X_apses_SCR_n = [X_apses_SCR_n; X_eventApsis];
    time_apses = [time_apses; time_eventApsis];
end

 %%% Propagating the post-burn States with Impact event
[time_n_Impact, X_BCR_n_Impact, time_eventImpact, X_eventImpact, index_eventImpact] = ode45(@Int_CR3Bn, [t_dV:dt:t_f], X0_postdV_n,...
                                                       options_Impact, prms);
  

% ========================================================================
%%% Post-Integration Work
% ========================================================================
% -------------------------------------------------
% Calculating Lat/Lon impact locations
% -------------------------------------------------
impactLatLons = [];
for kk = 1:size(X_eventImpact,1)
    %%% Creating SCR position
    rImpact_SCR = X_eventImpact(kk,1:3) - [1-secondary.MR, 0, 0];
    
    %%% Finding lat/lon of impact site
    [lat, lon] = ECEF2latlon(rImpact_SCR,'degrees');
    
    %%% Storing lat/lon of impact site
    impactLatLons = [impactLatLons; lat, lon];
end

% -------------------------------------------------
% Calculating impact velocities and angles
% -------------------------------------------------
impactVelocities = zeros(size(X_eventImpact,1));
impactAngles     = zeros(size(X_eventImpact,1));

for kk = 1:size(X_eventImpact,1)
    impactVelocities(kk) = norm(X_eventImpact(kk,4:6));
    
    rImpact_SCR_n = X_eventImpact(kk,1:3) - [1-secondary.MR,0,0];
    rImpact_SCR_n = rImpact_SCR_n./norm(rImpact_SCR_n);
    
    vImpact_n = X_eventImpact(kk,4:6)./norm(X_eventImpact(kk,4:6));
    
    A = R3(rImpact_SCR_n,pi/2);
    B = vImpact_n;
    impactAngles(kk) = acos(dot(A,B)/(norm(A)*norm(B)));
    if impactAngles(kk) > pi/2
        impactAngles(kk) = pi/2 - impactAngles(kk);
    end
end

% ========================================================================
%%% Plots
% ========================================================================
% -------------------------------------------------
% figure(1), subplot(2,3,[1 2]) - system, apses, and impacts
% -------------------------------------------------
figure(1); set(gcf,'Position',[266 54 986 678]); subplot(2,3,[1 2 3]); hold all
PlotBoi3('x$_n$','y$_n$','z$_n$',16,'Latex')

%%% Creating ZV contours
% Contour bounds
xCont_min = L123(1,1)-2*secondary.R_n;
xCont_max = L123(2,1)+2*secondary.R_n;
yCont_min = -6*secondary.R_n;
yCont_max = 6*secondary.R_n;

% Creating x-y grid
xs = linspace(xCont_min,xCont_max,600);
ys = linspace(yCont_min,yCont_max,100);
[X, Y] = meshgrid(xs,ys);
clear xs ys

% Calculating JCs across grid
Z = zeros(size(X));
for xk = 1:size(X,1)
    for yk = 1:size(X,2)
        %%% Zero-Velocity Curve
        zv = JacobiConstantCalculator(secondary.MR,[X(xk,yk), Y(xk,yk), 0] ,[0, 0, 0]);
        Z(xk,yk) = zv;
    end
end

%%% Plotting Lagrange Points
plot3(L123(1:2,1), L123(1:2,2), L123(1:2,3), '^', 'markeredgecolor', colors.std.black,...
                                       'markerfacecolor',colors.std.mag,'markersize',10);

%%% Plotting pre-burn contours
[C1,h1] = contour(X,Y,Z,[JC_scDesired, JC_scDesired],'color',colors.sch.d3_1(3,:),'linewidth',1.5);
% colors.sch.r9(9,:)

%%% Plotting post-burn contours
JCf = JacobiConstantCalculator(secondary.MR,[X_BCR_n_Impact(1,1),X_BCR_n_Impact(1,2),X_BCR_n_Impact(1,3)],[X_BCR_n_Impact(1,4),X_BCR_n_Impact(1,5),X_BCR_n_Impact(1,6)]);
[C2,h2] = contour(X,Y,Z,[JCf, JCf],'color',colors.std.black,'linewidth',1);

%%% Plotting pre-burn trajectory
plot3(X_BCR_n_predV(:,1),X_BCR_n_predV(:,2),X_BCR_n_predV(:,3),'linewidth',2,'color',colors.sch.d3_1(3,:))


%%% Plotting dV vector
qS = 1.5; % quiverScale
quiver3(X_BCR_n_predV(end,1),X_BCR_n_predV(end,2),X_BCR_n_predV(end,3),...
    dV_n(1)*qS,dV_n(2)*qS,dV_n(3)*qS,'color',colors.std.black)

%%% Plotting dV marker
[dvMarker] = plot3(X_BCR_n_predV(end,1),X_BCR_n_predV(end,2),X_BCR_n_predV(end,3),'*',...
    'markersize',8,'color',colors.std.black);

%%% Plotting post-burn IMPACT trajectory
plot3(X_BCR_n_Impact(:,1),X_BCR_n_Impact(:,2),X_BCR_n_Impact(:,3),'linewidth',1,'color',colors.sch.d3_1(1,:))

% %%% Plotting post-burn FULL trajectory
% plot3(X_BCR_n_Apsis(:,1),X_BCR_n_Apsis(:,2),X_BCR_n_Apsis(:,3),'linewidth',1,'color',colors.sch.d3_1(1,:))


%%% Plotting apses
if run_apsesAnalysis == 1
    if isempty(X_apses_SCR_n) == 0
        plot3(X_apses_SCR_n(:,1),X_apses_SCR_n(:,2),X_apses_SCR_n(:,3),'o','markersize',5,...
                    'markeredgecolor',colors.std.black,'markerfacecolor',colors.std.ltgrn,'linewidth',2)
    end
end

%%% Plotting secondary
plotBody2(secondary.R_n,[1-secondary.MR,0,0],[1,1,1],colors.std.black,1.5,0)

%%% Plotting impact locations
if isempty(X_eventImpact) == 0
    plot3(X_eventImpact(:,1),X_eventImpact(:,2),X_eventImpact(:,3),'o','markersize',8,...
                    'markeredgecolor',colors.std.black,'markerfacecolor',colors.sch.d3_1(1,:),'linewidth',2)
end

%%% Setting axes limits
xlim([xCont_min xCont_max])
ylim([yCont_min yCont_max])

%%% legend
legend([dvMarker],{'$\Delta V$'},'Interpreter','latex')


% -------------------------------------------------
% figure(1), subplot(2,3,4) - impact longitude vs impact speed
% -------------------------------------------------
subplot(2,3,4); hold all
xlim([-180 180])
PlotBoi2('Longitude, $^\circ$','Velocity Magnitude, $m/s$',12,'Latex')

%%% If there was an impact
if isempty(X_eventImpact) == 0
    %%% Plotting impact longitude vs speed
    plot(impactLatLons(:,2),impactVelocities.*(vNorm*1000),'o','markersize',8,...
                    'markeredgecolor',colors.std.black,'markerfacecolor',colors.sch.d3_1(1,:),'linewidth',2)
    
    %%% Plotting anti-primary line (0-deg longitude)
    [apLine] = plot([0 0],[0 max(abs(impactVelocities))].*(vNorm*1000*1.1),'color',colors.std.blue,'linewidth',1.5);

    %%% legend
    legend([apLine],{'Anti-Primary'},'Interpreter','latex','location','best')
end

% -------------------------------------------------
% figure(1), subplot(2,3,5) - impact longitude vs impact time
% -------------------------------------------------
subplot(2,3,5); hold all
xlim([-180 180])
PlotBoi2('Longitude, $^\circ$','Time to Impact, $\pi$',12,'Latex')

%%% Plotting impact angle vs impact time
if isempty(X_eventImpact) == 0
    %%% Plotting impact longitude vs speed
    plot(impactLatLons(:,2),time_eventImpact,'o','markersize',8,...
                    'markeredgecolor',colors.std.black,'markerfacecolor',colors.sch.d3_1(1,:),'linewidth',2)
    
    %%% Plotting anti-primary line (0-deg longitude)
    [apLine] = plot([0 0],[0 max(abs(time_eventImpact))].*(1.1),'color',colors.std.blue,'linewidth',1.5);

    %%% legend
    legend([apLine],{'Anti-Primary'},'Interpreter','latex','location','best')
end

% -------------------------------------------------
% figure(1), subplot(2,3,6) - impact angle vs impact speed
% -------------------------------------------------
subplot(2,3,6); hold all
PlotBoi2('Impact Angle, $^\circ$', 'Velocity Magnitude, $m/s$',12,'Latex')

%%% Plotting impact angle vs speed
if isempty(X_eventImpact) == 0
    plot(impactAngles.*(180/pi),impactVelocities.*(vNorm*1000),'o','markersize',8,...
                    'markeredgecolor',colors.std.black,'markerfacecolor',colors.sch.d3_1(1,:),'linewidth',2)
end

end
end
toc
















