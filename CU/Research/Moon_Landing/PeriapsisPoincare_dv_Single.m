clear
clc
close all
addpath(genpath('/Users/CU_Google_Drive/lukebury/Documents/MATLAB/mbin'))
tic
% ========================================================================
%%% Run/Plot Switches
% ========================================================================
run_apsesAnalysis = 0;

plot_referenceTrajectory  = 1;
plot_postBurnTrajectories = 0;
plot_burnQuivers          = 1;

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
% -------------------------------------------------
% Equilibrium Points and Jacobi Constants
% -------------------------------------------------
%%% Acquire Collinear Lagrange points
L123 = EquilibriumPoints(secondary.MR,1:3); % [3x3] of L1, L2, L3 normalized BCR coordinates

%%% Jacobi constant of Lagrange point
[JC_Lp] = JacobiConstantCalculator(secondary.MR,L123(Lpoint,:),[0,0,0]);

% -------------------------------------------------
% Setting Jacobi value of spacecraft relative to L-point
% -------------------------------------------------
%%% Choose how fast the s/c should be coming through the neck
dJC_vel_mps = 50; % Meters per second (56.5 -> 2e-5 on Europa .... 28.2 -> 5e-6)

%%% Converting this velocity to a jacobi constant value to be
%%% differenced from the value of the chosen L-point
dJC_vel_kps = dJC_vel_mps/1000;
dJC_Lp = (dJC_vel_kps/vNorm)^2;

% % %%% Difference in JC value between Lagrange point and s/c initial state.
% % %%% Basically controls where we want our ZV curves
% % dJC_Lp = 0.000005; % J0 = JC_Lp - dJC_Lp
    
%%% s/c starting JC (JC_scDesired) is lower than JC_Lp because the dJC_Lp is accounted for
%%% in velocity
JC_scInitial = JC_Lp-dJC_Lp;

% -------------------------------------------------
% Finding upper-y-value of neck at L-point
% -------------------------------------------------
%%% Set values
c = JC_scInitial;
u = secondary.MR;
x = L123(Lpoint,1);
z = 0;

%%% JC function equal to zero
f = @(y) x^2 + y.^2 - c + 2*(1-u)./sqrt((x+u)^2+y.^2+z^2) + 2*u./sqrt((x-1+u)^2+y.^2+z^2);

%%% Find the root of this function in the appropriate range
y_neck = fzero(f,[0 6*secondary.R_n]);

%%% Clear variables
clear c u x f z


% -------------------------------------------------
% Finding upper-z-value of neck at L-point
% -------------------------------------------------
%%% Set values
c = JC_scInitial;
u = secondary.MR;
x = L123(Lpoint,1);
y = 0;

%%% JC function equal to zero
f = @(z) x^2 + y.^2 - c + 2*(1-u)./sqrt((x+u)^2+y^2+z.^2) + 2*u./sqrt((x-1+u)^2+y^2+z.^2);

%%% Find the root of this function in the appropriate range
z_neck = fzero(f,[0 6*secondary.R_n]);

%%% Clear variables
clear c u x f y

% -------------------------------------------------
% Setting initial state
% -------------------------------------------------
%%% Initial position
r0_n = [L123(Lpoint,1), -0.8*y_neck, 0.2*secondary.R_n];

%%% With JC0 defined, starting velocity is a function of position. So first
%%% we must calculate the JC of the stationary starting position
JC_initialPos = JacobiConstantCalculator(secondary.MR,r0_n,[0,0,0]);

%%% Starting velocity is found from difference between s/c JC (JC_scDesired) and the
%%% JC of the stationary starting position (JC_initialPos)
dJC_forInitialVelocity = JC_initialPos - JC_scInitial;

%%% Find necessary velocity magnitude
if dJC_forInitialVelocity < 0
    warning('Spacecraft starting in a forbidden zone')
    return
elseif dJC_forInitialVelocity > 0
    v0i = sqrt(abs(dJC_forInitialVelocity));

    % Find direction for velocity
    vHat = R3([-1, 0, 0],0);

    % Create initial velocity vector
    v0_n = vHat .* v0i;
end

% v0_n = [0, 0, 0];
% % v0_n
% % v0_n = v0_n .* 0.5531
% % 989
%%% Store initial conditions
X0_n = [r0_n, v0_n];

% -------------------------------------------------
% Integration options
% -------------------------------------------------
%%% Setting time vector
t_i = 0; % sec
t_f = 4*pi;
dt = t_f/10000;
time0_n = t_i:dt:t_f;

%%% Choosing ode45 tolerance
tol = 1e-10;

%%% Setting integrator options
options = odeset('RelTol',tol,'AbsTol',tol);
options_Apsis = odeset('Events',@event_Apsis,'RelTol',tol,'AbsTol',tol);
options_Impact = odeset('Events',@event_Impact_CR3Bn,'RelTol',tol,'AbsTol',tol);

% -------------------------------------------------
% Integrating reference trajectory
% -------------------------------------------------
%%% Propagating Pre-Burn States
[time_n_ref, X_BCR_n_ref] = ode45(@Int_CR3Bn, time0_n, X0_n, options_Impact, secondary.MR, secondary.R_n);

if plot_referenceTrajectory == 1
    %%% Creating ZV contours
    % Contour bounds
    xCont_min = L123(1,1)-2*secondary.R_n;
    xCont_max = L123(2,1)+2*secondary.R_n;
    yCont_min = -6*secondary.R_n;
    yCont_max = 6*secondary.R_n;
    zCont_min = -3*secondary.R_n;
    zCont_max = 3*secondary.R_n;

    % Creating x-y grid
    xs = linspace(xCont_min,xCont_max,600);
    ys = linspace(yCont_min,yCont_max,200);
    zs = linspace(zCont_min,zCont_max,100);
    [X_xy, Y_xy] = meshgrid(xs,ys);
    [Y_yz, Z_yz] = meshgrid(ys,zs);
    clear xs ys zs

    % Calculating JCs across x-y grid
    JCs_xy = zeros(size(X_xy));
    for xk = 1:size(X_xy,1)
        for yk = 1:size(X_xy,2)
            %%% Zero-Velocity Curve
            zv = JacobiConstantCalculator(secondary.MR,[X_xy(xk,yk), Y_xy(xk,yk), 0] ,[0, 0, 0]);
            JCs_xy(xk,yk) = zv;
        end
    end
    
    % Calculating JCs across y-z grid
    JCs_yz_Lpoint = zeros(size(Y_yz));
    for yk = 1:size(Y_yz,1)
        for zk = 1:size(Y_yz,2)
            %%% Zero-Velocity Curve
            zv = JacobiConstantCalculator(secondary.MR,[L123(Lpoint,1), Y_yz(yk,zk), Z_yz(yk,zk)] ,[0, 0, 0]);
            JCs_yz_Lpoint(yk,zk) = zv;
        end
    end
    
    % ---------------------
    %%% x-y plane
    % ---------------------
    figure; hold all
    plot3(X_BCR_n_ref(:,1),X_BCR_n_ref(:,2),X_BCR_n_ref(:,3),'color',colors.std.blue,'linewidth',1.5)
    plot3(X_BCR_n_ref(1,1),X_BCR_n_ref(1,2),X_BCR_n_ref(1,3),'o','color',colors.std.blue,'linewidth',1.5)
    
    [Cref,href] = contour(X_xy,Y_xy,JCs_xy,[JC_scInitial, JC_scInitial],'color',colors.std.black,'linewidth',3);
    plot3(L123(Lpoint,1),L123(Lpoint,2),L123(Lpoint,3),'m.','markersize',10)
    if v0_n(3) == 0 && r0_n(3) == 0
        plotBody2(secondary.R_n,[1-secondary.MR,0,0],[1 1 1], [0 0 0], 1, .3)
    else
        plotBodyTexture3(secondary.R_n, [1-secondary.MR, 0, 0], secondary.img)
    end
    view(0,90)
    xlim([xCont_min xCont_max])
    ylim([yCont_min yCont_max])
    zlim([zCont_min zCont_max])
    PlotBoi3('X','Y','Z',14)
    axis equal
    
    [Cref,href] = contour(Y_yz,Z_yz,JCs_yz_Lpoint,[JC_scInitial, JC_scInitial],'Visible','off');
    xVals_temp = ones(size(Cref,2)-1).*L123(Lpoint,1);
    plot3(xVals_temp, Cref(1,2:end),Cref(2,2:end),'k','linewidth',3)

    view(79,12)
    return
end

% -------------------------------------------------
% dV Specifications
% -------------------------------------------------
% t_dV = 0.395*pi;
% % t_dV = 0.1*pi;
dVMag_metersPerSecond = 50; % m/s
dVMag_n = dVMag_metersPerSecond/1000/vNorm;
% % dV_headingAngle = (.5)*2*pi; % radians from direction of motion

% -------------------------------------------------
% Loop specifications
% -------------------------------------------------
nAngles = 11; % 21
nDVs = 10; % 150 .... 100 is good at enceladus
% t_dV_f = 0.40; % (pi)
% nAngles = 21;
% nDVs = 150;
nContoursPerBin = 4;

% -------------------------------------------------
% Find indices of reference trajectory to burn at
% -------------------------------------------------
% % % %%% Find 'nDVs' evenly spaced times;
% % % indicesFromImpact = 3;
% % % times_test = linspace(0,time_n_ref(end-indicesFromImpact),nDVs);
% % % 
% % % %%% For each of these times, find closest index of reference
% % % %%% trajectory (since they won't match exactly)
% % % burnIndices = [];
% % % for t = times_test
% % %     ti = find(abs(time_n_ref - t) == min(abs(time_n_ref - t)));
% % %     burnIndices = [burnIndices, ti];
% % % end

%%% Calculate distance traveled at each state
distTraveled = [0];
sum = 0;
for kk = 2:size(X_BCR_n_ref,1)
    distTraveled(kk) = distTraveled(kk-1) + norm(X_BCR_n_ref(kk,1:3)-X_BCR_n_ref(kk-1,1:3));
end

%%% Set target distances to burn at
burnDistances = linspace(0,distTraveled(end-1),nDVs);

%%% For each of these distances, find closest index of reference
%%% trajectory (since they won't match exactly)
burnIndices = [];
for kk = 1:length(burnDistances)
    ti = find(abs(distTraveled - burnDistances(kk)) == min(abs(distTraveled - burnDistances(kk))));
    burnIndices = [burnIndices, ti];
end

% -------------------------------------------------
% Set various heading angles for burns
% -------------------------------------------------
headingAngles = 0:(2*pi/nAngles):(2*pi)-(2*pi/nAngles);

% ========================================================================
%%% Looping through integrations and storing periapses
% ========================================================================
trajCount = 0;
for bi = burnIndices
for dV_headingAngle = headingAngles

trajCount = trajCount + 1;
% -------------------------------------------------
% Preallocation
% -------------------------------------------------
X_apses_SCR_n = [];
time_apses = [];

% -------------------------------------------------
% Integration
% -------------------------------------------------
%%% Assigning time of burn
t_dV = time_n_ref(bi);

%%% Assigning pre-dV state
X_predV_n = X_BCR_n_ref(bi,:);

%%% Adding tangent dV to states
% Finding vHat
vHat = X_predV_n(4:6)./norm(X_predV_n(4:6));
vHat = R3(vHat,dV_headingAngle);
% Determining dV
dV_n = vHat .* dVMag_n;
% Creating post-dV state
X0_postdV_n = [X_predV_n(1:3), X_predV_n(4:6)+dV_n];

%%% Propagating the post-burn States with Apsis event
if run_apsesAnalysis == 1
    [time_n_Apsis, X_BCR_n_Apsis, time_eventApsis, X_eventApsis, index_eventApsis] = ode45(@Int_CR3Bn, [t_dV:dt:t_f], X0_postdV_n,...
                                                           options_Apsis, secondary.MR, secondary.R_n);
    %%% Storing periapses and apoapses
    X_apses_SCR_n = [X_apses_SCR_n; X_eventApsis];
    time_apses = [time_apses; time_eventApsis];
end

 %%% Propagating the post-burn States with Impact event
[time_n_postBurn, X_BCR_n_postBurn, time_eventImpact, X_eventImpact, index_eventImpact] = ode45(@Int_CR3Bn, [t_dV:dt:t_f], X0_postdV_n,...
                                                       options_Impact, secondary.MR, secondary.R_n);
% % figure(1)
% % subplot(2,4,[1 2]); hold all
% % plot3(X_BCR_n_postBurn(:,1),X_BCR_n_postBurn(:,2),X_BCR_n_postBurn(:,3),'color',colors.std.black)

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
    %%% Norm of impact velocity
    impactVelocities(kk) = norm(X_eventImpact(kk,4:6));
    
    %%% Creating SCR position vector
    rImpact_SCR_n = X_eventImpact(kk,1:3) - [1-secondary.MR,0,0];
    rImpact_SCR_n = rImpact_SCR_n./norm(rImpact_SCR_n);
    
    %%% Velocity unit vector at impact
    vHatImpact_n = X_eventImpact(kk,4:6)./norm(X_eventImpact(kk,4:6));
    
    %%% Angle between velocity and surface
    A = R3(rImpact_SCR_n,pi/2);
    B = vHatImpact_n;
    impactAngles(kk) = acos(dot(A,B)/(norm(A)*norm(B)));
    if impactAngles(kk) > pi/2
        impactAngles(kk) = pi - impactAngles(kk);
    end
end

% -------------------------------------------------
% Post-burn Jacobi constant
% -------------------------------------------------
JC_postBurn = JacobiConstantCalculator(secondary.MR,[X_BCR_n_postBurn(1,1),X_BCR_n_postBurn(1,2),X_BCR_n_postBurn(1,3)],[X_BCR_n_postBurn(1,4),X_BCR_n_postBurn(1,5),X_BCR_n_postBurn(1,6)]);

% -------------------------------------------------
% Storing data 
% -------------------------------------------------
if isempty(X_eventImpact) == 0 % If there was an impact within the time
    impactTraj(trajCount).X_eventImpact    = X_eventImpact;
    impactTraj(trajCount).impactLatLons    = impactLatLons;
    impactTraj(trajCount).impactSpeeds     = impactVelocities;
    impactTraj(trajCount).time_eventImpact = time_eventImpact;
    impactTraj(trajCount).impactAngles     = impactAngles;
    impactTraj(trajCount).JC_postBurn      = JC_postBurn;
    impactTraj(trajCount).time_burn        = t_dV; 
    impactTraj(trajCount).burnQuivers      = [X0_postdV_n(1:3), dV_n];
    impactTraj(trajCount).burnHat          = vHat;
    if plot_postBurnTrajectories == 1
        impactTraj(trajCount).X_BCR_n_postBurn = X_BCR_n_postBurn;
    end
elseif isempty(X_eventImpact) == 1 % No impact within the time
% %     impactTraj(trajCount).X_eventImpact    = NaN(1,6);
% %     impactTraj(trajCount).impactLatLons    = NaN(1,2);
% %     impactTraj(trajCount).impactSpeeds     = NaN;
% %     impactTraj(trajCount).time_eventImpact = NaN;
% %     impactTraj(trajCount).impactAngles     = NaN;
% %     impactTraj(trajCount).JC_postBurn      = JC_postBurn;
% %     impactTraj(trajCount).time_burn        = t_dV; 
% %     impactTraj(trajCount).burnQuivers      = [X0_postdV_n(1:3), dV_n];
% %     impactTraj(trajCount).burnHat          = vHat;
% %     if plot_postBurnTrajectories == 1
% %         impactTraj(trajCount).X_BCR_n_postBurn = X_BCR_n_postBurn;
% %     end
end


end
end

% ========================================================================
%%% Post-Loop Work
% ========================================================================
% -------------------------------------------------
% Creating JC bins
% -------------------------------------------------
%%% Determine minimum and maximum JC after burn
JC_min = min([impactTraj(:).JC_postBurn]);
JC_max = max([impactTraj(:).JC_postBurn]);

%%% Bin count for JC values
binCount_JC = 4; 

%%% Create colors for bins
colorMatrix = [colors.sch.d4_1(2,:);colors.sch.d4_1(1,:);colors.sch.d4_1(4,:);colors.sch.d4_1(3,:)];

%%% Create values for bins
bins_JC = linspace(JC_min, JC_max, binCount_JC+1);
bins_JC(end) = bins_JC(end)+1; % make last bin big so everything fits in first set

%%% Assign bins and colors
for kk = 1:size(impactTraj,2)
    impactTraj(kk).bin_JC = discretize(impactTraj(kk).JC_postBurn,bins_JC);
end

% ========================================================================
%%% Sorting plot details by bin for faster plotting
% ========================================================================
%%% Predefining structure 
binData(1).X_eventImpact    = [];
binData(1).impactLatLons    = [];
binData(1).impactSpeeds     = [];
binData(1).impactSpeeds_mps = [];
binData(1).time_eventImpact = [];
binData(1).impactAngles     = [];
binData(1).JC_postBurn      = [];
binData(1).time_burn        = [];
binData(1).burnQuivers      = [];
binData(1).burnHat          = [];
if plot_postBurnTrajectories == 1
    binData(1).X_BCR_n_postBurn   = [];
end

%%% Assign bin colors
for kk = 1:binCount_JC
    binData(kk).bin_col = colorMatrix(kk,:);
end

cc = 0;
for kk = 1:size(impactTraj,2)
    %%% Find bin index (1:binCount_JC)
    bi = impactTraj(kk).bin_JC;
    
    %%% Store data into correct bin
    if isempty(impactTraj(kk).X_eventImpact) == 0
        binData(bi).X_eventImpact    = [binData(bi).X_eventImpact;    impactTraj(kk).X_eventImpact];
        binData(bi).impactLatLons    = [binData(bi).impactLatLons;    impactTraj(kk).impactLatLons];
        binData(bi).impactSpeeds     = [binData(bi).impactSpeeds;     impactTraj(kk).impactSpeeds];
        binData(bi).impactSpeeds_mps = [binData(bi).impactSpeeds_mps; impactTraj(kk).impactSpeeds*vNorm*1000];
        binData(bi).time_eventImpact = [binData(bi).time_eventImpact; impactTraj(kk).time_eventImpact];
        binData(bi).impactAngles     = [binData(bi).impactAngles;     impactTraj(kk).impactAngles];
        binData(bi).JC_postBurn      = [binData(bi).JC_postBurn;      impactTraj(kk).JC_postBurn];
        binData(bi).time_burn        = [binData(bi).time_burn;        impactTraj(kk).time_burn];
        binData(bi).burnQuivers      = [binData(bi).burnQuivers;      impactTraj(kk).burnQuivers];
        binData(bi).burnHat          = [binData(bi).burnHat;          impactTraj(kk).burnHat];
        if plot_postBurnTrajectories == 1
            binData(bi).X_BCR_n_postBurn   = [binData(bi).X_BCR_n_postBurn; NaN(1,6); impactTraj(kk).X_BCR_n_postBurn];
        end
        
    elseif isempty(impactTraj(kk).X_eventImpact) == 1
%         binData(bi).X_eventImpact    = [binData(bi).X_eventImpact;    impactTraj(kk).X_eventImpact];
%         binData(bi).impactLatLons    = [binData(bi).impactLatLons;    impactTraj(kk).impactLatLons];
%         binData(bi).impactSpeeds     = [binData(bi).impactSpeeds;     impactTraj(kk).impactSpeeds];
%         binData(bi).time_eventImpact = [binData(bi).time_eventImpact; impactTraj(kk).time_eventImpact];
%         binData(bi).impactAngles     = [binData(bi).impactAngles;     impactTraj(kk).impactAngles];
% %         binData(bi).JC_postBurn      = [binData(bi).JC_postBurn;      impactTraj(kk).JC_postBurn];
% %         binData(bi).time_burn        = [binData(bi).time_burn;        impactTraj(kk).time_burn];
% %         binData(bi).burnQuivers      = [binData(bi).burnQuivers;      impactTraj(kk).burnQuivers];
% %         binData(bi).burnHat          = [binData(bi).burnHat;          impactTraj(kk).burnHat];
%         if plot_postBurnTrajectories == 1
%             binData(bi).X_BCR_n_postBurn   = [binData(bi).X_BCR_n_postBurn;   impactTraj(kk).X_BCR_n_postBurn];
%         end
    end
    
end

% ========================================================================
%%% Plots
% ========================================================================
% -------------------------------------------------
% Set up plots
% -------------------------------------------------
markerSize = 13;
% ------------------
% figure(1), subplot(2,4,[1 2]) - system, apses, and impacts
% ------------------
figure(1); 
% set(gcf,'Position',[266 54 986 678]); % big
set(gcf,'Position',[1 1 1440 804]); % fullscreen mac
subplot(2,4,[1 2]); hold all
axis equal
PlotBoi3('x$_n$','y$_n$','z$_n$',16,'Latex') 

%%% Title
title_body = secondary.name;
title_J0 = ['J$_0$ = J',sprintf('$_{L%1.0d}$ + %2.2e',Lpoint,-dJC_Lp)];
title_v0 = ['$\mathopen|v_0 \mathclose|$ = ',sprintf('%2.1f',sqrt(dJC_Lp)*vNorm*1000), ' $m/s$'];
title_R_n = sprintf('$R_n$ = %1.2e',secondary.R_n);
title_dV = ['$\mathopen|\Delta V \mathclose|$ = ',sprintf('%2.0f',dVMag_metersPerSecond), ' $m/s$'];

title({[title_body, ', \hspace{0.5cm}', title_J0, ', \hspace{0.5cm}', title_v0,...
         ', \hspace{0.5cm}', title_R_n, ', \hspace{0.5cm}', title_dV], ''},'Interpreter','Latex');

%%% Creating ZV contours
% Contour bounds
xCont_min = L123(1,1)-2*secondary.R_n;
xCont_max = L123(2,1)+2*secondary.R_n;
yCont_min = -6*secondary.R_n;
yCont_max = 6*secondary.R_n;

% Creating x-y grid
xs = linspace(xCont_min,xCont_max,600);
ys = linspace(yCont_min,yCont_max,200);
[X_xy, Y_xy] = meshgrid(xs,ys);
clear xs ys

% Calculating JCs across grid
JCs_xy = zeros(size(X_xy));
for xk = 1:size(X_xy,1)
    for yk = 1:size(X_xy,2)
        %%% Zero-Velocity Curve
        zv = JacobiConstantCalculator(secondary.MR,[X_xy(xk,yk), Y_xy(xk,yk), 0] ,[0, 0, 0]);
        JCs_xy(xk,yk) = zv;
    end
end

%%% Plot limit
xlim([xCont_min xCont_max])
ylim([yCont_min yCont_max])

%%% Plotting secondary
plotBody2(secondary.R_n,[1-secondary.MR,0,0],[1,1,1],colors.std.black,1.5,0)

% ------------------
% Colorbar for whole figure
% ------------------
h8 = get(subplot(2,4,8),'Position');
cbar = colorbar('Position', [h8(1)+h8(3)+0.01  h8(2)  0.025  h8(2)+h8(3)*4.6]);
colormap(colorMatrix)
caxis([0 size(colorMatrix,1)]);
cbar.FontName     = 'Arial';
cbar.FontSize     = 10;
cbar.Ticks        = sort([0:binCount_JC]);
TickLabels = num2cell(bins_JC(1:end-1)-JC_Lp);
for kk = 1:length(TickLabels)
    TickLabels{kk} = sprintf('%1.2e',TickLabels{kk});
end
cbar.TickLabels   = [TickLabels, sprintf('J_{L%1.0d} +',Lpoint)];
cbar.Label.String = {'Post-Burn Jacobi Constant'};
cbar.Label.Rotation = 0; % 0 = horizontal, 90 = vertical
cbar.Label.Position = [h8(4)+.2, 4.15, 0];

% -------------------------------------------------
% Loop through bins and plot data
% -------------------------------------------------
% for kk = 4
for kk = 1:binCount_JC
    %%% Check if this bin has any impact
    if isempty(binData(kk).X_eventImpact) == 1
        continue
    end
    % -------------------------------------------------
    % figure(1), subplot(2,4,[1 2]) - system, apses, and impacts
    % -------------------------------------------------
    figure(1); subplot(2,4,[1 2]); hold all
%     axis equal
%     PlotBoi3('x$_n$','y$_n$','z$_n$',16,'Latex')    
    
    %%% Plotting post-burn contours
    orderedPostBurnJCs = sort(binData(kk).JC_postBurn);
    contoursForPlotting = linspace(orderedPostBurnJCs(1),orderedPostBurnJCs(end),nContoursPerBin);
    for currentContour = contoursForPlotting
%     for jj = 1:length(binData(kk).JC_postBurn)
%         if binData(kk).JC_postBurn(jj) == min(binData(kk).JC_postBurn) || binData(kk).JC_postBurn(jj) == max(binData(kk).JC_postBurn)
%         if rem(jj,10) == 0
%             [C2,h2] = contour(X,Y,Z,[binData(kk).JC_postBurn(jj), binData(kk).JC_postBurn(jj)],'color',binData(kk).bin_col,'linewidth',2);
%         end
        [C2,h2] = contour(X_xy,Y_xy,JCs_xy,[currentContour, currentContour],'color',binData(kk).bin_col,'linewidth',3);
    end
    
    % %%% Plotting post-burn FULL trajectory
    if plot_postBurnTrajectories == 1
        plot3(binData(kk).X_BCR_n_postBurn(:,1),binData(kk).X_BCR_n_postBurn(:,2),binData(kk).X_BCR_n_postBurn(:,3),'linewidth',1,'color',binData(kk).bin_col)
    end
    
    %%% Plotting apses
    if run_apsesAnalysis == 1
        if isempty(X_apses_SCR_n) == 0
            plot3(X_apses_SCR_n(:,1),X_apses_SCR_n(:,2),X_apses_SCR_n(:,3),'o','markersize',5,...
                        'markeredgecolor',colors.std.black,'markerfacecolor',colors.std.ltgrn,'linewidth',2)
        end
    end

    %%% Plotting impact locations
    if isempty(binData(kk).X_eventImpact) == 0
%         plot3(binData(kk).X_eventImpact(:,1),binData(kk).X_eventImpact(:,2),binData(kk).X_eventImpact(:,3),'o','markersize',8,...
%                         'markeredgecolor',colors.std.black,'markerfacecolor',binData(kk).bin_col,'linewidth',2)
        plot3(binData(kk).X_eventImpact(:,1),binData(kk).X_eventImpact(:,2),binData(kk).X_eventImpact(:,3),'.','markersize',markerSize,...
                        'color',binData(kk).bin_col,'linewidth',2)
    end
    
    %%% Plotting burn quivers
    if plot_burnQuivers == 1
        qS = 0.002; % quiverScale
        quiver3(binData(kk).burnQuivers(:,1),binData(kk).burnQuivers(:,2),binData(kk).burnQuivers(:,3),binData(kk).burnHat(:,1).*qS,...
                binData(kk).burnHat(:,2).*qS,binData(kk).burnHat(:,3).*qS,0,'color',binData(kk).bin_col,'linewidth',0.75,'maxheadsize',0.4);
    end
    
    % -------------------------------------------------
    % figure(1), subplot(2,4,3) - impact longitude vs post burn Jacobi
    % value
    % -------------------------------------------------
    subplot(2,4,3); hold all
    xlim([-180 180])
    ylim([min(vertcat(binData(:).JC_postBurn))*.9999 max(vertcat(binData(:).JC_postBurn))*1.0001])
    PlotBoi2('Longitude, $^\circ$','Post-Burn Jacobi Value',12,'Latex')
    
    %%% If there was an impact
    if isempty(binData(kk).X_eventImpact) == 0
        %%% Plotting impact longitude vs post burn jacobi value
        plot(binData(kk).impactLatLons(:,2),binData(kk).JC_postBurn,'.','markersize',markerSize,...
                        'color',binData(kk).bin_col,'linewidth',2)

        %%% Plotting anti-primary line (0-deg longitude)
        [antiPrimaryLine] = plot([0 0],[0 max(vertcat(binData(:).JC_postBurn))].*(1.1),'color',colors.std.blue,'linewidth',1.5);

        %%% legend
        legend([antiPrimaryLine],{'Anti-Primary'},'Interpreter','latex','location','best')
    end
    
    % -------------------------------------------------
    % figure(1), subplot(2,4,4) - impact longitude vs impact angle
    % -------------------------------------------------
    subplot(2,4,4); hold all
    xlim([-180 180])
    ylim([0 90])
    PlotBoi2('Longitude, $^\circ$','Impact Angle, $^\circ$',12,'Latex')
    
    %%% If there was an impact
    if isempty(binData(kk).X_eventImpact) == 0
        %%% Plotting impact longitude vs post burn jacobi value
        plot(binData(kk).impactLatLons(:,2),binData(kk).impactAngles.*(180/pi),'.','markersize',markerSize,...
                        'color',binData(kk).bin_col,'linewidth',2)
        %%% Plotting anti-primary line (0-deg longitude)
        [antiPrimaryLine] = plot([0 0],[0 max(vertcat(binData(:).impactAngles))].*(1.1*(180/pi)),'color',colors.std.blue,'linewidth',1.5);

        %%% legend
        legend([antiPrimaryLine],{'Anti-Primary'},'Interpreter','latex','location','best')
    
    end
    
    % -------------------------------------------------
    % figure(1), subplot(2,4,5) - impact longitude vs impact speed
    % -------------------------------------------------
    subplot(2,4,5); hold all
    xlim([-180 180])
    ylim([min(vertcat(binData(:).impactSpeeds_mps))*.999 max(vertcat(binData(:).impactSpeeds_mps))*1.001])
    PlotBoi2('Longitude, $^\circ$','Impact Speed, $m/s$',12,'Latex')

    %%% If there was an impact
    if isempty(binData(kk).X_eventImpact) == 0
        %%% Plotting impact longitude vs speed
        plot(binData(kk).impactLatLons(:,2),binData(kk).impactSpeeds_mps,'.','markersize',markerSize,...
                        'color',binData(kk).bin_col,'linewidth',2)

        %%% Plotting anti-primary line (0-deg longitude)
        [antiPrimaryLine] = plot([0 0],[0 max(vertcat(binData(:).impactSpeeds_mps))].*1.1,'color',colors.std.blue,'linewidth',1.5);

        %%% legend
        legend([antiPrimaryLine],{'Anti-Primary'},'Interpreter','latex','location','best')
    end

    % -------------------------------------------------
    % figure(1), subplot(2,4,6) - impact longitude vs impact time
    % -------------------------------------------------
    subplot(2,4,6); hold all
    xlim([-180 180])
    PlotBoi2('Longitude, $^\circ$','Time to Impact, $\pi$',12,'Latex')

    %%% Plotting impact angle vs impact time
    if isempty(binData(kk).X_eventImpact) == 0
        %%% Plotting impact longitude vs speed
        plot(binData(kk).impactLatLons(:,2),binData(kk).time_eventImpact./pi,'.','markersize',markerSize,...
                        'color',binData(kk).bin_col,'linewidth',2)
        
        %%% Plotting time of reference impact
        [referenceImpactLine] = plot([-180 180],[time_n_ref(end) time_n_ref(end)]./pi,'color',colors.std.black,'linewidth',1.5);

        %%% Plotting anti-primary line (0-deg longitude)
        [antiPrimaryLine] = plot([0 0],[0 max(vertcat(binData(:).time_eventImpact))].*(1.1/pi),'color',colors.std.blue,'linewidth',1.5);

        %%% legend
        legend([antiPrimaryLine, referenceImpactLine],{'Anti-Primary','Reference Impact'},'Interpreter','latex','location','best')
    end
    
    % -------------------------------------------------
    % figure(1), subplot(2,4,7) - impact longitude vs time of dV
    % -------------------------------------------------
    subplot(2,4,7); hold all
    xlim([-180 180])
    ylim([0 max(vertcat(binData(:).time_burn))*1.1/pi])
    PlotBoi2('Longitude, $^\circ$','Time of $\Delta V$, $\pi$',12,'Latex')

    %%% Plotting impact angle vs impact time
    if isempty(binData(kk).X_eventImpact) == 0
        %%% Plotting impact longitude vs speed
        plot(binData(kk).impactLatLons(:,2),binData(kk).time_burn./pi,'.','markersize',markerSize,...
                        'color',binData(kk).bin_col,'linewidth',2)

        %%% Plotting anti-primary line (0-deg longitude)
        [antiPrimaryLine] = plot([0 0],[0 max(vertcat(binData(:).time_eventImpact))].*(1.1),'color',colors.std.blue,'linewidth',1.5);

        %%% legend
        legend([antiPrimaryLine],{'Anti-Primary'},'Interpreter','latex','location','best')
    end
    
    % -------------------------------------------------
    % figure(1), subplot(2,4,8) - impact angle vs impact speed
    % -------------------------------------------------
    subplot(2,4,8); hold all
    xlim([0 90])
    ylim([min(vertcat(binData(:).impactSpeeds_mps))*.999 max(vertcat(binData(:).impactSpeeds_mps))*1.001])
    PlotBoi2('Impact Angle, $^\circ$', 'Impact Speed, $m/s$',12,'Latex')

    %%% Plotting impact angle vs speed
    if isempty(binData(kk).X_eventImpact) == 0
        plot(binData(kk).impactAngles.*(180/pi),binData(kk).impactSpeeds_mps,'.','markersize',markerSize,...
                        'color',binData(kk).bin_col,'linewidth',2)
    end
    
end

% -------------------------------------------------
% Plot things to be on top
% -------------------------------------------------
% ------------------
% figure(1), subplot(2,4,[1 2 3 4]) - system, apses, and impacts
% ------------------
figure(1); subplot(2,4,[1 2]); hold all
%%% Plotting pre-burn trajectory
plot3(X_BCR_n_ref(:,1),X_BCR_n_ref(:,2),X_BCR_n_ref(:,3),'linewidth',2,'color',colors.std.black)

%%% Plotting Lagrange Points
plot3(L123(1:2,1), L123(1:2,2), L123(1:2,3), '^', 'markeredgecolor', colors.std.black,...
                                       'markerfacecolor',colors.std.mag,'markersize',10);

%%% Plotting pre-burn contour
[C1,h1] = contour(X_xy,Y_xy,JCs_xy,[JC_scInitial, JC_scInitial],'color',colors.std.black,'linewidth',2);

toc






