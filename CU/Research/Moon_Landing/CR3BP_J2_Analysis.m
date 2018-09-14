clear
clc
close all
addpath(genpath('/Users/lukebury/Documents/MATLAB/mbin'))

% ========================================================================
%%% Run Switches
% ========================================================================
run_J2Traj         = 1;
run_JCverification = 0;
run_zeroVelocity   = 0;
run_J2Potential    = 0;
run_EquiMovement   = 0;
% ========================================================================
%%% Importing Data
% ========================================================================
%%% General data on solar system bodies
bodies = getBodyData();

%%% Color options/schemes
colors = get_colors();

% ========================================================================
%%% Constant / Conditions
% ========================================================================
%%% Selecting bodies
primary = bodies.jupiter;
secondary = bodies.europa;

%%% Normalizing constants
rNorm = secondary.a;         % n <-> km
tNorm = 1/secondary.meanMot; % n <-> sec
vNorm = rNorm / tNorm;       % n <-> km/sec

%%% Acquire Lagrange points of non-J2 case
L_points_n = EquilibriumPoints(secondary.MR);
L1_n = L_points_n(1,:); % [1x3] normalized barycentric coordinates
L2_n = L_points_n(2,:); % [1x3] normalized barycentric coordinates
L3_n = L_points_n(3,:); % [1x3] normalized barycentric coordinates

% ========================================================================
%%% Integrating with and without J2
% ========================================================================
if run_J2Traj == 1
%%% Initial State
% positionNudge = [(1-secondary.MR)*rNorm+3*secondary.R, 0, 0]./rNorm; % enter in km
% velocityNudge = [0, 0, 800]./1000./vNorm; % enter in m/s
% velocityNudge = [0, 550, 550]./1000./vNorm; % enter in m/s
% 
% positionNudge = [0, 0, 0]./rNorm; % enter in km
% velocityNudge = [-1, -20, 130]./1000./vNorm; % enter in m/s

positionNudge = [.5-secondary.MR, sqrt(3)/2, 0]
velocityNudge = [0,0,0]

%%% Initial State
% x0_n = L2_n + positionNudge;
x0_n = positionNudge;
xdot0_n = [0, 0, 0] + velocityNudge;
X0_n = [x0_n xdot0_n]';

%%% Setting time vector and normalizing 
ti = 0; % sec
dt = 60; % sec
tf = 3600*24*3.65; % sec
time0_n = [ti:dt:tf] ./ tNorm;

%%% Choosing ode45 tolerance
tol = 1e-9;

%%% Setting integrator options
options = odeset('Events',@impactEvent_CR3Bn,'RelTol',tol,'AbsTol',tol);
options_J2 = odeset('Events',@impactEvent_CR3Bn_J2,'RelTol',tol,'AbsTol',tol);

%%% Propagating the States without J2
[time_n, X_BCR_n] = ode45(@Int_CR3Bn, time0_n, X0_n, options, secondary.MR, secondary.R_n);
    
%%% Propagating the States with J2
[time_J2_n, X_BCR_J2_n] = ode45(@Int_CR3Bn_J2, time0_n, X0_n, options_J2, secondary.MR, primary.R/rNorm, secondary.R_n, primary.J2, 0);

%%% Acquiring modified jacobi constant values
[JCs] = JacobiConstantCalculator_J2(secondary.MR,X_BCR_J2_n(:,1:3),X_BCR_J2_n(:,4:6), primary.R/rNorm, secondary.R_n, primary.J2, 0);

%%% Determining shorter simulation, for plotting/comparison purposes
if length(time_n) < length(time_J2_n)
    timeShort = time_n;
else
    timeShort = time_J2_n;
end
% -----------------------------
% Rotating Plot
% -----------------------------
figure; hold all
plot3(X_BCR_n(:,1),X_BCR_n(:,2),X_BCR_n(:,3),'linewidth',2)
plot3(X_BCR_J2_n(:,1),X_BCR_J2_n(:,2),X_BCR_J2_n(:,3),'linewidth',2)
plotBodyTexture3(secondary.R_n, [1-secondary.MR, 0, 0],secondary.img);
% plotBody3(secondary.R_n, [1-secondary.MR, 0, 0], secondary.color)
PlotBoi3('X_n','Y_n','Z_n',14)
legend('No J2','W/ J2')
axis equal
view(0,90)
grid on

% -----------------------------
% Plotting Differences in Position
% -----------------------------
figure;
subplot(3,1,1);
plot(timeShort.*tNorm./3600,(X_BCR_n(1:length(timeShort),1)-X_BCR_J2_n(1:length(timeShort),1)).*rNorm,'linewidth',2);
PlotBoi2('','x-Diff, km',12)
subplot(3,1,2);
plot(timeShort.*tNorm./3600,(X_BCR_n(1:length(timeShort),2)-X_BCR_J2_n(1:length(timeShort),2)).*rNorm,'linewidth',2);
PlotBoi2('','y-Diff, km',12)
subplot(3,1,3);
plot(timeShort.*tNorm./3600,(X_BCR_n(1:length(timeShort),3)-X_BCR_J2_n(1:length(timeShort),3)).*rNorm,'linewidth',2);
PlotBoi2('Time, hr','z-Diff, km',12)

if run_JCverification == 1
    figure
    plot(percentchange(JCs))
    PlotBoi2('','% Change',14)
end


end % end run_J2Traj

% ========================================================================
%%% Zero Velocity Curves
% ========================================================================
if run_zeroVelocity == 1
%%% Test system
xs = linspace(-.1,1.1,1000);
% ys = linspace(-2*secondary.R_n,2*secondary.R_n,100);
ys = linspace(-.2,.2,1000);
[X, Y] = meshgrid(xs,ys);
clear xs ys

    Z_J2 = zeros(size(X));
    Z = zeros(size(X));
    for xk = 1:size(X,1)
        for yk = 1:size(X,2)
            
            %%% J2 Zero-Velocity Curve
            zv_J2 = JacobiConstantCalculator_J2(secondary.MR,[X(xk,yk), Y(xk,yk), 0] ,[0, 0, 0] , primary.R/rNorm, secondary.R_n, primary.J2, 0);
            Z_J2(xk,yk) = zv_J2;
            
            %%% Zero-Velocity Curve
            zv = JacobiConstantCalculator(secondary.MR,[X(xk,yk), Y(xk,yk), 0] ,[0, 0, 0]);
            Z(xk,yk) = zv;
        end
    end

    
    figure;
    subplot(2,1,1); hold all; title('J2')
    [C_J2,h_J2] = contour(X,Y,Z_J2,[3.001, 3.00361, 3.0038, 3.005,3.006,3.007, 3.009, 3.4, 9, 15],'color','k');
    plotBody2( secondary.R_n, [1-secondary.MR,0,0], secondary.color )
    plotBody2( primary.R/rNorm, [-secondary.MR,0,0], primary.color )
    axis equal; clabel(C_J2,h_J2)
    xlim([0.93, 1.07]); ylim([-0.02, 0.02])
    PlotBoi2('','y_n - BCR',14)
    
    subplot(2,1,2); hold all; title('No J2')
    [C,h] = contour(X,Y,Z,[3.001, 3.00361, 3.0038, 3.005,3.006,3.007, 3.009, 3.4, 9, 15],'color','k');
    plotBody2( secondary.R_n, [1-secondary.MR,0,0], secondary.color )
    plotBody2( primary.R/rNorm, [-secondary.MR,0,0], primary.color )
    axis equal; clabel(C,h)
    xlim([0.93, 1.07]); ylim([-0.02, 0.02])
    PlotBoi2('x_n - BCR','y_n - BCR',14)
end

% ========================================================================
%%% J2 Potential / Equillibrium points
% ========================================================================
if run_J2Potential == 1
     
%%% Acquire Lagrange points of both models
L_points123_n = EquilibriumPoints(secondary.MR,1:3);
L_points123_J2n = EquilibriumPoints_J2(secondary.MR,primary.J2,0,primary.R/rNorm,secondary.R_n);
LsDiff = L_points123_J2n-L_points123_n;
LsDiff.*rNorm


%%% Initial State
X0_n   = [L_points123_n(1,1:3),0,0,0]';
X0_J2n = [L_points123_J2n(1,1:3),0,0,0]';


%%% Setting time vector and normalizing 
ti = 0; % sec
dt = 10; % sec
tf = 3600*24; % sec
time0_n = [ti:dt:tf] ./ tNorm;

%%% Choosing ode45 tolerance
tol = 1e-10;

%%% Setting integrator options
options = odeset('Events',@impactEvent_CR3Bn,'RelTol',tol,'AbsTol',tol);
options_J2 = odeset('Events',@impactEvent_CR3Bn_J2,'RelTol',tol,'AbsTol',tol);

%%% Propagating the States without J2
[time_n, X_BCR_n_eqTest] = ode45(@Int_CR3Bn, time0_n, X0_n, options, secondary.MR, secondary.R_n);
    
%%% Propagating the States with J2
[time_J2_n, X_BCR_J2_n_eqTest] = ode45(@Int_CR3Bn_J2, time0_n, X0_J2n, options_J2, secondary.MR, primary.R/rNorm, secondary.R_n, primary.J2, 0);


figure
subplot(3,2,1)
plot(time_n,(X_BCR_n_eqTest(:,1)-X_BCR_n_eqTest(1,1)).*rNorm)
PlotBoi2('','No J2, x-dev',10)
subplot(3,2,3)
plot(time_n,(X_BCR_n_eqTest(:,2)-X_BCR_n_eqTest(1,2)).*rNorm)
PlotBoi2('','No J2, y-dev',10)
subplot(3,2,5)
plot(time_n,(X_BCR_n_eqTest(:,3)-X_BCR_n_eqTest(1,3)).*rNorm)
PlotBoi2('','No J2, z-dev',10)
subplot(3,2,2)
plot(time_J2_n,(X_BCR_J2_n_eqTest(:,1)-X_BCR_J2_n_eqTest(1,1)).*rNorm)
PlotBoi2('','w/ J2, x-dev',10)
subplot(3,2,4)
plot(time_J2_n,(X_BCR_J2_n_eqTest(:,2)-X_BCR_J2_n_eqTest(1,2)).*rNorm)
PlotBoi2('','w/ J2, y-dev',10)
subplot(3,2,6)
plot(time_J2_n,(X_BCR_J2_n_eqTest(:,3)-X_BCR_J2_n_eqTest(1,3)).*rNorm)
PlotBoi2('','w/ J2, z-dev',10)
% % Setting up system
% u = secondary.MR;
% R1 = primary.R/rNorm;
% J21 = primary.J2;
% J22 = 0; R2 = 0;
% 
% x = linspace(-2,2,1000000);
% U1 = .5*x.^2 + (1-u)./abs(x+u) + u./abs(x-1+u);
% U2 = .5*x.^2 + (1-u)./abs(x+u) + u./abs(x-1+u) + (1-u)*J21*R1*R1./(2*abs(x+u).^3) + u*J22*R2*R2./(2*abs(x-1+u).^3);
% 
% figure;
% subplot(2,1,1)
% plot(x,U1); ylim([0 5])
% PlotBoi2('','U(x), no J2',12)
% subplot(2,1,2)
% plot(x,U2); ylim([0 5])
% PlotBoi2('','U(x), w/ J2',12)
% 
% % x positions of L3, L1, L2 in BCR
% U1_E = [x(find(U1==min(U1(1:500000)))), x(find(U1==min(U1(500000:750000)))),x(find(U1==min(U1(750000:1000000))))]
% U2_E = [x(find(U2==min(U2(1:500000)))), x(find(U2==min(U2(500000:750000)))),x(find(U2==min(U2(750000:1000000))))]
% warning('Not a dynamic search')
% % Difference between L3, L1, and L2 in km
% U_diff = (U1_E - U2_E).*rNorm

end

% ========================================================================
%%% Analysis of Movement of Linear Equilibrium Points for Different Systems
% ========================================================================
if run_EquiMovement == 1
%%% Setting 3B systems to loop through
systems = {bodies.earth,bodies.moon;...
          bodies.mars,bodies.phobos;...
          bodies.mars,bodies.deimos;...
          bodies.jupiter,bodies.io;... 
          bodies.jupiter,bodies.europa;...
          bodies.saturn,bodies.enceladus;...
          bodies.saturn,bodies.titan;...
          bodies.neptune,bodies.triton};

%%% Preallocating results
dataTable = {'Primary','Secondary','J2,2 Considered','L1 Diff_n','L2 Diff_n','L3 Diff_n','L1 Diff','L2 Diff','L3 Diff'};

%%% Looping through systems and building table of L123 differences
for kk = 1:size(systems,1)
    %%% Reassigning data for clarity
    primary   = systems{kk,1};
    secondary = systems{kk,2};
    
    %%% Normalizing constants
    rNorm = secondary.a; % n <-> km
    
    %%% Acquire Lagrange points of both models
    L_points123_n = EquilibriumPoints(secondary.MR,1:3);
%     if isfield(secondary,'J2') == 0
%         secondary.J2 = 0;
%     end
%     L_points123_J2n = EquilibriumPoints_J2(secondary.MR,primary.J2,secondary.J2,primary.R/rNorm,secondary.R_n);
    % Using J2,2 = 0
    L_points123_J2n = EquilibriumPoints_J2(secondary.MR,primary.J2,0,primary.R/rNorm,secondary.R_n);
    LsDiff = L_points123_J2n-L_points123_n;

    %%% Storing data
    dataTable{kk+1,1} = primary.name;
    dataTable{kk+1,2} = secondary.name;
    dataTable{kk+1,4} = LsDiff(1);
    dataTable{kk+1,5} = LsDiff(2);
    dataTable{kk+1,6} = LsDiff(3);
    dataTable{kk+1,7} = LsDiff(1)*rNorm;
    dataTable{kk+1,8} = LsDiff(2)*rNorm;
    dataTable{kk+1,9} = LsDiff(3)*rNorm;
    if isfield(secondary,'J2') == 0 || secondary.J2 == 0
        dataTable{kk+1,3} = 'No';
    else
        dataTable{kk+1,3} = 'Yes';
    end
end
dataTable
end

















