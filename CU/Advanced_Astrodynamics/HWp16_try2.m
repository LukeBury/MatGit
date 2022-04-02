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
primary = bodies.earth;
secondary = bodies.moon;

%%% Normalizing constants
rNorm = secondary.a;         % n <-> km
tNorm = 1/secondary.meanMot; % n <-> sec
vNorm = rNorm / tNorm;       % n <-> km/sec

%%% Choosing ode tolerance
tol = 1e-13;

%%% Setting integrator options
options = odeset('RelTol',tol,'AbsTol',tol);

%%% Setting necessary parameters for integration
prms.u = secondary.MR;

% ========================================================================
%%% Integrating in Rotating Frame
% ========================================================================
% ------------------------------------
%%% Preparing for integration
% ------------------------------------
%%% Initial conditions
r0_BCRn = [1-secondary.MR - 100*secondary.R_n; 0; 0];
v0_BCRn = [0; 0.6; 0];
X0_BCRn = [r0_BCRn; v0_BCRn];

%%% Setting time
t_i_n = 0; % sec
t_f_n = 2*pi; % Long bc events are watching for impact or escape
n_dt = 1000;
time0_n = linspace(t_i_n,t_f_n,n_dt);

% ------------------------------------
%%% Integrating
% ------------------------------------
%%% Integrating
[time_n, X_BCR_n] = ode113(@Int_CR3Bn, time0_n, X0_BCRn, options, prms);

% ------------------------------------
%%% Build rotation matrix
% ------------------------------------
%%% Setting up rotation matrix
T = @(t) [cos(t), -sin(t), 0;...
          sin(t), cos(t),  0;...
          0,      0,       1];
T0 = T(0);

% ------------------------------------
%%% Post-processing
% ------------------------------------
%%% Checking Jacobi constant
[JC_sc] = JacobiConstantCalculator(secondary.MR,X_BCR_n(:,1:3),X_BCR_n(:,4:6));

figure
plot(time_n, percentchange(JC_sc),'b.')
PlotBoi2('Time','\% Change, $J_C$',18,'LaTex')


[ r_BCIn ] = r_BCR2BCI( X_BCR_n(:,1:3), time_n, 1 );
[ v_BCIn ] = v_BCR2BCI( X_BCR_n(:,4:6), r_BCIn, time_n, 1 );
X_BCIn_rotating = [r_BCIn, v_BCIn];
% ------------------------------------
%%% Plotting
% ------------------------------------
%%% Plotting in rotating frame
figure; hold all
plot3(X_BCR_n(:,1),X_BCR_n(:,2),X_BCR_n(:,3),'b','linewidth',1.5)
plotBody2(primary.R/rNorm,[-secondary.MR,0,0],colors.std.blue,colors.std.black,1.5,1)
plotBody2(secondary.R_n,[1-secondary.MR,0,0],colors.std.grey,colors.std.black,1.5,1)
axis equal
PlotBoi2('$X_n$','$Y_n$',18,'LaTex')
title('Trajectory in normalized CR3BP')

figure; hold all
% plot3(X_BCI_convert(:,1),X_BCI_convert(:,2),X_BCI_convert(:,3),'g','linewidth',1.5)
plot3(X_BCIn_rotating(:,1),X_BCIn_rotating(:,2),X_BCIn_rotating(:,3),'b','linewidth',1.5)
plotBody2(primary.R/rNorm,[-secondary.MR,0,0],colors.std.blue,colors.std.black,1.5,1)
plotBody2(secondary.R_n,[1-secondary.MR,0,0],colors.std.grey,colors.std.black,1.5,1)
axis equal
PlotBoi2('$X_n$','$Y_n$',18,'LaTex')
title('CR3BP trajectory rotated into inertial frame')
% ========================================================================
%%% Integrating in Rotating Frame
% ========================================================================
% ------------------------------------
%%% Converting initial conditions to inertial non-normalized
% ------------------------------------
r0_BCIn = T0 * (r0_BCRn);
v0_BCIn = T0 * ([v0_BCRn(1) - r0_BCRn(2);...
               v0_BCRn(2) + r0_BCRn(1);...
               v0_BCRn(3)]);
X0_BCIn = [r0_BCIn; v0_BCIn];

% ------------------------------------
%%% Integrating
% ------------------------------------
%%% Integrating
[time_n, X_BCIn_inertial] = ode113(@Int_3BI, time0_n, X0_BCIn, options, secondary.MR);

% ------------------------------------
%%% Plotting
% ------------------------------------
%%% Plotting in rotating frame
figure; hold all
plot3(X_BCIn_inertial(:,1),X_BCIn_inertial(:,2),X_BCIn_inertial(:,3),'c','linewidth',1.5)
plotBody2(primary.R/rNorm,[-secondary.MR,0,0],colors.std.blue,colors.std.black,1.5,1)
plotBody2(secondary.R_n,[1-secondary.MR,0,0],colors.std.grey,colors.std.black,1.5,1)
axis equal
PlotBoi2('$X_n$','$Y_n$',18,'LaTex')
title('Trajectory integrated with inertial dynamics')

% ========================================================================
%%% Finding difference
% ========================================================================
diff_X = X_BCIn_inertial - X_BCIn_rotating;
% figure
% subplot(3,2,1)
% plot(time./86400,diff_X(:,1),'r','linewidth',2)
% PlotBoi2('','$\Delta x$',18,'LaTex')
% subplot(3,2,3)
% plot(time./86400,diff_X(:,2),'r','linewidth',2)
% PlotBoi2('','$\Delta y$',18,'LaTex')
% subplot(3,2,5)
% plot(time./86400,diff_X(:,3),'r','linewidth',2)
% PlotBoi2('Time, days','$\Delta z$',18,'LaTex')
% subplot(3,2,2)
% plot(time./86400,diff_X(:,4),'r','linewidth',2)
% PlotBoi2('','$\Delta \dot{x}$',18,'LaTex')
% subplot(3,2,4)
% plot(time./86400,diff_X(:,5),'r','linewidth',2)
% PlotBoi2('','$\Delta \dot{y}$',18,'LaTex')
% subplot(3,2,6)
% plot(time./86400,diff_X(:,6),'r','linewidth',2)
% PlotBoi2('Time, days','$\Delta \dot{z}$',18,'LaTex')
figure
subplot(2,2,1)
plot(time_n./8640,diff_X(:,1),'r','linewidth',2)
PlotBoi2('','$\Delta x$',18,'LaTex')
subplot(2,2,3)
plot(time_n,diff_X(:,2),'r','linewidth',2)
PlotBoi2('Time','$\Delta y$',18,'LaTex')
subplot(2,2,2)
plot(time_n,diff_X(:,4),'r','linewidth',2)
PlotBoi2('','$\Delta \dot{x}$',18,'LaTex')
subplot(2,2,4)
plot(time_n,diff_X(:,5),'r','linewidth',2)
PlotBoi2('Time','$\Delta \dot{y}$',18,'LaTex')


% ========================================================================
%%% 16-c, Looking at orbital elements
% ========================================================================
as        = zeros(size(X_BCIn_rotating,1),1);
es        = zeros(size(X_BCIn_rotating,1),1);
% is_deg    = zeros(size(X_BCI_rotating,1),1);
% raans_deg = zeros(size(X_BCI_rotating,1),1);
ws_deg    = zeros(size(X_BCIn_rotating,1),1);
tas_deg   = zeros(size(X_BCIn_rotating,1),1);

for kk = 1:size(X_BCIn_rotating,1)
    [a,e,i_deg,raan_deg,w_deg,ta_deg] = RV2COE(X_BCIn_rotating(kk,1:3),X_BCIn_rotating(kk,4:6),1);
    
    as(kk)        = a; 
    es(kk)        = e;
%     is_deg(kk)    = i_deg;
%     raans_deg(kk) = raan_deg;
    ws_deg(kk)    = w_deg;
    tas_deg(kk)   = ta_deg;
end

% figure
% plot(as)
% figure
% plot(es)
% % figure
% % plot(is_deg)
% % figure
% % plot(raans_deg)
% figure
% plot(ws_deg)
% figure
% plot(tas_deg)

figure('position',[440 368 706 430])
subplot(2,2,1)
plot(time_n,as,'b.')
PlotBoi2('','$a$, km', 18, 'LaTex')
xlim([0 time_n(end)])
subplot(2,2,3)
plot(time_n,es,'b.')
PlotBoi2('Time','$e$', 18, 'LaTex')
xlim([0 time_n(end)])
subplot(2,2,2)
plot(time_n,ws_deg,'b.')
PlotBoi2('','$\omega$, $^\circ$', 18, 'LaTex')
xlim([0 time_n(end)])
subplot(2,2,4)
plot(time_n,tas_deg,'b.')
PlotBoi2('Time','$f$, $^\circ$', 18, 'LaTex')
xlim([0 time_n(end)])


% ========================================================================
%%% 16-c, Looking at orbital elements for 2B motion
% ========================================================================
% % prms.u = 0;
% % [time_n, X2_BCR_n] = ode113(@Int_CR3Bn, time0_n, X0_BCRn, options, prms);
% % %%% Converting X_BCR_n to X_BCI
% % X2_BCR = [X2_BCR_n(:,1:3).*rNorm, X2_BCR_n(:,4:6).*vNorm];
% % [ r2_BCI ] = r_BCR2BCI( X2_BCR(:,1:3), time_n.*tNorm, secondary.meanMot );
% % [ v2_BCI ] = v_BCR2BCI( X2_BCR(:,4:6), r2_BCI, time_n.*tNorm, secondary.meanMot );
% % X2_BCI_rotating = [r2_BCI, v2_BCI];

[time, X2_BCI_inertial] = ode113(@Int_2BI, time0_n.*tNorm, [X0_BCIn(1:3).*rNorm; X0_BCIn(4:6).*vNorm], options, primary.u);

%%% Plotting in rotating frame
figure; hold all
plot3(X2_BCI_inertial(:,1)./rNorm,X2_BCI_inertial(:,2)./rNorm,X2_BCI_inertial(:,3)./rNorm,'m','linewidth',1.5)
plotBody2(primary.R/rNorm,[0,0,0],colors.std.blue,colors.std.black,1.5,1)
axis equal
PlotBoi2('$X_n$','$Y_n$',18,'LaTex')
title('Trajectory integrated with 2-body dynamics around Earth')

as2        = zeros(size(X2_BCI_inertial,1),1);
es2        = zeros(size(X2_BCI_inertial,1),1);
% is_deg    = zeros(size(X2_BCI_rotating,1),1);
% raans_deg = zeros(size(X2_BCI_rotating,1),1);
ws2_deg    = zeros(size(X2_BCI_inertial,1),1);
tas2_deg   = zeros(size(X2_BCI_inertial,1),1);

for kk = 1:size(X2_BCI_inertial,1)
    [a,e,i_deg,raan_deg,w_deg,ta_deg] = RV2COE(X2_BCI_inertial(kk,1:3),X2_BCI_inertial(kk,4:6),(primary.u));
    
    as2(kk)        = a; 
    es2(kk)        = e;
%     is_deg(kk)    = i_deg;
%     raans_deg(kk) = raan_deg;
    ws2_deg(kk)    = w_deg;
    tas2_deg(kk)   = ta_deg;
end

figure('position',[440 368 706 430])
subplot(2,2,1)
plot(time_n,percentchange(as2),'m.')
PlotBoi2('','\% Change - $a$', 18, 'LaTex')
xlim([0 time_n(end)])
subplot(2,2,3)
plot(time_n,percentchange(es2),'m.')
PlotBoi2('Time','\% Change - $e$', 18, 'LaTex')
xlim([0 time_n(end)])
subplot(2,2,2)
plot(time_n,percentchange(ws2_deg),'m.')
PlotBoi2('','\% Change - $\omega$', 18, 'LaTex')
xlim([0 time_n(end)])
subplot(2,2,4)
plot(time_n,tas2_deg,'m.')
PlotBoi2('Time','$f$, $^\circ$', 18, 'LaTex')
xlim([0 time_n(end)])



% ========================================================================
% ========================================================================
%%% Functions
% ========================================================================
% ========================================================================
function [dX] = Int_3BI(t,X,mu)
%%% For numerical integration of ECI state under influence of standard
%%% 3-body gravity
%%% Inputs:
%          t - time vector
%          X - initial state [6x1]
%          u - gravitational parameter of primary body
% ------------------------------------
%%% Preallocate
% ------------------------------------
%%% Preallocate state output
dX = zeros(6,1);

% ------------------------------------
%%% Compute positions of bodies
% ------------------------------------
r = [X(1); X(2); X(3)]; % km
% 
% r_B1 = R3([-mu*a;    0; 0], t*n);
% r_B2 = R3([(1-mu)*a; 0; 0], t*n);
r_B1 = [-mu*cos(t); -mu*sin(t); 0];
r_B2 = [(1-mu)*cos(t); (1-mu)*sin(t); 0];

% ------------------------------------
%%% Find relative positions
% ------------------------------------
r1 = r - r_B1;
r1_mag = norm(r1);
r2 = r - r_B2;
r2_mag = norm(r2);

% ------------------------------------
%%% EOM
% ------------------------------------
%%% 3B accelerations
% a_3B = -u_B1*r1 / (r1_mag^3) - u_B2*r2 / (r2_mag^3);
a_3B = -(1-mu)*r1 / (r1_mag^3) - mu*r2 / (r2_mag^3);

%%% Output the derivative of the state
dX(1:3) = [X(4); X(5); X(6)];    % km/s
dX(4:6) = [a_3B(1); a_3B(2); a_3B(3)]; % km/s^2

end



