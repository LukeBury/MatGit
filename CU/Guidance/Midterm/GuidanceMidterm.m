clear
clc
close all
mbinPath = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/mbin';
figSavePath = '/Users/lukebury/CU_Google_Drive/Documents/School/CU/Courses/5-6519-Space_Vehicle_Guidance_and_Control/Midterm/Midterm_LaTex/Figures/';
addpath(genpath(mbinPath))
tic
% ========================================================================
%%% Run/Plot Switches
% ========================================================================

% ========================================================================
%%% Importing Data
% ========================================================================
%%% General data on solar system bodies
bodies = getBodyData(mbinPath);

%%% Color options/schemes
colors = get_colors();

% ========================================================================
%%% 1 - Natural Dynamics
% ========================================================================
% -------------------------------------------------
% System parameters
% -------------------------------------------------
earth = bodies.earth;

%%% parameters
earth.u = 398600.4418; % km^3/s^2

ISS.a = 6728; % km

%%% Initial conditions
x0  = -1.5; % km
y0  = -0.5; % km
dx0 = 3e-3; % km/s
dy0 = 0;    % km/s

%%% Final conditions
xf  = -250e-3; % km
yf  = 0;       % km
dxf = 0.2e-3;  % km/s
dyf = 0;       % km/s

%%% Time vector for integration
t0 = 0;     % sec
dt = .01;     % sec
tf = 20*60; % sec
timeVec = t0:dt:tf; % sec

%%% Choosing ode45 tolerance
tol = 1e-10;

%%% Setting integrator options
options = odeset('RelTol',tol,'AbsTol',tol);

% -------------------------------------------------
% Intregrating states with no control and plotting results
% -------------------------------------------------
%%% Calculating mean motion
ISS.n = sqrt(earth.u/(ISS.a^3));

%%% Setting initial and final state
X0 = [x0; y0; dx0; dy0];
Xf = [xf; yf; dxf; dyf];

%%% Integrating
[t_out, X_nc] = ode45(@Int_GuidanceMidterm_State_noControl, timeVec, X0, options, ISS.n);

%%% Looking at minimum distance to origin w/ natural dynamics
rNorm_nc = rowNorm(X_nc(:,1:2));
rNorm_nc_min = min(rNorm_nc);
rNorm_nc_minIndx = find(rNorm_nc == rNorm_nc_min);


figure('position',[440 464 803 334])
subplot(1,2,1); hold all
p1 = plot(X_nc(:,1),X_nc(:,2),'b','linewidth',2);
plot(X_nc(1,1),X_nc(1,2),'ro','linewidth',1.5,'markersize',10)
plot(X_nc(end,1),X_nc(end,2),'rx','linewidth',1.5,'markersize',10)
plot(X_nc(rNorm_nc_minIndx,1),X_nc(rNorm_nc_minIndx,2),'ko','linewidth',1.5,'markersize',10,'markerfacecolor','c')
PlotBoi2('x, $km$','y, $km$',16,'LaTex')
legend([p1], 'Nominal Uncontrolled','location','best')
%
subplot(1,2,2); hold all
plot(timeVec./60, rNorm_nc,'b','linewidth',2)
p1 = plot(timeVec(1)./60,rNorm_nc(1),'ro','linewidth',1.5,'markersize',10);
p2 = plot(timeVec(end)./60,rNorm_nc(end),'rx','linewidth',1.5,'markersize',10);
p3 = plot(timeVec(rNorm_nc_minIndx)./60,rNorm_nc_min,'ko','linewidth',1.5,'markersize',10,'markerfacecolor','c');
legend([p1, p2, p3],'r_0','r_f','r_{min}','location','best')
PlotBoi2('Time, $min$','Distance from ISS (origin), $km$',16,'LaTex')
%%% Saving figure
saveas(gcf,[figSavePath,'NaturalDynamicsProp.png'])


% ========================================================================
%%% 2 - Optimal Trajectory
% ========================================================================
%%% Creating A and B matrices
A = [0          0   1         0;...
     0          0   0         1;...
     3*ISS.n^2  0   0         2*ISS.n;...
     0          0   -2*ISS.n  0];

B = [0 0;...
     0 0;...
     1 0;...
     0 1];


%%% Creating backwards time vector for propagating backwards
timeVec_bkwd = fliplr(timeVec);

%%% Initializing K with its final conditions
K_tf = eye(4);
K_tf = reshape(K_tf,16,1);
% K_tf = zeros(16,1);

%%% Integrating K
[t_bkwds, K_bkwds] = ode45(@Int_GuidanceMidterm_K, timeVec_bkwd, K_tf, options, A, B);
K = flipud(K_bkwds);

%%% Integrating s
s_tf = -Xf; % s_tf = -S*ref;
[t_bkwds, s_bkwd] = ode45(@Int_GuidanceMidterm_s, timeVec_bkwd,s_tf,options,K,A,B,Xf,timeVec);
s = flipud(s_bkwd);

%%% Integrating controlled state
[t_out, X_ctrl] = ode45(@Int_GuidanceMidterm_State_optControl, timeVec, X0, options, B, K, s, ISS.n, timeVec);

%%% Learning what the controls were
us = zeros(size(X_ctrl,1),2);
J = zeros(size(X_ctrl,1),1);
for kk = 1:size(X_ctrl,1)
    K_now = reshape(K(kk,:),4,4);
    s_now = s(kk,:)';

    %%% Setting controls
    us(kk,:) = -B'*K_now*X_ctrl(kk,:)' - B'*s_now;
    
    %%% Looking at cost
    J(kk) = us(kk,:) * us(kk,:)';
end
J = trapz(timeVec,J);
J


figure('position',[440 56 803 334]);
subplot(1,2,1); hold all
p1 = plot(X_ctrl(:,1),X_ctrl(:,2),'b','linewidth',2);
plot(X_ctrl(1,1),X_ctrl(1,2),'ro','linewidth',1.5,'markersize',10);
p3 = plot(X_nc(:,1),X_nc(:,2),'--b','linewidth',2);
p2 = plot(Xf(1),Xf(2),'k^','linewidth',1.5,'markersize',10,'markerfacecolor','g','markeredgecolor','k');
legend([p3 p1 p2],'Nominal Uncontrolled','Nominal Controlled','Target','location','best')
PlotBoi2('x, $km$','y, $km$',16,'LaTex')
subplot(1,2,2); hold all
p1 = plot(timeVec,us(:,1),'r','linewidth',2);
p2 = plot(timeVec,us(:,2),'k','linewidth',2);
legend([p1 p2],'u_x','u_y','location','best')
PlotBoi2('Time, $sec$','Control, $km/s^2$',16,'LaTex')
%%% Saving figure
saveas(gcf,[figSavePath,'ControlledTrajectory.png'])


% -------------------------------------------------
% Testing something
% -------------------------------------------------
%%
lambda_true = zeros(size(us,1),4);
for kk = 1:size(us,1)
    lambda_true(kk,:) = (-B')\us(kk,:)';
    
    
end
% 
% lam_tf = -Xf;
% % lam_tf = reshape(K(end,:),4,4)*Xf + s(end,:)';
% %%% Integrating
% [t_out_bkwd, lambda_bkwd] = ode45(@Int_GuidanceMidterm_lambda, timeVec_bkwd, lam_tf, options, A);
% lambda_test = flipud(lambda_bkwd);
% figure
% plot(lambda_test)
% lam_t0 = reshape(K(1,:),4,4)*Xf;
% [t_out, lambda1] = ode45(@Int_GuidanceMidterm_lambda, timeVec, lam_t0, options, A);
% 
% lam_t0 = reshape(K(1,:),4,4)*X0;
% [t_out, lambda2] = ode45(@Int_GuidanceMidterm_lambda, timeVec, lam_t0, options, A);
% 
% lam_t0 = -Xf*1e-6;
% [t_out, lambda3] = ode45(@Int_GuidanceMidterm_lambda, timeVec, lam_t0, options, A);
% 
% 
% 
% 
% [t_out, X_text] = ode45(@Int_GuidanceMidterm_test, timeVec, X0, options, B, lambda_test, ISS.n, timeVec);
% 
% 
% figure; hold all
% plot(X_text(:,1),X_text(:,2),'b')
% plot(X_ctrl(:,1),X_ctrl(:,2),'r')

% -------------------------------------------------
% What happens if you get to the target then coast? - 2d
% -------------------------------------------------
X0_coast = X_ctrl(end,:)';
%%% Integrating
[t_out, X_coast] = ode45(@Int_GuidanceMidterm_State_noControl, timeVec, X0_coast, options, ISS.n);

figure('position',[440 464 433 334]); hold all
plot(X_ctrl(1,1),X_ctrl(1,2),'ro','linewidth',1.5,'markersize',10);
p1 = plot(X_ctrl(:,1),X_ctrl(:,2),'b','linewidth',2);
p2 = plot(Xf(1),Xf(2),'^','linewidth',1.5,'markersize',10,'markerfacecolor','g','markeredgecolor','k');
p3 = plot(X_coast(:,1),X_coast(:,2),'m','linewidth',2);
plot(X_coast(end,1),X_coast(end,2),'rx','linewidth',1.5,'markersize',10);
legend([p1 p2 p3],'Nominal Controlled','Target','Coasting','location','best')
PlotBoi2('x, $km$','y, $km$',16,'LaTex')
%%% Saving figure
saveas(gcf,[figSavePath,'CoastTrajectory.png'])


% ========================================================================
%%% 3 - Guidance
% ========================================================================
% -------------------------------------------------
% Using the same control - 3a
% -------------------------------------------------
%%% Creating perturbed state
x0_p  = -0.98;  % km
y0_p  = -0.4;   % km
dx0_p = 3.7e-3; % km/s
dy0_p = 0.3e-3; % km/s
X0_p = [x0_p; y0_p; dx0_p; dy0_p];

%%% Integrating
[t_out, X_p] = ode45(@Int_GuidanceMidterm_State_uGiven, timeVec, X0_p, options, us, ISS.n, timeVec);

figure('position',[440 464 433 334]); hold all
plot(X_ctrl(1,1),X_ctrl(1,2),'ro','linewidth',1.5,'markersize',10);
p1 = plot(X_ctrl(:,1),X_ctrl(:,2),'b','linewidth',2);
plot(Xf(1),Xf(2),'^','linewidth',1.5,'markersize',10,'markerfacecolor','g','markeredgecolor','k');
p3 = plot(X_p(:,1),X_p(:,2),'--b','linewidth',2);
plot(X_p(1,1),X_p(1,2),'ro','linewidth',1.5,'markersize',10)
plot(X_p(end,1),X_p(end,2),'rx','linewidth',1.5,'markersize',10)
legend([p1 p3],'Nominal Controlled','Perturbed, old controls','location','best')
PlotBoi2('x, $km$','y, $km$',16,'LaTex')
%%% Saving figure
saveas(gcf,[figSavePath,'PerturbedTrajectory.png'])


% -------------------------------------------------
% Neighboring State - 3b
% -------------------------------------------------

%%% Integrating
[t_out, X_nbr] = ode45(@Int_GuidanceMidterm_neighboringState, timeVec, X0_p,options,B,K,X_ctrl,us,ISS.n,timeVec);

%%% Calculating ranges of each trajectory from target over time
range_ctrl = rowNorm(X_ctrl(:,1:2) - Xf(1:2)');
range_nbr  = rowNorm(X_nbr(:,1:2) - Xf(1:2)');
vel_ctrl = rowNorm(X_ctrl(:,3:4) - Xf(3:4)');
vel_nbr  = rowNorm(X_nbr(:,3:4) - Xf(3:4)');

figure('position',[440 464 803 334])
subplot(2,2,[1 3]); hold all
plot(X_ctrl(1,1),X_ctrl(1,2),'ro','linewidth',1.5,'markersize',10);
p1 = plot(X_ctrl(:,1),X_ctrl(:,2),'b','linewidth',2);
plot(Xf(1),Xf(2),'^','linewidth',1.5,'markersize',10,'markerfacecolor','g','markeredgecolor','k');
p2 = plot(X_nbr(:,1),X_nbr(:,2),'--b','linewidth',2);
plot(X_nbr(1,1),X_nbr(1,2),'ro','linewidth',1.5,'markersize',10)
PlotBoi2('x, $km$','y, $km$',16,'LaTex')
legend([p1 p2],'Nominal Controlled','Perturbed, new controls','location','best')
subplot(2,2,2); hold all
plot(timeVec./60, range_ctrl,'b','linewidth',1.5);
plot(timeVec./60, range_nbr, '--b','linewidth',1.5);
PlotBoi2('Time, $min$','Target position error, $km$',12,'LaTex')
subplot(2,2,4); hold all
plot(timeVec./60, vel_ctrl,'b','linewidth',1.5);
plot(timeVec./60, vel_nbr, '--b','linewidth',1.5);
PlotBoi2('Time, $min$','Target velocity error, $km/s$',12,'LaTex')
%%% Saving figure
saveas(gcf,[figSavePath,'NeighborTrajectory.png'])














% ========================================================================
%%% Functions
% ========================================================================

function [dX] = Int_GuidanceMidterm_State_noControl(t,X,n)
%%% Preallocating
dX = zeros(4,1);

%%% Unpacking
x  = X(1);
dx = X(3);
dy = X(4);

%%% Setting dynamics
dX(1:2,1) = [dx;dy];
dX(3) = 2*n*dy + 3*n*n*x;
dX(4) = -2*n*dx;

end


function [dK] = Int_GuidanceMidterm_K(t,K,A,B)
K = reshape(K,4,4);
dK = -A'*K - K*A + K*B*(B'*K);
dK = reshape(dK,16,1);
end


function [ds] = Int_GuidanceMidterm_s(t,s,K,A,B,ref,fullTime)
%%% Interpollating to get current K
K_now = interp1(fullTime,K,t);
K = reshape(K_now,4,4);

ds = -(A' - (K * B)*B')*s;
end


% function [ds] = Int_GuidanceHW2P2_LittleS(t,s,A,B,Q,R,S,fullTime)
% S_now = interp1(fullTime,S,t);
% S = reshape(S_now,5,5);
% 
% ref = [0; 0; sin(pi*t/5); 0; 0];
% 
% ds = -[A' - S * B * inv(R) *B']*s + Q*ref;
% 
% end



function [dX] = Int_GuidanceMidterm_State_optControl(t,X,B,K,s,n,fullTime)
%%% Preallocating
dX = zeros(4,1);

%%% Unpacking
x  = X(1);

dx = X(3);
dy = X(4);

%%% Interpollating to get current K
K_now = interp1(fullTime,K,t);
K = reshape(K_now,4,4);

%%% Interpollating to get current s
s_now = interp1(fullTime,s,t);
s_now = s_now';

%%% Setting controls
u = -B'*K*X - B'*s_now;

%%% Setting dynamics
dX(1:2,1) = X(3:4);
dX(3) = 2*n*dy + 3*n*n*x + u(1);
dX(4) = -2*n*dx + u(2);

end



function [dLambda] = Int_GuidanceMidterm_lambda(t,Lambda,A)
dLambda = -A' * Lambda;
end


function [dX] = Int_GuidanceMidterm_State_uGiven(t,X,u,n,fullTime)
%%% Preallocating
dX = zeros(4,1);

%%% Unpacking
x  = X(1);

dx = X(3);
dy = X(4);

%%% Interpollating to get current u
u = interp1(fullTime,u,t);
u = u';

%%% Setting dynamics
dX(1:2,1) = X(3:4);
dX(3) = 2*n*dy + 3*n*n*x + u(1);
dX(4) = -2*n*dx + u(2);

end


function [dX] = Int_GuidanceMidterm_neighboringState(t,X,B,K,XRef,uRef,n,fullTime)
%%% Preallocating
dX = zeros(4,1);

%%% Unpacking
x  = X(1);

dx = X(3);
dy = X(4);

%%% Interpollating to get current K
K_now = interp1(fullTime,K,t);
K = reshape(K_now,4,4);

%%% Interpollating to get current reference state
XRef_now = interp1(fullTime,XRef,t);
XRef = XRef_now';

%%% Interpollating to get current u
uRef = interp1(fullTime,uRef,t);
uRef = uRef';

%%% Determining neighboring u component
du = -(B'*K)*(X - XRef);

%%% Determining total u
u = uRef + du;

%%% Setting dynamics
dX(1:2,1) = X(3:4);
dX(3) = 2*n*dy + 3*n*n*x + u(1);
dX(4) = -2*n*dx + u(2);

end


function [dX] = Int_GuidanceMidterm_test(t,X,B,lambdas,n,fullTime)
%%% Preallocating
dX = zeros(4,1);

%%% Unpacking
x  = X(1);

dx = X(3);
dy = X(4);


%%% Interpollating to get current lambda
lam_now = interp1(fullTime,lambdas,t);
lam = lam_now';

%%% Interpollating to get current u
u = -B'*lam;

%%% Setting dynamics
dX(1:2,1) = X(3:4);
dX(3) = 2*n*dy + 3*n*n*x + u(1);
dX(4) = -2*n*dx + u(2);
end









