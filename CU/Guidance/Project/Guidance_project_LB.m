clear
clc
close all
mbinPath = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/mbin';
addpath(genpath(mbinPath))
tic

% ========================================================================
%%% Run Switches
% ========================================================================

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
% -------------------------------
% System info
% -------------------------------
%%% bodies
primary = bodies.jupiter;
secondary = bodies.europa;

% primary = bodies.saturn;
% secondary = bodies.enceladus;

%%% Normalizing constants
rNorm = secondary.a;         % n <-> km
tNorm = 1/secondary.meanMot; % n <-> sec
vNorm = rNorm / tNorm;       % n <-> km/sec

%%% Lagrange Points
L123 = EquilibriumPoints(secondary.MR,1:3); % [3x3] of L1, L2, L3 normalized BCR coordinates


% -------------------------------
% Creating Rerence Trajectory and Choosing Perturbation
% -------------------------------
%%% Setting time vector
t_0 = 0;
% t_f = 4*pi;
% t_f = 2.677618658082257; % Europa impact
% t_f = 2.287449542144822; % Enceladus impact
t_f = 5;
n_dt = 5000;
time0_n = linspace(t_0,t_f,n_dt);
dt = time0_n(2)-time0_n(1);

% dt = 0.001;
% time0_n = [0:dt:5]';

%%% To impact Europa plume site at low angle
r0ref_n = [1.0204617015266166, -0.0054921659260262, -0.0099525155929432]';
v0ref_n = [-0.0111591144905944, -0.0104200955172558, 0.0145837924677663]';

% %%% Enceladus south pole
% r0ref_n = [1.0039914454138370, 0.0006268291981895, 0.0032971880321676]';
% v0ref_n = [-0.0012896112019572, -0.0043048383232825, -0.0030710153996081]';

%%% Perturbation
dX0 = [0; 100/rNorm;0;-0.05/vNorm;0;0];
% dX0 = [10/rNorm; 0;0;0;0;0];

%%% Setting initial state of reference
X0ref_n = [r0ref_n; v0ref_n];

%%% extra parameters
prms = struct();
prms.u = secondary.MR;
prms.R2_n = secondary.R_n;
prms.L1x = L123(1,1);
prms.L2x = L123(2,1);

%%% Choosing ode tolerance
% tol = 2.22045e-14;
tol = 1e-12;

%%% Setting Options
options = odeset('RelTol',tol,'AbsTol',tol);
options_ImpactEscape = odeset('Events',@event_ImpactEscape_CR3Bn,'RelTol',tol,'AbsTol',tol);

%%% Integrating reference
[time_n, Xref_n, time_eventImpact, X_eventImpact, index_eventImpact] = ode113(@Int_CR3Bn,...
            time0_n, X0ref_n, options_ImpactEscape, prms);
 
%%% Creating a backwards time vector for integration of K
time_bkwds_n = flipud(time_n);

% -------------------------------
% Creating B matrix and selecting weights for R, Q, and S
% -------------------------------
%%% B matrix (control inputs)
B = [zeros(3,3); eye(3)];

%%% Weights
% a_rf = 1000; % Final position error weights (S,1-3)
% a_vf = 1000; % Final velocity error weights (S,4-6)
% a_r  = 1; % Position error weights (Q,1-3)
% a_v  = 1; % Velocity error weights (Q,4-6)
% a_u  = 1e-6; % Control weights (R)

a_rf = 1; % Final position error weights (S,1-3)
a_vf = 1; % Final velocity error weights (S,4-6)
a_r  = 1000; % Position error weights (Q,1-3)
a_v  = 10; % Velocity error weights (Q,4-6)
a_u  = 1; % Control weights (R)

%%% Weights on controls (u' R u)
R = eye(3).*a_u;
%%% Weights on states (x' Q x)
Q = blkdiag(eye(3).*a_r, eye(3).*a_v);
%%% Weights on final states (xf' S xf
S = blkdiag(eye(3).*a_rf, eye(3).*a_vf);


% Final Position Error: 210.8 meters
% Final Velocity Error: 0.190 m/s
% Total delta-V: 1.665 m/s
% Elapsed time is 45.552687 seconds.
% -------------------------------
% Creating time-history of S
% -------------------------------

%%% Choosing initial (final) conditions
K_tf = S;
K_tf = reshape(K_tf,6,6);

%%% Integrating
[t_bkwds, K_out_bkwds] = ode45(@Int_K, time_bkwds_n, K_tf, options,B,R,Q,Xref_n,time_n,secondary.MR);
K_out_fwds = flipud(K_out_bkwds);

% -------------------------------
% Integrating state with guidance algorithm
% -------------------------------

%%% Creating perturbed X0
X0_n = X0ref_n + dX0;

%%% Integrating
[t_out, X_guidance] = ode45(@Int_State, time_n,X0_n,options,B,R,K_out_fwds,Xref_n,time_n,secondary.MR,dt);

% -------------------------------
% Recreating controls
% -------------------------------
%%% Preallocating
us_mpss = zeros(length(t_out)-1,3);
for kk = 1:(length(t_out)-1)
    %%% Current K matrix
    K = reshape(K_out_fwds(kk,:),6,6);
    
    %%% Computing Controls
    u_t = -(inv(R)*B'*K)*(X_guidance(kk,:)' - Xref_n(kk,:)');
    
    %%% In meters/second^2
    us_mpss(kk,:) = u_t' .*(rNorm/(tNorm^2))*1000;
end


figure('position',[723 141 560 420])
subplot(3,1,1)
plot(t_out(1:end-1).*tNorm./3600,us_mpss(:,1),'linewidth',2)
PlotBoi2('','$u_x$, $m/s^2$',16,'LaTex')
subplot(3,1,2)
plot(t_out(1:end-1).*tNorm./3600,us_mpss(:,2),'linewidth',2)
PlotBoi2('','$u_y$, $m/s^2$',16,'LaTex')
subplot(3,1,3)
plot(t_out(1:end-1).*tNorm./3600,us_mpss(:,3),'linewidth',2)
PlotBoi2('Time, hr','$u_z$, $m/s^2$',16,'LaTex')

figure('position',[155 137 560 420]); hold all
plot3(Xref_n(:,1),Xref_n(:,2),Xref_n(:,3),'b','linewidth',1.5)
plot3(X_guidance(:,1),X_guidance(:,2),X_guidance(:,3),'r','linewidth',1.5)
plotBodyTexture3(secondary.R_n,[1-secondary.MR,0,0],secondary.img)
PlotBoi3('$x_n$','$y_n$','$z_n$',16,'LaTex')
legend('Reference','Perturbed')
axis equal

figure
subplot(4,2,1)
plot(t_out.*tNorm./3600, (X_guidance(:,1)-Xref_n(:,1)).*rNorm,'r','linewidth',2)
PlotBoi2('','$\Delta x$, $km$',16,'LaTex')
subplot(4,2,3)
plot(t_out.*tNorm./3600, (X_guidance(:,2)-Xref_n(:,2)).*rNorm,'r','linewidth',2)
PlotBoi2('','$\Delta y$, $km$',16,'LaTex')
subplot(4,2,5)
plot(t_out.*tNorm./3600, (X_guidance(:,3)-Xref_n(:,3)).*rNorm,'r','linewidth',2)
PlotBoi2('','$\Delta z$, $km$',16,'LaTex')
subplot(4,2,2)
plot(t_out.*tNorm./3600, (X_guidance(:,4)-Xref_n(:,4)).*vNorm*1000,'r','linewidth',2)
PlotBoi2('','$\Delta \dot{x}$, $m/s$',16,'LaTex')
subplot(4,2,4)
plot(t_out.*tNorm./3600, (X_guidance(:,5)-Xref_n(:,5)).*vNorm*1000,'r','linewidth',2)
PlotBoi2('','$\Delta \dot{y}$, $m/s$',16,'LaTex')
subplot(4,2,6)
plot(t_out.*tNorm./3600, (X_guidance(:,6)-Xref_n(:,6)).*vNorm*1000,'r','linewidth',2)
PlotBoi2('','$\Delta \dot{z}$, $m/s$',16,'LaTex')
subplot(4,2,7)
plot(t_out.*tNorm./3600, (rowNorm(X_guidance(:,1:3) - [1-secondary.MR,0,0])-secondary.R_n).*rNorm,'b','linewidth',2)
PlotBoi2('Time, hr','Altitude, $km$',14,'LaTex')
ylim([0 max((rowNorm(X_guidance(:,1:3) - [1-secondary.MR,0,0])-secondary.R_n).*rNorm)])
subplot(4,2,8)
plot(t_out.*tNorm./3600, (rowNorm(X_guidance(:,1:3) - [1-secondary.MR,0,0])-secondary.R_n).*rNorm,'b','linewidth',2)
PlotBoi2('Time, hr','Altitude, $km$',14,'LaTex')
ylim([0 max((rowNorm(X_guidance(:,1:3) - [1-secondary.MR,0,0])-secondary.R_n).*rNorm)])


fprintf('Final Position Error: %1.1f meters\n', norm(X_guidance(end,1:3) - Xref_n(end,1:3)).*rNorm*1000)
fprintf('Final Velocity Error: %1.3f m/s\n',norm(X_guidance(end,4:6) - Xref_n(end,4:6)).*vNorm*1000)
fprintf('Total delta-V: %1.3f m/s\n', trapz(rowNorm(us_mpss)))


toc
% ========================================================================
% ========================================================================
%%% Functions
% ========================================================================
% ========================================================================
function [dK] = Int_K(t,K,B,R,Q,Xref,tref,u)

%%% Interpolating reference state for current values
Xref_now = interp1(tref,Xref,t);

%%% Getting A-Matrix
[A] = get_Amat_CR3BP(u, Xref_now(1), Xref_now(2), Xref_now(3));

K = reshape(K,6,6);

%%% Preallocate state output
dK = -K*A + K*B*inv(R)*B'*K - Q - A'*K;
dK = reshape(dK,36,1);
end


function [dX] = Int_State(t,X,B,R,K,Xref,tref,u,dt)

if t < (tref(end)-dt)
    % -------------------------------
    %%% Interpolating reference and calculating A matrix
    % -------------------------------
    %%% Interpolating reference state for current values
    Xref_now = interp1(tref,Xref,t);
    Xref_now = Xref_now';
    % 
    % %%% Getting A-Matrix
    % [A] = get_Amat_CR3BP(u, Xref_now(1), Xref_now(2), Xref_now(3));

    % -------------------------------
    %%% Creating controls
    % -------------------------------
    K_now = interp1(tref,K,t);
    K = reshape(K_now,6,6);


    deltaX = X - Xref_now;
    u_t = -(inv(R)*B'*K)*deltaX;
%     u_t = [0;0;0];
    % -------------------------------
    %%% Integrating state
    % -------------------------------
    %%% Preallocate state output
    dX = zeros(6,1);

    %%% Distances to primary (1) and secondary (2) bodies
    r1 = sqrt((X(1)+u)^2 + X(2)^2 + X(3)^2);
    r2 = sqrt((X(1)+u-1)^2 + X(2)^2 + X(3)^2);

    %%% Equations of Motion
    ddx = 2*X(5) + X(1) - (1-u)*(X(1)+u)/(r1^3) - u*(X(1)+u-1)/(r2^3) + u_t(1);
    ddy = -2*X(4) + X(2) -((1-u)/(r1^3) + u/(r2^3))*X(2) + u_t(2);
    ddz = -((1-u)/(r1^3) + u/(r2^3))*X(3) + u_t(3);

    %%% Output the derivative of the state
    dX(1:3) = [X(4); X(5); X(6)]; % km/s
    dX(4:6) = [ddx; ddy; ddz]; % km/s^2
else
    % -------------------------------
    %%% Integrating state
    % -------------------------------
    %%% Preallocate state output
    dX = zeros(6,1);

    %%% Distances to primary (1) and secondary (2) bodies
    r1 = sqrt((X(1)+u)^2 + X(2)^2 + X(3)^2);
    r2 = sqrt((X(1)+u-1)^2 + X(2)^2 + X(3)^2);

    %%% Equations of Motion
    ddx = 2*X(5) + X(1) - (1-u)*(X(1)+u)/(r1^3) - u*(X(1)+u-1)/(r2^3);
    ddy = -2*X(4) + X(2) -((1-u)/(r1^3) + u/(r2^3))*X(2);
    ddz = -((1-u)/(r1^3) + u/(r2^3))*X(3);

    %%% Output the derivative of the state
    dX(1:3) = [X(4); X(5); X(6)]; % km/s
    dX(4:6) = [ddx; ddy; ddz]; % km/s^2
end
end













