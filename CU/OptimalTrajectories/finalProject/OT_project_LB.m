clear
clc
close all
mbinPath = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/mbin';
addpath(genpath(mbinPath))
tic

% ========================================================================
%%% Run Switches
% ========================================================================
run_symbollicWork = 0;
% ========================================================================
%%% Importing Data
% ========================================================================
%%% General data on solar system bodies
bodies = getBodyData(mbinPath);

%%% Color options/schemes
colors = get_colors();

% ========================================================================
%%% Symbollic work
% ========================================================================
if run_symbollicWork == 1
syms u x y z dx dy dz real

r1 = sqrt((-u-x)^2+y^2+z^2);
r2 = sqrt((x-1+u)^2+y^2+z^2);

%%% Defining Equations of Motion, State Variables, and A Matrix
ddx = 2*dy + x - (1-u)*(x+u)/(r1^3) - u*(x+u-1)/(r2^3);
ddy = -2*dx + y -((1-u)/(r1^3) + u/(r2^3))*y;
ddz = -((1-u)/(r1^3) + u/(r2^3))*z;

EOM = [dx; dy; dz; ddx; ddy; ddz];
state = [x; y; z; dx; dy; dz];
A = jacobian(EOM, state);

989
return

end % run_symbollic work
% ========================================================================
%%% Setup
% ========================================================================

% -------------------------------
% System info
% -------------------------------
%%% bodies
primary = bodies.jupiter;
secondary = bodies.europa;

%%% Normalizing constants
rNorm = secondary.a;         % n <-> km
tNorm = 1/secondary.meanMot; % n <-> sec
vNorm = rNorm / tNorm;       % n <-> km/sec

%%% Lagrange Points
L123 = EquilibriumPoints(secondary.MR,1:3); % [3x3] of L1, L2, L3 normalized BCR coordinates
% -------------------------------
% Test integration
% -------------------------------
%%% Setting time vector
t_0 = 0;
t_f = (10*3600)/tNorm;
% t_f = 0.025944599748468;
n_dt = 5000;
time0_n = linspace(1/tNorm,t_f,n_dt);

%%% extra parameters
prms = struct();
prms.u = secondary.MR;
prms.R2_n = secondary.R_n;

%%% Choosing ode tolerance
tol = 1e-13;


% X0_n = [1.0204617015266166,-0.0030844497599572,-0.0085995858200228,-0.0028498125830546,-0.0189072710942121,0.0000000000000000]';
X0_n = [L123(2,1)-10/rNorm, 0, 0, 0, 0, 0]';
% Xref_n = [L123(2,:), 0, 0, 0]';
Xtgt_n = [L123(2,1)-(25/rNorm), 0, 0, 0, 0, 0]';
% Xref_n = [L123(2,1)-(0.005/rNorm), (0.02/rNorm), 0, 0, 0, 0]';
Xref_n = X0_n;
dX_f = Xtgt_n - X0_n;
dX_f = dX_f + [0.148351025408555,0.000000006683836,0,-0.001875221345827,-0.000091612877502,0]'*1e-4 - [-0.218177103050721,0.000023543081683,0,0.287277126022772,-0.137782950149768,0]'*1e-7 - [-0.027911895017496,-0.002010232497417,0,0.189472686734755,-0.130509616057475,0]'*1e-9;

% 
% [time_n, X_BCR_n, time_eventImpact, X_eventImpact, index_eventImpact] = ode113(@Int_CR3Bn,...
%                 time0_n, X0_n, options_ImpactEscape, prms);
%             
% figure; hold all
% plot3(X_BCR_n(:,1),X_BCR_n(:,2),X_BCR_n(:,3),'b','linewidth',2)
% plotBodyTexture3(secondary.R_n, [1-secondary.MR, 0, 0], secondary.img)
% axis equal
% PlotBoi3('$x$','$y$','$z$',16,'LaTex')


% -------------------------------
% A-matrix
% -------------------------------
[A] = get_Amat_CR3BP(secondary.MR, 1.2,.3,.13);


% dX = CR3Bn_EOM(X0_n,secondary.MR)



% -------------------------------
% 
% -------------------------------

% dXf_r_tol = 1e-7;
% dXf_v_tol = 1e-3;
dXf_r_tol_meters = 0.1;
dXf_v_tol_mps = 1e-6;
dXf_r_tol = (dXf_r_tol_meters/1000)/rNorm;
dXf_v_tol = (dXf_v_tol_mps/1000)/vNorm;

umax_mpss = 20; % m/s^2
umax = umax_mpss * (1/1000) * ((tNorm^2)/rNorm);
% umax = 100000 * ((tNorm^2)/rNorm);

% options_OTconvergence = odeset('Events',@event_OTproject_convergence,'RelTol',tol,'AbsTol',tol);
% [time_n, X_BCR_n, time_eventImpact, X_eventImpact, index_eventImpact] = ode113(@Int_CR3Bn_OTproject,...
%                 time0_n, X0_n, options_OTconvergence,Xref_n, t_0,t_f, dX_f, dXf_r_tol,dXf_v_tol,umax,Xf_n, prms);
% time_eventImpact

options = odeset('RelTol',tol,'AbsTol',tol);
[time_n, X_BCR_n] = ode113(@Int_CR3Bn_OTproject,...
                time0_n(1:end-1), X0_n, options,Xref_n, t_0,t_f, dX_f, dXf_r_tol,dXf_v_tol,umax,Xtgt_n, prms);

%%% Recreate controls
us = zeros(length(time_n),3);
for kk = 1:length(time_n)
    [u_i] = getControls(time_n(kk),X_BCR_n(kk,:)',Xref_n,t_0,t_f,dX_f,umax,prms);
    us(kk,:) = u_i';
end

fprintf('Final position error: %1.2f m\n',norm(X_BCR_n(end,1:3)' - Xtgt_n(1:3))*rNorm*1000)
fprintf('Final velcotiy error: %1.2f mm/s\n',norm(X_BCR_n(end,4:6)' - Xtgt_n(4:6))*vNorm*(1e6))

% figure('position',[-648 549 560 420]); hold all
figure; hold all
plot3(X_BCR_n(:,1),X_BCR_n(:,2),X_BCR_n(:,3),'b','linewidth',1.5)
plot3(X0_n(1),X0_n(2),X0_n(3),'ro','markersize',10)
plot3(Xtgt_n(1),Xtgt_n(2),Xtgt_n(3),'rx','markersize',10)
plot3(L123(2,1),0,0,'^','markersize',10,'markeredgecolor',colors.std.black,'markerfacecolor',colors.std.grn)

% figure('position',[-1573 132 886 847])
figure
subplot(3,2,1)
plot(time_n.*tNorm, (X_BCR_n(:,1) - Xtgt_n(1)).*rNorm*1000,'r','linewidth',1.5)
PlotBoi2('','$x_{error}, m$',20,'LaTex')
subplot(3,2,3)
plot(time_n.*tNorm, (X_BCR_n(:,2) - Xtgt_n(2)).*rNorm*1000,'r','linewidth',1.5)
PlotBoi2('','$y_{error}, m$',20,'LaTex')
subplot(3,2,5)
plot(time_n.*tNorm, (X_BCR_n(:,3) - Xtgt_n(3)).*rNorm*1000,'r','linewidth',1.5)
PlotBoi2('Time, $sec$','$z_{error}, m$',20,'LaTex')

subplot(3,2,2)
plot(time_n.*tNorm, (X_BCR_n(:,4) - Xtgt_n(4)).*vNorm*1000,'r','linewidth',1.5)
PlotBoi2('','$\dot{x}_{error}, m/s$',20,'LaTex')
subplot(3,2,4)
plot(time_n.*tNorm, (X_BCR_n(:,5) - Xtgt_n(5)).*vNorm*1000,'r','linewidth',1.5)
PlotBoi2('','$\dot{y}_{error}, m/s$',20,'LaTex')
subplot(3,2,6)
plot(time_n.*tNorm, (X_BCR_n(:,6) - Xtgt_n(6)).*vNorm*1000,'r','linewidth',1.5)
PlotBoi2('Time, $sec$','$\dot{z}_{error}$, m/s',20,'LaTex')


% figure('position',[-647 42 560 420])
figure
subplot(3,1,1)
plot(time_n.*tNorm,us(:,1).*(1000*rNorm/(tNorm^2)),'b','linewidth',1.5)
PlotBoi2('','$u_x, m/s^2$',20,'LaTex')
subplot(3,1,2)
plot(time_n.*tNorm,us(:,2).*(1000*rNorm/(tNorm^2)),'b','linewidth',1.5)
PlotBoi2('','$u_y, m/s^2$',20,'LaTex')
subplot(3,1,3)
plot(time_n.*tNorm,us(:,3).*(1000*rNorm/(tNorm^2)),'b','linewidth',1.5)
PlotBoi2('Time, $sec$','$u_z, m/s^2$',20,'LaTex')



















toc
% ========================================================================
% ========================================================================
%%% Functions
% ========================================================================
% ========================================================================
%%% Function for returning Xdot in CR3BP
function [dX] = CR3Bn_EOM(X,u)

%%% Preallocate state output
dX = zeros(6,1);

%%% Unpack the barycentric state vector
x = X(1); y = X(2); z = X(3); % Position
dx = X(4); dy = X(5); dz = X(6); % Velocity

%%% Distances to primary (1) and secondary (2) bodies
r1 = sqrt((x+u)^2 + y^2 + z^2);
r2 = sqrt((x+u-1)^2 + y^2 + z^2);

%%% Equations of Motion
ddx = 2*dy + x - (1-u)*(x+u)/(r1^3) - u*(x+u-1)/(r2^3);
ddy = -2*dx + y -((1-u)/(r1^3) + u/(r2^3))*y;
ddz = -((1-u)/(r1^3) + u/(r2^3))*z;

%%% Output the derivative of the state
dX(1:3) = [dx; dy; dz];
dX(4:6) = [ddx; ddy; ddz];

end

%%% Function for integrating X in CR3BP
function [dX] = Int_CR3Bn_OTproject(t,X,Xref,t0,tf,dX_f,dXf_r_tol,dXf_v_tol,umax,Xtgt_n,prms)
%%% For numerical integration in the normalized CR3BP
%%% Inputs:
%          t - normalized time vector
%          X - initial state [6x1]
%          prms - (u)

%%% Preallocate state output
dX = zeros(6,1);

%%% Distances to primary (1) and secondary (2) bodies
r1 = sqrt((X(1)+prms.u)^2 + X(2)^2 + X(3)^2);
r2 = sqrt((X(1)+prms.u-1)^2 + X(2)^2 + X(3)^2);


% dX_f = Xtgt_n - Xref;
[u_i] = getControls(t,X,Xref,t0,tf,dX_f,umax,prms);

% norm(u_i)*1000*2.814322571059422e-04





%%% Equations of Motion
ddx = 2*X(5) + X(1) - (1-prms.u)*(X(1)+prms.u)/(r1^3) - prms.u*(X(1)+prms.u-1)/(r2^3) + u_i(1);
ddy = -2*X(4) + X(2) -((1-prms.u)/(r1^3) + prms.u/(r2^3))*X(2) + u_i(2);
ddz = -((1-prms.u)/(r1^3) + prms.u/(r2^3))*X(3) + u_i(3);

%%% Output the derivative of the state
dX(1:3) = [X(4); X(5); X(6)]; % km/s
dX(4:6) = [ddx; ddy; ddz]; % km/s^2
end


function [u_i] = getControls(t,X,Xref,t0,tf,dX_f,umax,prms)

[A] = get_Amat_CR3BP(prms.u, X(1), X(2), X(3));

dX_i = X - Xref;

%%% building M
M_UL = A;
M_UR = blkdiag(zeros(3,3),-eye(3));
M_BL = zeros(6,6);
M_BR = -A';
M = [M_UL, M_UR;...
     M_BL, M_BR];
 
f = CR3Bn_EOM(Xref,prms.u);

% warning('His f is just the 3 accelerations')
N = [f; zeros(6,1)];
% warning('He might have implied that some of these^ zeros are vectors')

Q = inv(M)*N;
Q_X = Q(1:6,1);
% warning('try phi_dot = M phi')
% warning('try pre-computing phi')

phi_t = expm(M.*(tf-t)); %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
phi_XX = phi_t(1:6,1:6);
phi_XP = phi_t(1:6,7:12);
 
P_i = inv(phi_XP) * (dX_f - phi_XX*dX_i + Q_X); 

Pv_i = P_i(4:6);
u_i = -Pv_i;

if norm(u_i) > umax
    uHat = u_i ./ norm(u_i);
    u_i = uHat .* umax;
end

end

% 
% function [value, isterminal, direction] = event_OTproject_convergence(t,X,Xref,t0,tf,dX_f,dXf_r_tol,dXf_v_tol,umax,Xtgt_n,prms)
% %%% Designed for standard normalized CR3BP
% %%% Event function watching for when "value" = 0
% %%% Inputs:
% %          t - normalized time vector
% %          X - initial state [6x1]
% %          prms - (u, R2_n, L1x, L2x)
% 
% % value = [norm((norm(X(1:3) - Xf_n(1:3))-dXf_r_tol)) + norm((norm(X(4:6) - Xf_n(4:6))-dXf_v_tol))]; 
% value = [(norm(X(1:3) - Xf_n(1:3))-dXf_r_tol) + (norm(X(4:6) - Xf_n(4:6))-dXf_v_tol)]; 
% % value = norm(X - Xf_n);
% isterminal = [1]; % stops the integration
% direction = [0]; % negative direction only
% 
% % dX_i = X - Xref;
% % 
% % value = [norm(dX_i(1:3) - dX_f(1:3)) - dXf_r_tol, norm(dX_i(4:6) - dX_f(4:6)) - dXf_v_tol]; 
% % isterminal = [1, 1]; % stops the integration
% % direction = [0, 0]; % negative direction only
% end

% 
% 
% function [ stm_dot_2 ] = OTproject_STMInt(stm_in_1,M)
% %     %%% Unpack state
% %     x = stm_in(1);
% %     y = stm_in(2);
% %     z = stm_in(3);
% %     dx = stm_in(4);
% %     dy = stm_in(5);
% %     dz = stm_in(6);
% 
%     %%% Reshape (n^2,1) stm to (n,n)
%     stm_in_2 = reshape(stm_in_1,12,12);
% 
% %     %%% Build A matrix and evaluate at current state
% %     A = zeros(6,6);
% %     A(1:3,4:6) = eye(3,3);
% %     % Asym(4,1)
% %     A(4,1) = -(u/(x^2 + y^2 + z^2)^(3/2) - (3*u*x^2)/(x^2 + y^2 + z^2)^(5/2) - (3*J2*RE^2*u)/(2*(x^2 + y^2 + z^2)^(5/2)) - (15*J3*RE^3*u*z)/(2*(x^2 + y^2 + z^2)^(7/2)) + (15*J2*RE^2*u*x^2)/(2*(x^2 + y^2 + z^2)^(7/2)) + (15*J2*RE^2*u*z^2)/(2*(x^2 + y^2 + z^2)^(7/2)) + (35*J3*RE^3*u*z^3)/(2*(x^2 + y^2 + z^2)^(9/2)) + (105*J3*RE^3*u*x^2*z)/(2*(x^2 + y^2 + z^2)^(9/2)) - (105*J2*RE^2*u*x^2*z^2)/(2*(x^2 + y^2 + z^2)^(9/2)) - (315*J3*RE^3*u*x^2*z^3)/(2*(x^2 + y^2 + z^2)^(11/2)));
% %     % Asym(4,2)
% %     A(4,2) = -((15*J2*RE^2*u*x*y)/(2*(x^2 + y^2 + z^2)^(7/2)) - (3*u*x*y)/(x^2 + y^2 + z^2)^(5/2) + (105*J3*RE^3*u*x*y*z)/(2*(x^2 + y^2 + z^2)^(9/2)) - (105*J2*RE^2*u*x*y*z^2)/(2*(x^2 + y^2 + z^2)^(9/2)) - (315*J3*RE^3*u*x*y*z^3)/(2*(x^2 + y^2 + z^2)^(11/2)));
% %     % Asym(4,3)
% %     A(4,3) = -((45*J2*RE^2*u*x*z)/(2*(x^2 + y^2 + z^2)^(7/2)) - (15*J3*RE^3*u*x)/(2*(x^2 + y^2 + z^2)^(7/2)) - (3*u*x*z)/(x^2 + y^2 + z^2)^(5/2) - (105*J2*RE^2*u*x*z^3)/(2*(x^2 + y^2 + z^2)^(9/2)) + (105*J3*RE^3*u*x*z^2)/(x^2 + y^2 + z^2)^(9/2) - (315*J3*RE^3*u*x*z^4)/(2*(x^2 + y^2 + z^2)^(11/2)));
% %     % Asym(5,1)
% %     A(5,1) = -((15*J2*RE^2*u*x*y)/(2*(x^2 + y^2 + z^2)^(7/2)) - (3*u*x*y)/(x^2 + y^2 + z^2)^(5/2) + (105*J3*RE^3*u*x*y*z)/(2*(x^2 + y^2 + z^2)^(9/2)) - (105*J2*RE^2*u*x*y*z^2)/(2*(x^2 + y^2 + z^2)^(9/2)) - (315*J3*RE^3*u*x*y*z^3)/(2*(x^2 + y^2 + z^2)^(11/2)));
% %     % Asym(5,2)
% %     A(5,2) = -(u/(x^2 + y^2 + z^2)^(3/2) - (3*u*y^2)/(x^2 + y^2 + z^2)^(5/2) - (3*J2*RE^2*u)/(2*(x^2 + y^2 + z^2)^(5/2)) - (15*J3*RE^3*u*z)/(2*(x^2 + y^2 + z^2)^(7/2)) + (15*J2*RE^2*u*y^2)/(2*(x^2 + y^2 + z^2)^(7/2)) + (15*J2*RE^2*u*z^2)/(2*(x^2 + y^2 + z^2)^(7/2)) + (35*J3*RE^3*u*z^3)/(2*(x^2 + y^2 + z^2)^(9/2)) + (105*J3*RE^3*u*y^2*z)/(2*(x^2 + y^2 + z^2)^(9/2)) - (105*J2*RE^2*u*y^2*z^2)/(2*(x^2 + y^2 + z^2)^(9/2)) - (315*J3*RE^3*u*y^2*z^3)/(2*(x^2 + y^2 + z^2)^(11/2)));
% %     % Asym(5,3)
% %     A(5,3) = -((45*J2*RE^2*u*y*z)/(2*(x^2 + y^2 + z^2)^(7/2)) - (15*J3*RE^3*u*y)/(2*(x^2 + y^2 + z^2)^(7/2)) - (3*u*y*z)/(x^2 + y^2 + z^2)^(5/2) - (105*J2*RE^2*u*y*z^3)/(2*(x^2 + y^2 + z^2)^(9/2)) + (105*J3*RE^3*u*y*z^2)/(x^2 + y^2 + z^2)^(9/2) - (315*J3*RE^3*u*y*z^4)/(2*(x^2 + y^2 + z^2)^(11/2)));
% %     % Asym(6,1)
% %     A(6,1) = -((45*J2*RE^2*u*x*z)/(2*(x^2 + y^2 + z^2)^(7/2)) - (15*J3*RE^3*u*x)/(2*(x^2 + y^2 + z^2)^(7/2)) - (3*u*x*z)/(x^2 + y^2 + z^2)^(5/2) - (105*J2*RE^2*u*x*z^3)/(2*(x^2 + y^2 + z^2)^(9/2)) + (105*J3*RE^3*u*x*z^2)/(x^2 + y^2 + z^2)^(9/2) - (315*J3*RE^3*u*x*z^4)/(2*(x^2 + y^2 + z^2)^(11/2)));
% %     % Asym(6,2)
% %     A(6,2) = -((45*J2*RE^2*u*y*z)/(2*(x^2 + y^2 + z^2)^(7/2)) - (15*J3*RE^3*u*y)/(2*(x^2 + y^2 + z^2)^(7/2)) - (3*u*y*z)/(x^2 + y^2 + z^2)^(5/2) - (105*J2*RE^2*u*y*z^3)/(2*(x^2 + y^2 + z^2)^(9/2)) + (105*J3*RE^3*u*y*z^2)/(x^2 + y^2 + z^2)^(9/2) - (315*J3*RE^3*u*y*z^4)/(2*(x^2 + y^2 + z^2)^(11/2)));
% %     % Asym(6,3)
% %     A(6,3) = -(u/(x^2 + y^2 + z^2)^(3/2) - (3*u*z^2)/(x^2 + y^2 + z^2)^(5/2) - (9*J2*RE^2*u)/(2*(x^2 + y^2 + z^2)^(5/2)) - (75*J3*RE^3*u*z)/(2*(x^2 + y^2 + z^2)^(7/2)) + (45*J2*RE^2*u*z^2)/(x^2 + y^2 + z^2)^(7/2) - (105*J2*RE^2*u*z^4)/(2*(x^2 + y^2 + z^2)^(9/2)) + (175*J3*RE^3*u*z^3)/(x^2 + y^2 + z^2)^(9/2) - (315*J3*RE^3*u*z^5)/(2*(x^2 + y^2 + z^2)^(11/2)));
% 
%     %%% Calculate new STM
%     stm_dot_1 = M*stm_in_2;
% 
% %     %%% Creating X-dot
% %     % Spacecraft velocities
% %     dX(1:3) = [dx; dy; dz];
% %     
% %     
% %     %%% Using u, J2, and J3 terms as dynamics
% %     % -EQM(4)
% %     dX(4) = (3*J2*RE^2*u*x)/(2*(x^2 + y^2 + z^2)^(5/2)) - (u*x)/(x^2 + y^2 + z^2)^(3/2) + (15*J3*RE^3*u*x*z)/(2*(x^2 + y^2 + z^2)^(7/2)) - (15*J2*RE^2*u*x*z^2)/(2*(x^2 + y^2 + z^2)^(7/2)) - (35*J3*RE^3*u*x*z^3)/(2*(x^2 + y^2 + z^2)^(9/2)) - aDrag(1);
% %     % -EQM(5)
% %     dX(5) = (3*J2*RE^2*u*y)/(2*(x^2 + y^2 + z^2)^(5/2)) - (u*y)/(x^2 + y^2 + z^2)^(3/2) + (15*J3*RE^3*u*y*z)/(2*(x^2 + y^2 + z^2)^(7/2)) - (15*J2*RE^2*u*y*z^2)/(2*(x^2 + y^2 + z^2)^(7/2)) - (35*J3*RE^3*u*y*z^3)/(2*(x^2 + y^2 + z^2)^(9/2)) - aDrag(2);
% %     % -EQM(6)
% %     dX(6) = (9*J2*RE^2*u*z)/(2*(x^2 + y^2 + z^2)^(5/2)) - (3*J3*RE^3*u)/(2*(x^2 + y^2 + z^2)^(5/2)) - (u*z)/(x^2 + y^2 + z^2)^(3/2) - (15*J2*RE^2*u*z^3)/(2*(x^2 + y^2 + z^2)^(7/2)) + (15*J3*RE^3*u*z^2)/(x^2 + y^2 + z^2)^(7/2) - (35*J3*RE^3*u*z^4)/(2*(x^2 + y^2 + z^2)^(9/2)) - aDrag(3);
% 
%     % Filling in reshaped (6^2,1) STM to state
% %     dX(7:end) = reshape(stm_dot,36,1);
%     stm_dot_2 = reshape(stm_dot_1,12*12,1);
% end
% 
% 







