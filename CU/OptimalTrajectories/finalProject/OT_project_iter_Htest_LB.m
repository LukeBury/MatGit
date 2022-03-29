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

plot_secondary = 1;

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
%%%%%%%%%%%%%%%%%%%%%%%%%% Europa 1
%%% Setting time vector


%%% ------------------
% X0_n - Starting Location (fixed)
% Xtgt_n - Final, target location
% Xref_n - State to be linearized about
% dX_f - State vector changed needed to hit target
%%% ------------------

X0_n = [L123(2,1)-10/rNorm, 0, 100/rNorm, 0, 0, 0]';
% Xref_n = [L123(2,:), 0, 0, 0]';
% Xtgt_n = [L123(2,1)-(11500/rNorm), 0, 0, 0, 0, 0]';
Xtgt_n = [1-secondary.MR+secondary.R_n+1000/rNorm, 0, 0, 0, 0, 0]';
Xlin_n = X0_n;


%%%%%%%%%%%%%%%%%%%%%%%%%%% Europa 2
% %%% Setting time vector
% t_0 = 0;
% t_f = (1*3600)/tNorm;
% % t_f = 0.025944599748468;
% n_dt = 5000;
% time0_n = linspace(1/tNorm,t_f,n_dt);
% 
% X0_n = [L123(2,1)-(10000/rNorm), 0, 0, 0, 0, 0]';
% % Xref_n = [L123(2,:), 0, 0, 0]';
% [rBCR] = latlon2surfCR3BP('secondary', 0,-90, secondary.MR, secondary.R_n);
% Xtgt_n = [rBCR', 0, 0, 0]';
% Xref_n = X0_n;
% dX_f = Xtgt_n - X0_n;

%%%%%%%%%%%%%%%%%%%%%%%%%%% Enceladus 1
% %%% Setting time vector
% t_0 = 0;
% t_f = (20*3600)/tNorm;
% % t_f = 0.025944599748468;
% n_dt = 5000;
% time0_n = linspace(1/tNorm,t_f,n_dt);
% 
% X0_n = [L123(2,1)-(300/rNorm), 0, 0, 0, 0, 0]';
% % Xref_n = [L123(2,:), 0, 0, 0]';
% [rBCR] = latlon2surfCR3BP('secondary', 0,-90, secondary.MR, secondary.R_n);
% rBCR = rBCR + [0; 100/rNorm; 0];
% Xtgt_n = [rBCR', 0, 0, 0]';
% Xref_n = X0_n;
% dX_f = Xtgt_n - X0_n;

% umax_mpss = 10; % m/s^2

umax_mpss = 100000000; % m/s^2


%%% extra parameters
prms = struct();
prms.u = secondary.MR;
prms.R2_n = secondary.R_n;

%%% Choosing ode tolerance
tol = 1e-13;

%%% Iteration stuff
n_iter = 6;
iterColors = colorScale([colors.std.ltred; colors.std.blue],n_iter);

%%% Preallocating
plotHandles{n_iter} = [];
legendNames{n_iter} = [];



% -------------------------------
% 
% -------------------------------

% % dXf_r_tol = 1e-7;
% % dXf_v_tol = 1e-3;
% dXf_r_tol_meters = 0.1;
% dXf_v_tol_mps = 1e-6;
% dXf_r_tol = (dXf_r_tol_meters/1000)/rNorm;
% dXf_v_tol = (dXf_v_tol_mps/1000)/vNorm;

umax_n = umax_mpss * (1/1000) * ((tNorm^2)/rNorm);
% umax = 100000 * ((tNorm^2)/rNorm);

% options_OTconvergence = odeset('Events',@event_OTproject_convergence,'RelTol',tol,'AbsTol',tol);
% [time_n, X_BCR_n, time_eventImpact, X_eventImpact, index_eventImpact] = ode113(@Int_CR3Bn_OTproject,...
%                 time0_n, X0_n, options_OTconvergence,Xref_n, t_0,t_f, dX_f, dXf_r_tol,dXf_v_tol,umax,Xf_n, prms);
% time_eventImpact

n_t_iter = 10;

t_0 = 0;
t_fs = linspace(0.5, 5, n_t_iter).*(3600/tNorm)./10;
% [5 6 7 8 9 10 15 20 35 50].*(3600/tNorm)./10;
% t_f = 0.025944599748468;
n_dt = 5000;

Hams{n_t_iter}.Hfs = zeros(length(n_dt),1);
for t_iter = 1:n_t_iter
dX_f = Xtgt_n - X0_n;

for iter = 1:n_iter
t_f = t_fs(t_iter);
time0_n = linspace(1/tNorm,t_f,n_dt);
    
options = odeset('RelTol',tol,'AbsTol',tol);
[time_n, X_BCR_n] = ode113(@Int_CR3Bn_OTproject,...
                time0_n(1:end-1), X0_n, options,Xlin_n, t_0,t_f, dX_f,umax_n, prms);

%%% Parameters for H calculation
drf = dX_f(1:3,1);
dvf = dX_f(4:6,1);
Af = get_Amat_CR3BP(prms.u, X_BCR_n(end,1),X_BCR_n(end,2),X_BCR_n(end,3));
Af_BL = Af(4:6,1:3);
Af_BR = Af(4:6,4:6);
Cf = CR3Bn_EOM(Xlin_n,prms.u);
Cr = Cf(1:3,1);
Cv = Cf(4:6,1);
    
%%% Recreate controls
us = zeros(length(time_n),3);

for kk = 1:length(time_n)
    [u_i, Pf, P_i, dX_i] = getControls(time_n(kk),X_BCR_n(kk,:)',Xlin_n,t_0,t_f,dX_f,umax_n,prms);
    us(kk,:) = u_i';
    
    %%% Computing Hf and Hi
    Prf = Pf(1:3,1);
    Pvf = Pf(4:6,1);
    Pri = P_i(1:3,1);
    Pvi = P_i(4:6,1);
    dri = dX_i(1:3,1);
    dvi = dX_i(4:6,1);
    Ai = get_Amat_CR3BP(prms.u, X_BCR_n(kk,1),X_BCR_n(kk,2),X_BCR_n(kk,3));
    Ai_BL = Ai(4:6,1:3);
    Ai_BR = Ai(4:6,4:6);


%     Hams{t_iter}.Hs1(kk) = 0.5*(u_i'*u_i) + Pri'*dvi + Pvi'*Ai_BL*dri + Pvi'*Ai_BR*dvi + Pvi'*u_i + Pri'*Cr + Pvf'*Cv;
    Hams{t_iter}.Hs(kk) = Pri'*dvi + Pvi'*Ai_BL*dri + Pvi'*Ai_BR*dvi - 0.5*(Pvi'*Pvi) + Pri'*Cr + Pvi'*Cv;
    
    Hams{t_iter}.Hfs(kk) = Prf'*dvf + Pvf'*Af_BL*drf + Pvf'*Af_BR*dvf - 0.5*(Pvf'*Pvf) + Prf'*Cr + Pvf'*Cv;
end
% figure(t_iter); hold all
% plot(Hams{t_iter}.Hs(:))
fprintf('t_iter: %1d, iter: %1d\n',t_iter,iter)
if iter == n_iter
    t_f
    Hams{t_iter}.Hfs(1)
end



% figure(1)
% subplot(1,2,1); hold all
% plotHandles{iter} = plot3(X_BCR_n(:,1),X_BCR_n(:,2),X_BCR_n(:,3),'linewidth',1.5,'color',iterColors(iter,:));
% plot3(X0_n(1),X0_n(2),X0_n(3),'bo','markersize',10)
% plot3(Xtgt_n(1),Xtgt_n(2),Xtgt_n(3),'bx','markersize',10)
% plot3(L123(2,1),0,0,'^','markersize',10,'markeredgecolor',colors.std.black,'markerfacecolor',colors.std.grn)
% PlotBoi3('$x_n$','$y_n$','$z_n$',20,'LaTex')
% axis equal
% legendNames{iter} = sprintf('Iteration %1d',iter);
% 
% subplot(1,2,2); hold all
% plotHandles{iter} = plot3(X_BCR_n(:,1),X_BCR_n(:,2),X_BCR_n(:,3),'linewidth',1.5,'color',iterColors(iter,:));
% plot3(X0_n(1),X0_n(2),X0_n(3),'bo','markersize',10)
% plot3(Xtgt_n(1),Xtgt_n(2),Xtgt_n(3),'bx','markersize',10)
% % plot3(L123(2,1),0,0,'^','markersize',10,'markeredgecolor',colors.std.black,'markerfacecolor',colors.std.grn)
% PlotBoi3('$x_n$','$y_n$','$z_n$',20,'LaTex')
% axis equal
% 
% % figure('position',[-1573 132 886 847])
% figure(2)
% subplot(3,2,1); hold all
% plot(time_n.*tNorm, (X_BCR_n(:,1) - Xtgt_n(1)).*rNorm*1000,'linewidth',1.5,'color',iterColors(iter,:))
% PlotBoi2('','$x_{error}, m$',20,'LaTex')
% subplot(3,2,3); hold all
% plot(time_n.*tNorm, (X_BCR_n(:,2) - Xtgt_n(2)).*rNorm*1000,'linewidth',1.5,'color',iterColors(iter,:))
% PlotBoi2('','$y_{error}, m$',20,'LaTex')
% subplot(3,2,5); hold all
% plot(time_n.*tNorm, (X_BCR_n(:,3) - Xtgt_n(3)).*rNorm*1000,'linewidth',1.5,'color',iterColors(iter,:))
% PlotBoi2('Time, $sec$','$z_{error}, m$',20,'LaTex')
% 
% subplot(3,2,2); hold all
% plot(time_n.*tNorm, (X_BCR_n(:,4) - Xtgt_n(4)).*vNorm*1000,'linewidth',1.5,'color',iterColors(iter,:))
% PlotBoi2('','$\dot{x}_{error}, m/s$',20,'LaTex')
% subplot(3,2,4); hold all
% plot(time_n.*tNorm, (X_BCR_n(:,5) - Xtgt_n(5)).*vNorm*1000,'linewidth',1.5,'color',iterColors(iter,:))
% PlotBoi2('','$\dot{y}_{error}, m/s$',20,'LaTex')
% subplot(3,2,6); hold all
% plot(time_n.*tNorm, (X_BCR_n(:,6) - Xtgt_n(6)).*vNorm*1000,'linewidth',1.5,'color',iterColors(iter,:))
% PlotBoi2('Time, $sec$','$\dot{z}_{error}$, m/s',20,'LaTex')

% 
% % figure('position',[-647 42 560 420])
% figure(3)
% subplot(3,1,1); hold all
% plot(time_n.*tNorm,us(:,1).*(1000*rNorm/(tNorm^2)),'linewidth',1.5,'color',iterColors(iter,:))
% PlotBoi2('','$u_x, m/s^2$',20,'LaTex')
% subplot(3,1,2); hold all
% plot(time_n.*tNorm,us(:,2).*(1000*rNorm/(tNorm^2)),'linewidth',1.5,'color',iterColors(iter,:))
% PlotBoi2('','$u_y, m/s^2$',20,'LaTex')
% subplot(3,1,3); hold all
% plot(time_n.*tNorm,us(:,3).*(1000*rNorm/(tNorm^2)),'linewidth',1.5,'color',iterColors(iter,:))
% PlotBoi2('Time, $sec$','$u_z, m/s^2$',20,'LaTex')



dX_f = dX_f - (X_BCR_n(end,:)' - Xtgt_n);

end % iter
end % t_iter = 1:n_t_iter

if plot_secondary == 1
    figure(1)
    subplot(1,2,2)
    plotBodyTexture3(secondary.R_n,[1-secondary.MR,0,0],secondary.img)
    view(0,90)
end

%%% Creating legend
legend([plotHandles{:}],legendNames)










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
function [dX] = Int_CR3Bn_OTproject(t,X,Xlin_n,t0,tf,dX_f,umax_n,prms)
%%% Inputs:
%          t - normalized time vector
%          X - initial state [6x1]
%          Xlin_n - State being linearized about [6x1]
%          t0 - t0 of integration time
%          tf - tf of integration time
%          dX_f - State vector changed needed to hit target
%          umax - maximum magnitude of control vector
%          prms - (u)
% ======================================================================
%%% Preallocate state output
dX = zeros(6,1);

%%% Distances to primary (1) and secondary (2) bodies
r1 = sqrt((X(1)+prms.u)^2 + X(2)^2 + X(3)^2);
r2 = sqrt((X(1)+prms.u-1)^2 + X(2)^2 + X(3)^2);

% dX_f = Xtgt_n - Xref;

%%% Calculate control vector
[u_i, ~, ~, ~] = getControls(t,X,Xlin_n,t0,tf,dX_f,umax_n,prms);

%%% Equations of Motion
ddx = 2*X(5) + X(1) - (1-prms.u)*(X(1)+prms.u)/(r1^3) - prms.u*(X(1)+prms.u-1)/(r2^3) + u_i(1);
ddy = -2*X(4) + X(2) -((1-prms.u)/(r1^3) + prms.u/(r2^3))*X(2) + u_i(2);
ddz = -((1-prms.u)/(r1^3) + prms.u/(r2^3))*X(3) + u_i(3);

%%% Output the derivative of the state
dX(1:3) = [X(4); X(5); X(6)]; % km/s
dX(4:6) = [ddx; ddy; ddz]; % km/s^2
end


function [u_i, Pf, P_i, dX_i] = getControls(t,X,Xlin_n,t0,tf,dX_f,umax_n,prms)

% [A] = get_Amat_CR3BP(prms.u, X(1), X(2), X(3));
[A] = get_Amat_CR3BP(prms.u, Xlin_n(1), Xlin_n(2), Xlin_n(3));

dX_i = X - Xlin_n;

%%% building M
M_UL = A;
M_UR = blkdiag(zeros(3,3),-eye(3));
M_BL = zeros(6,6);
M_BR = -A';
M = [M_UL, M_UR;...
     M_BL, M_BR];

% % % %%% Calculating Xdot of linearization state (Xlin_n)
% % % f = CR3Bn_EOM(Xlin_n,prms.u);
% % % 
% % % %%% Building N
% % % N = [f; zeros(6,1)];
% % % % N = [zeros(12,1)];
% % % 
% % % %%% Building Q
% % % Q = inv(M)*N;
% % % Q_X = Q(1:6,1);
% % % Q_P = Q(7:12,1);
% % % % warning('try phi_dot = M phi')
% % % % warning('try pre-computing phi')

phi_tf_t = expm(M.*(tf-t)); 
phi_XX = phi_tf_t(1:6,1:6);
phi_XP = phi_tf_t(1:6,7:12);
phi_PP = phi_tf_t(7:12,7:12);

% % % P_i = inv(phi_XP) * (dX_f - phi_XX*dX_i + Q_X); 
% % % 
% % % Pf = phi_PP*P_i - Q_P;
 
P_i = inv(phi_XP) * (dX_f - phi_XX*dX_i); 

Pf = phi_PP*P_i;

Pv_i = P_i(4:6);
u_i = -Pv_i;

if norm(u_i) > umax_n
    uHat = u_i ./ norm(u_i);
    u_i = uHat .* umax_n;
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
% function [u_i, Pf, P_i, dX_i] = GetControls(t,X,Xlin_n,tf,dX_f,prms)
% %%% Inputs
% % t      - current time
% % X      - current state [6x1]
% % Xlin_n - state that I'm linearizing about [6x1]
% % tf     - given final time
% % dX_f   - target dX
% % prms.u - CR3BP system mass ratio
% % ==============================================================
% %%% Current A-matrix (dF/dX evaluated at X)
% [A] = get_Amat_CR3BP(prms.u, X(1), X(2), X(3));
% 
% %%% Current dX
% dX_i = X - Xlin_n;
% 
% %%% building M
% M_UL = A;
% M_UR = blkdiag(zeros(3,3),-eye(3));
% M_BL = zeros(6,6);
% M_BR = -A';
% M = [M_UL, M_UR;...
%      M_BL, M_BR];
% 
% %%% State transition matrix (tf,t)
% phi_tf_t = expm(M.*(tf-t)); 
% phi_XX = phi_tf_t(1:6,1:6);
% phi_XP = phi_tf_t(1:6,7:12);
% phi_PP = phi_tf_t(7:12,7:12);
% 
% %%% Current co-state
% P_i = inv(phi_XP) * (dX_f - phi_XX*dX_i); 
% 
% %%% Projected final co-state
% Pf = phi_PP*P_i;
% 
% %%% Creating controls
% Pv_i = P_i(4:6);
% u_i = -Pv_i;
% 
% end
% 
% 






