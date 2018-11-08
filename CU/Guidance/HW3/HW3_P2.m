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
%%% P2-a
% ========================================================================
%%% Gravitational parameter
u = bodies.earth.u;

%%% States
r0 = [6578; 0; 0]; % km
rf = [-42164*cosd(5); 42164*sind(5); 0]; % km

%%% Circular orbit velocity
v0 = [0; sqrt(u/norm(r0)); 0]; % km/s

%%% Flight time
tf = 37864; % sec

%%% Lambert Solver
[v1, v2, exitflag] = lambertSolver(r0,rf,tf,0,0,u);
v1
dv = v1-v0
fprintf('Iterations = 5\n')
v2


% % -------------------------------------------------
% % Testing solution
% % -------------------------------------------------
% %%% Propagation time
% t0 = 0;
% dt = 1;
% tf = 37864;
% time = t0:dt:tf;
% 
% %%% IC
% X0 = [r0; v1];
% 
% %%% Choosing ode tolerance
% tol = 1e-10;
% 
% %%% Setting integrator options
% options = odeset('RelTol',tol,'AbsTol',tol);
% 
% %%% Integrating
% [time, Xs_2a] = ode113(@Int_2BI, time, X0, options, u);
% 
% figure; hold all
% plot(r0(1),r0(2),'bo','linewidth',2,'markersize',8)
% plot(rf(1),rf(2),'bx','linewidth',2,'markersize',8)
% plot(Xs_2a(:,1),Xs_2a(:,2))
% axis equal
% PlotBoi2('X, $km$','Y, $km$',16,'LaTex')
% 


% ========================================================================
%%% P2-b
% ========================================================================
% -------------------------------------------------
% Setup
% -------------------------------------------------
%%% Propagation time
t0 = 0;
dt = .001;
% tf = 10000
% warning('wrong tf')
% tf = 1000;
time = t0:dt:tf;
t_tgt = 37864;

%%% Initial state
X0 = [r0; v0];

%%% Desired velocity tolerance
vDesTol = 1e-8;

%%% Vehicle acceleration
aMag = 30e-3; % km/s^2

%%% Choosing ode tolerance
tol = 1e-12;

%%% Setting integrator options
options = odeset('RelTol',tol,'AbsTol',tol);
% [time, Xs] = ode45(@Int_2BI_HW3_P2, time, X0, options, u, rf, t_tgt, aMag, vDesTol);

options_event1 = odeset('Events',@event1_2BI_HW3_P2,'RelTol',tol,'AbsTol',tol);
[time1, Xburn, t_ev, X_ev] = ode113(@Int_2BI_HW3_P2, time, X0, options_event1, u, rf, t_tgt, aMag, vDesTol);
t_ev
% options_event2 = odeset('Events',@event2_2BI_HW3_P2,'RelTol',tol,'AbsTol',tol);
[time2, Xs] = ode113(@Int_2BI, [t_ev, tf], X_ev', options, u);

% % options_event2 = odeset('Events',@event2_2BI_HW3_P2,'RelTol',tol,'AbsTol',tol);
% [time3, Xs3] = ode113(@Int_2BI, [0:0.0001:0.5], Xs(end,:)', options, u);

% Xs = [Xburn; Xs];
time = [time1; time2];


% %%% Preallocating
% Xs = zeros(length(time),6);
% nC_APNb = zeros(length(time),1);
% Xs(1,:) = X0';
% 
% %%% Propagating
% engineCutoff = 0;
% for ti = 1:length(time)-1
%     %%% Determining required velocity
%     [vReq, v2, exitflag] = lambertSolver(Xs(ti,1:3)',rf,tf-ti,0,0,u);
%     
%     %%% Determining desired velocity
%     vDes = (vReq - Xs(ti,4:6)');
%     
%     aT = [0; 0; 0];
%     %%% Determining steering direction
%     if norm(vDes) > vDesTol && engineCutoff == 0
%         vDesHat = vDes./norm(vDes);
%         aT = vDesHat.*aMag;
%     elseif norm(vDes) <= vDesTol && engineCutoff == 0
%         989
%         engineCutoff = 1;
%     end
% 
%     %%% Integrating
%     [X_p] = rk4_2BI_u(Xs(ti,:)',dt,u,aT);
% 
%     %%% Storing
%     Xs(ti+1,:) = X_p';
% end

figure; hold all
p1 = plot(Xburn(:,1),Xburn(:,2),'linewidth',4,'color',colors.std.red);
p2 = plot(Xs(:,1),Xs(:,2),'linewidth',2,'color',colors.std.blue);
% plot(Xs3(:,1),Xs3(:,2),'linewidth',2,'color',colors.std.grn);
plot(r0(1),r0(2),'m.','linewidth',2,'markersize',35)
plot(rf(1),rf(2),'mx','linewidth',2,'markersize',14)
axis equal
PlotBoi2('X, $km$','Y, $km$',16,'LaTex')
plotBodyTexture3(bodies.earth.R,[0,0,0],bodies.earth.img)
xlim([-5 1].*1e4)
ylim([-0.6 3.8].*1e4)
% xlim([-4.201 -4.200].*1e4)
% ylim([3.674 3.676].*1e4)
view(0,90)
legend([p1 p2],'Burn', 'Coast')

missDistance = norm(Xs(end,1:3)' - rf)

toc


% ========================================================================
% ========================================================================
%%% Functions
% ========================================================================
% ========================================================================

function [X_p] = rk4_2BI_u(X,dt,u,a)
%%%
%%% Inputs:
%           1) X  - State, [6x1]
%           2) dt - step size, [1x1]
%           3) u - gravitational parameter
%           4) aV - vehicle acceleration
%           5) 
%           6)
%%% Outputs:
%           1) X_p - State at time t+dt, [3x1]
% ========================================================================
%%% EOM
dXdt = @(X, u, a) [X(4);...
                   X(5);...
                   X(6);...
                   -u*X(1)/(sqrt(X(1)^2 + X(2)^2 + X(3)^2)^3) + a(1);...
                   -u*X(2)/(sqrt(X(1)^2 + X(2)^2 + X(3)^2)^3) + a(2);...
                   -u*X(3)/(sqrt(X(1)^2 + X(2)^2 + X(3)^2)^3) + a(3)];     
%                        
%%% Runge Kutta Time!
k1 = dXdt(X, u, a);
X1 = X + k1.*(dt/2);
k2 = dXdt(X1, u, a);
X2 = X + k2.*(dt/2);
k3 = dXdt(X2, u, a);
X3 = X + k3.*dt;
k4 = dXdt(X3, u, a);

X_p = X + (k1 + 2*k2 + 2*k3 + k4)*(dt/6);

end




function [dX] = Int_2BI_HW3_P2(t,X,u,rf,t_tgt,aMag,vDesTol)
%%% For numerical integration of ECI state under influence of standard
%%% 2-body gravity
%%% Inputs:
%          t - time vector
%          X - initial state [6x1]
%          u - gravitational parameter of primary body

%%% Preallocate state output
dX = zeros(6,1);

%%% Distances from body to spacecraft
r = sqrt(X(1)^2 + X(2)^2 + X(3)^2);

%%% Determining required velocity
[vReq, v2, exitflag] = lambertSolver([X(1); X(2); X(3)],rf,t_tgt-t,0,0,u);

%%% Determining desired velocity
vDes = (vReq - [X(4); X(5); X(6)]);

aT = [0; 0; 0];
% norm(vDes)
%%% Determining acceleration
if norm(vDes) > vDesTol
    vDesHat = vDes./norm(vDes);
    aT = vDesHat.*aMag;
end

%%% Equations of Motion
ddx = X(1)*(-u/(r^3)) + aT(1);
ddy = X(2)*(-u/(r^3)) + aT(2);
ddz = X(3)*(-u/(r^3)) + aT(3);

%%% Output the derivative of the state
dX(1:3) = [X(4); X(5); X(6)]; % km/s
dX(4:6) = [ddx; ddy; ddz]; % km/s^2

end

function [dX] = Int_2BI_HW3_P22(t,X,u,rf)
%%% For numerical integration of ECI state under influence of standard
%%% 2-body gravity
%%% Inputs:
%          t - time vector
%          X - initial state [6x1]
%          u - gravitational parameter of primary body

%%% Preallocate state output
dX = zeros(6,1);

%%% Distances from body to spacecraft
r = sqrt(X(1)^2 + X(2)^2 + X(3)^2);

%%% Equations of Motion
ddx = X(1)*(-u/(r^3));
ddy = X(2)*(-u/(r^3));
ddz = X(3)*(-u/(r^3));

%%% Output the derivative of the state
dX(1:3) = [X(4); X(5); X(6)]; % km/s
dX(4:6) = [ddx; ddy; ddz]; % km/s^2

end



function [value, isterminal, direction] = event1_2BI_HW3_P2(t,X,u,rf,t_tgt,aMag,vDesTol)

%%% Determining required velocity
[vReq, v2, exitflag] = lambertSolver([X(1); X(2); X(3)],rf,t_tgt-t,0,0,u);

%%% Determining desired velocity
vDes = (vReq - [X(4); X(5); X(6)]);

value = norm(vDes) - vDesTol; % When the surface is impacted
isterminal = 1; % stops the integration
direction = -1; % negative direction only
end


% 
% function [value, isterminal, direction] = event2_2BI_HW3_P2(t,X,u,rf)
% dist = norm([X(1); X(2); X(3)] - rf) - 0.5;
% if t > 37000
%     dist
% end
% value = dist; % When the surface is impacted
% isterminal = 1; % stops the integration
% direction = -1; % negative direction only
% end
% 
% 



























