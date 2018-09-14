clear
clc
close all
addpath(genpath('/Users/lukebury/Documents/MATLAB/mbin'))

% ------------------------------------------------------------------------
%%% P3 - 8.13
% ------------------------------------------------------------------------
I1 = 10; I2 = 20; I3 = 30; I = blkdiag(I1,I2,I3);

%%% Finding initial values
P1 = I1/50;
K = P1^2/I1;
P2 = sqrt(K*I2);
P3 = sqrt(K*I3);
P = blkdiag(P1,P2,P3);

%%% Setting intial conditions
s0 = [0.8; 0.1; -0.1]; % MRPs
w0 = [0;0;0]; % rad/s

%%% Initial state
X0 = [s0;w0];

%%% Time vector
t0 = 0; % sec
dt = .1; % sec
tf = 1000; % sec
times = t0:dt:tf;

%%% Preallocating
states_MRPw_lin = zeros(6,length(times));

%%% Storing ICs
states_MRPw_lin(:,1) = X0;

%%% Integrating
for kk = 1:length(times)-1    
    %%% Checking MRPs for shadow set switch
    s = norm(states_MRPw_lin(1:3,kk));
    if s > 1
        states_MRPw_lin(1:3,kk) = -states_MRPw_lin(1:3,kk)./(s^2);
    end

    [X_p] = rk4_MRPw_lin(states_MRPw_lin(:,kk),dt,I,K,P);
    states_MRPw_lin(:,kk+1) = X_p;

end
K
P

figure('pos',[10 900 850 350]); 
subplot(1,2,1); hold all
plot(times,states_MRPw_lin(1,:),'.')
plot(times,states_MRPw_lin(2,:),'.')
plot(times,states_MRPw_lin(3,:),'.')
[legh,objh] = legend('\sigma1', '\sigma2', '\sigma3');
lineh = findobj(objh,'type','line');
set(lineh,'linestyle','-');
ylim([-1 1])
PlotBoi2('Time, sec','MRP Values',16,'Latex')

subplot(1,2,2); hold all
plot(times,colnorm(states_MRPw_lin(4:6,:)),'linewidth',2)
set(gca,'yscale','log')
PlotBoi2('Time, sec','Control Error, rad/s',16,'Latex')











%%% Preallocating
states_MRPw = zeros(6,length(times));
controls_MRPw = zeros(3,length(times));

%%% Storing ICs
states_MRPw(:,1) = X0;

%%% Integrating
for kk = 1:length(times)-1    
    %%% Checking MRPs for shadow set switch
    s = norm(states_MRPw(1:3,kk));
    if s > 1
        states_MRPw(1:3,kk) = -states_MRPw(1:3,kk)./(s^2);
    end
    
    %%% Computing controls
    w = states_MRPw(4:6,kk);
    wTilde = [0 -w(3) w(2); w(3) 0 -w(1); -w(2) w(1) 0];
    controls_MRPw(:,kk) = -K*states_MRPw(1:3,kk) - P*w + wTilde*I*w;
    
    [X_p] = rk4_MRPw(states_MRPw(:,kk),dt,I,controls_MRPw(:,kk));
    states_MRPw(:,kk+1) = X_p;

end
K
P

figure('pos',[10 900 850 350]); 
subplot(1,2,1); hold all
plot(times,states_MRPw(1,:),'.')
plot(times,states_MRPw(2,:),'.')
plot(times,states_MRPw(3,:),'.')
[legh,objh] = legend('\sigma1', '\sigma2', '\sigma3');
lineh = findobj(objh,'type','line');
set(lineh,'linestyle','-');
ylim([-1 1])
PlotBoi2('Time, sec','MRP Values',16,'Latex')

subplot(1,2,2); hold all
plot(times,colnorm(states_MRPw(4:6,:)),'linewidth',2)
set(gca,'yscale','log')
PlotBoi2('Time, sec','Control Error, rad/s',16,'Latex')


%%% Plotting differences
figure('pos',[10 900 850 350]); 
subplot(3,2,1); hold all
plot(times,states_MRPw(1,:),'b','linewidth',2)
plot(times,states_MRPw_lin(1,:),'b--','linewidth',2)
PlotBoi2('','$\sigma_1$',16,'Latex')
legend('Nonlinear','Linear')
subplot(3,2,3); hold all
plot(times,states_MRPw(2,:),'r','linewidth',2)
plot(times,states_MRPw_lin(2,:),'r--','linewidth',2)
PlotBoi2('','$\sigma_2$',16,'Latex')
legend('Nonlinear','Linear')
subplot(3,2,5); hold all
plot(times,states_MRPw(3,:),'color',[0.9290,0.6940,0.1250],'linewidth',2)
plot(times,states_MRPw_lin(3,:),'--','color',[0.9290,0.6940,0.1250],'linewidth',2)
PlotBoi2('Time, sec','$\sigma_3$',16,'Latex')
legend('Nonlinear','Linear')

subplot(3,2,[2,4,6]); hold all
plot(times,colnorm(states_MRPw(4:6,:)-states_MRPw_lin(4:6,:)),'linewidth',2)
set(gca,'yscale','log')
PlotBoi2('Time, sec','$\Delta$Control Error, rad/s',16,'Latex')




function [X_p] = rk4_MRPw_lin(X,dt,I,K,P)
%%%
%%% Inputs:
%           1) X  - State, [6x1]
%           2) dt - step size, [1x1]
%           3) I  - Inertia tensor, [3x3]
%           4) K  - [1x1]
%           5) P  - [3x3]
%%% Outputs:
%           1) X_p - State at time t+dt, [nx1]
% ========================================================================
%%% EOM
dXdt = @(s,w,I,K,P) [.25*w(1);... % dMRP-1
                   .25*w(2);... % dMRP-2
                   .25*w(3);... % dMRP-3
                   -K*s(1)/I(1,1) - P(1,1)*w(1)/I(1,1);...% dw-1
                   -K*s(2)/I(2,2) - P(2,2)*w(2)/I(2,2);...% dw-2
                   -K*s(3)/I(3,3) - P(3,3)*w(3)/I(3,3)];   % dw-3

%%% Runge Kutta Time!
k1 = dXdt(X(1:3),X(4:6),I,K,P);
X1 = X + k1.*(dt/2);
k2 = dXdt(X1(1:3),X1(4:6),I,K,P);
X2 = X + k2.*(dt/2);
k3 = dXdt(X2(1:3),X2(4:6),I,K,P);
X3 = X + k3.*dt;
k4 = dXdt(X3(1:3),X3(4:6),I,K,P);

X_p = X + (k1 + 2*k2 + 2*k3 + k4)*(dt/6);

end

function [X_p] = rk4_MRPw(X,dt,I,u)
%%%
%%% Inputs:
%           1) X  - State, [6x1]
%           2) dt - step size, [1x1]
%           3) I  - Inertia tensor, [3x3]
%%% Outputs:
%           1) X_p - State at time t+dt, [nx1]
% ========================================================================
%%% EOM
dXdt = @(s,w,I,u) [0.25 * [1-norm(s)^2+2*s(1)^2, 2*(s(1)*s(2)-s(3)), 2*(s(1)*s(3)+s(2))] * w;... % dMRP-1
                   0.25 * [2*(s(2)*s(1)+s(3)), 1-norm(s)^2+2*s(2)^2, 2*(s(2)*s(3)-s(1))] * w;... % dMRP-2
                   0.25 * [2*(s(3)*s(1)-s(2)), 2*(s(3)*s(2)+s(1)), 1-norm(s)^2+2*s(3)^2] * w;... % dMRP-3
                   -(I(3,3)-I(2,2))*w(2)*w(3)/I(1,1) + u(1)/I(1,1);...             % dw-1
                   -(I(1,1)-I(3,3))*w(1)*w(3)/I(2,2) + u(2)/I(2,2);...             % dw-2
                   -(I(2,2)-I(1,1))*w(1)*w(2)/I(3,3) + u(3)/I(3,3)];               % dw-3

%%% Runge Kutta Time!
k1 = dXdt(X(1:3),X(4:6),I,u);
X1 = X + k1.*(dt/2);
k2 = dXdt(X1(1:3),X1(4:6),I,u);
X2 = X + k2.*(dt/2);
k3 = dXdt(X2(1:3),X2(4:6),I,u);
X3 = X + k3.*dt;
k4 = dXdt(X3(1:3),X3(4:6),I,u);

X_p = X + (k1 + 2*k2 + 2*k3 + k4)*(dt/6);

end