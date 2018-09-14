clear
clc
close all
addpath(genpath('/Users/lukebury/Documents/MATLAB/mbin'))

% ------------------------------------------------------------------------
%%% P1 - 8.10
% ------------------------------------------------------------------------
torque_ON = 1;

run_MRP         = 0;
run_EulerAngles = 0;
run_Quaternions = 0;
run_IntegralMRP = 1; % torque always on ... has its own gains

%%% Set Torque
if torque_ON == 0
    L = [0;0;0];
%     K = 1;
%     P = eye(3);
    K = 5;
    P = eye(3).*5;
elseif torque_ON == 1
    L = [1;2;-1]; % Nm
%     K = 1;
%     P = eye(3);
    K = 5;
    P = eye(3).*5;
end
% ------------------------
%%% Initial Conditions
% ------------------------
%%% Inertia
I = blkdiag(10,20,30); % km m^2

%%% Angular Velocity
w0 = [0; 2; 0]; % rad/s

%%% Attitude
s0 = [0; 0.8; 0]; % MRPs
[ C0 ] = mrp2DCM( s0 );

%%% Euler Parameters
[ B0 ] = DCM2Quat( C0 );

%%% Euler Angles (Assuming 3-2-1)
[ EA0 ] = DCM2EulerAngles( C0 , '321'); % rad

% % % % ------------------------
% % % %%% Finding EOMs
% % % % ------------------------
% % % %%% EAs
% % % syms th1 th2 th3 w1 w2 w3 B0 B1 B2 B3 real
% % % w = [w1; w2; w3];
% % % 
% % % Bth = (1/cos(th2))*[0, sin(th3), cos(th3);...
% % %                     0, cos(th2)*cos(th3), -cos(th2)*sin(th3);...
% % %                     cos(th2), sin(th2)*sin(th3), sin(th2)*cos(th3)];
% % %                 
% % % EAdot = Bth*w
% % % 
% % % B = [B0, -B1, -B2, -B3;...
% % %      B1, B0, -B3, B2;...
% % %      B2, B3, B0, -B1;...
% % %      B3, -B2, B1, B0];
% % % Bdot = .5*B*[0;w]

% ------------------------
%%% Integrating w MRPs
% ------------------------
if run_MRP == 1
%%% Initial state
X0 = [s0;w0];

%%% Time vector
t0 = 0; % sec
dt = .1; % sec
tf = 300; % sec
times = t0:dt:tf;

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
    controls_MRPw(:,kk) = -K*states_MRPw(1:3,kk) - P*w + wTilde*I*w - L;
    
    [X_p] = rk4_MRPw(states_MRPw(:,kk),dt,I,controls_MRPw(:,kk),L);
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
title(sprintf('External Torque, L = [%1.0f, %1.0f, %1.0f]'' Nm',L(1),L(2),L(3)))
PlotBoi2('Time, sec','Control Error, rad/s',16,'Latex')

end % run_MRP

% ------------------------
%%% Integrating w Euler Parameters
% ------------------------
if run_Quaternions == 1
%%% Initial state
X0 = [B0;w0];

%%% Time vector
t0 = 0; % sec
dt = .1; % sec
tf = 300; % sec
times = t0:dt:tf;

%%% Preallocating
states_EPw = zeros(7,length(times));
controls_EPw = zeros(3,length(times));

%%% Storing ICs
states_EPw(:,1) = X0;

%%% Integrating
for kk = 1:length(times)-1   
    %%% Setting epsilon vector
    eps = states_EPw(2:4,kk);
    
    %%% Computing controls
    w = states_EPw(5:7,kk);
    wTilde = [0 -w(3) w(2); w(3) 0 -w(1); -w(2) w(1) 0];
    controls_EPw(:,kk) = -K*eps - P*w + wTilde*I*w - L;
    
    [X_p] = rk4_EPw(states_EPw(:,kk),dt,I,controls_EPw(:,kk),L);
    states_EPw(:,kk+1) = X_p;

end
K
P

figure('pos',[10 900 850 350]); 
subplot(1,2,1); hold all
plot(times,states_EPw(1,:),'.')
plot(times,states_EPw(2,:),'.')
plot(times,states_EPw(3,:),'.')
plot(times,states_EPw(4,:),'.')
[legh,objh] = legend('\beta0', '\beta1', '\beta2', '\beta3');
lineh = findobj(objh,'type','line');
set(lineh,'linestyle','-');
% ylim([-1 1])
PlotBoi2('Time, sec','Euler Parameter Values',16,'Latex')

subplot(1,2,2); hold all
plot(times,colnorm(states_EPw(5:7,:)),'linewidth',2)
set(gca,'yscale','log')
title(sprintf('External Torque, L = [%1.0f, %1.0f, %1.0f]'' Nm',L(1),L(2),L(3)))
PlotBoi2('Time, sec','Control Error, rad/s',16,'Latex')

end % run_Quaternions

% ------------------------
%%% Integrating w Euler Angles
% ------------------------
if run_EulerAngles == 1
%%% Initial state
X0 = [EA0;w0];

%%% Time vector
t0 = 0; % sec
dt = .1; % sec
tf = 300; % sec
times = t0:dt:tf;

%%% Preallocating
states_EAw = zeros(6,length(times));
controls_EAw = zeros(3,length(times));

%%% Storing ICs
states_EAw(:,1) = X0;

%%% Integrating
for kk = 1:length(times)-1   
    
    %%% Computing [B(theta)]
    Bth = (1/cos(states_EAw(2,kk)))*[0, sin(states_EAw(3,kk)), cos(states_EAw(3,kk));...
                        0, cos(states_EAw(2,kk))*cos(states_EAw(3,kk)), -cos(states_EAw(2,kk))*sin(states_EAw(3,kk));...
                        cos(states_EAw(2,kk)), sin(states_EAw(2,kk))*sin(states_EAw(3,kk)), sin(states_EAw(2,kk))*cos(states_EAw(3,kk))];
    
    %%% Computing controls
    w = states_EAw(4:6,kk);
    wTilde = [0 -w(3) w(2); w(3) 0 -w(1); -w(2) w(1) 0];
    controls_EAw(:,kk) = -Bth*K*states_EAw(1:3,kk) - P*w + wTilde*I*w - L;
    
    [X_p] = rk4_EA321w(states_EAw(:,kk),dt,I,controls_EAw(:,kk),L);
    states_EAw(:,kk+1) = X_p;

end
K
P

figure('pos',[10 900 850 350]); 
subplot(1,2,1); hold all
plot(times,states_EAw(1,:),'.')
plot(times,states_EAw(2,:),'.')
plot(times,states_EAw(3,:),'.')
[legh,objh] = legend('\theta1', '\theta2', '\theta3');
lineh = findobj(objh,'type','line');
set(lineh,'linestyle','-');
PlotBoi2('Time, sec','Euler Angle Values, rad',16,'Latex')

subplot(1,2,2); hold all
plot(times,colnorm(states_EAw(4:6,:)),'linewidth',2)
set(gca,'yscale','log')
title(sprintf('External Torque, L = [%1.0f, %1.0f, %1.0f]'' Nm',L(1),L(2),L(3)))
PlotBoi2('Time, sec','Control Error, rad/s',16,'Latex')

end % run_EulerAngles

% ------------------------
%%% Integrating MRPs w/ Integral control
% ------------------------
if run_IntegralMRP == 1
%%% Unmodeled torque
L = [1;2;-1]; % Nm
%%% Set initial z
z0 = zeros(3,1);

%%% Initial state
X0 = [s0;w0;z0];

%%% Time vector
t0 = 0; % sec
dt = .1; % sec
tf = 800; % sec
times = t0:dt:tf;

%%% Preallocating
states_MRPwz = zeros(9,length(times));
controls_MRPwz = zeros(3,length(times));

%%% Storing ICs
states_MRPwz(:,1) = X0;

%%% Integrating
for kk = 1:length(times)-1    
    %%% Checking MRPs for shadow set switch
    s = norm(states_MRPwz(1:3,kk));
    if s > 1
        states_MRPwz(1:3,kk) = -states_MRPwz(1:3,kk)./(s^2);
    end
    
    %%% Computing controls
    K = 5;
    P = eye(3).*5;
    KI = eye(3).*0.1;
    w = states_MRPwz(4:6,kk);
    wTilde = [0 -w(3) w(2); w(3) 0 -w(1); -w(2) w(1) 0];
    controls_MRPwz(:,kk) = -K*states_MRPwz(1:3,kk) - P*w + wTilde*I*w - P*KI*states_MRPwz(7:9,kk);
    [X_p] = rk4_MRPwz(states_MRPwz(:,kk),dt,I,controls_MRPwz(:,kk),L,K);
    states_MRPwz(:,kk+1) = X_p;

end
K
P
KI

figure('pos',[10 900 850 350]); 
subplot(1,2,1); hold all
plot(times,states_MRPwz(1,:),'.')
plot(times,states_MRPwz(2,:),'.')
plot(times,states_MRPwz(3,:),'.')
[legh,objh] = legend('\sigma1', '\sigma2', '\sigma3');
lineh = findobj(objh,'type','line');
set(lineh,'linestyle','-');
ylim([-1 1])
PlotBoi2('Time, sec','MRP Values',16,'Latex')

subplot(1,2,2); hold all
plot(times,colnorm(states_MRPwz(4:6,:)),'linewidth',2)
set(gca,'yscale','log')
title(sprintf('External Torque, L = [%1.0f, %1.0f, %1.0f]'' Nm',L(1),L(2),L(3)))
PlotBoi2('Time, sec','Control Error, rad/s',16,'Latex')

figure; hold all
plot(times,states_MRPwz(7,:),'linewidth',2)
plot(times,states_MRPwz(8,:),'linewidth',2)
plot(times,states_MRPwz(9,:),'linewidth',2)
title(sprintf('Unmodeled Torque, L = [%1.0f, %1.0f, %1.0f]'' Nm',L(1),L(2),L(3)))
legend('z_1','z_2','z_3')
PlotBoi2('Time, sec','Z Integral',16,'Latex')
end % run_IntegralMRP




function [X_p] = rk4_MRPw(X,dt,I,u,L)
%%%
%%% Inputs:
%           1) X  - State, [6x1]
%           2) dt - step size, [1x1]
%           3) I  - Inertia tensor, [3x3]
%           4) u  - control input, [3x1]
%           5) L  - External torque, [3x1]
%%% Outputs:
%           1) X_p - State at time t+dt, [nx1]
% ========================================================================
%%% EOM
dXdt = @(s,w,I,u) [0.25 * [1-norm(s)^2+2*s(1)^2, 2*(s(1)*s(2)-s(3)), 2*(s(1)*s(3)+s(2))] * w;... % dMRP-1
                   0.25 * [2*(s(2)*s(1)+s(3)), 1-norm(s)^2+2*s(2)^2, 2*(s(2)*s(3)-s(1))] * w;... % dMRP-2
                   0.25 * [2*(s(3)*s(1)-s(2)), 2*(s(3)*s(2)+s(1)), 1-norm(s)^2+2*s(3)^2] * w;... % dMRP-3
                   -(I(3,3)-I(2,2))*w(2)*w(3)/I(1,1) + u(1)/I(1,1) + L(1)/I(1,1);...             % dw-1
                   -(I(1,1)-I(3,3))*w(1)*w(3)/I(2,2) + u(2)/I(2,2) + L(2)/I(2,2);...             % dw-2
                   -(I(2,2)-I(1,1))*w(1)*w(2)/I(3,3) + u(3)/I(3,3) + L(3)/I(3,3)];               % dw-3

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


function [X_p] = rk4_EPw(X,dt,I,u,L)
%%%
%%% Inputs:
%           1) X  - State, [7x1]
%           2) dt - step size, [1x1]
%           3) I  - Inertia tensor, [3x3]
%           4) u  - control input, [3x1]
%           5) L  - External torque, [3x1]
%%% Outputs:
%           1) X_p - State at time t+dt, [nx1]
% ========================================================================
%%% EOM
dXdt = @(B,w,I,u) [-(B(2)*w(1))/2 - (B(3)*w(2))/2 - (B(4)*w(3))/2;...                % dQuat-1
                   (B(1)*w(1))/2 + (B(3)*w(3))/2 - (B(4)*w(2))/2;...                 % dQuat-2
                   (B(1)*w(2))/2 - (B(2)*w(3))/2 + (B(4)*w(1))/2;...                 % dQuat-3
                   (B(1)*w(3))/2 + (B(2)*w(2))/2 - (B(3)*w(1))/2;...                 % dQuat-4
                   -(I(3,3)-I(2,2))*w(2)*w(3)/I(1,1) + u(1)/I(1,1) + L(1)/I(1,1);... % dw-1
                   -(I(1,1)-I(3,3))*w(1)*w(3)/I(2,2) + u(2)/I(2,2) + L(2)/I(2,2);... % dw-2
                   -(I(2,2)-I(1,1))*w(1)*w(2)/I(3,3) + u(3)/I(3,3) + L(3)/I(3,3)];   % dw-3

%%% Runge Kutta Time!
k1 = dXdt(X(1:4),X(5:7),I,u);
X1 = X + k1.*(dt/2);
k2 = dXdt(X1(1:4),X1(5:7),I,u);
X2 = X + k2.*(dt/2);
k3 = dXdt(X2(1:4),X2(5:7),I,u);
X3 = X + k3.*dt;
k4 = dXdt(X3(1:4),X3(5:7),I,u);

X_p = X + (k1 + 2*k2 + 2*k3 + k4)*(dt/6);

end

function [X_p] = rk4_EA321w(X,dt,I,u,L)
%%%
%%% Inputs:
%           1) X  - State, [6x1]
%           2) dt - step size, [1x1]
%           3) I  - Inertia tensor, [3x3]
%           4) u  - control input, [3x1]
%           5) L  - External torque, [3x1]
%%% Outputs:
%           1) X_p - State at time t+dt, [nx1]
% ========================================================================
%%% EOM
dXdt = @(EA,w,I,u) [(w(3)*cos(EA(3)))/cos(EA(2)) + (w(2)*sin(EA(3)))/cos(EA(2));...                             % dEA-1
                   w(2)*cos(EA(3)) - w(3)*sin(EA(3));...                                                        % dEA-2
                   w(1) + (w(3)*cos(EA(3))*sin(EA(2)))/cos(EA(2)) + (w(2)*sin(EA(2))*sin(EA(3)))/cos(EA(2));... % dEA-3
                   -(I(3,3)-I(2,2))*w(2)*w(3)/I(1,1) + u(1)/I(1,1) + L(1)/I(1,1);...                            % dw-1
                   -(I(1,1)-I(3,3))*w(1)*w(3)/I(2,2) + u(2)/I(2,2) + L(2)/I(2,2);...                            % dw-2
                   -(I(2,2)-I(1,1))*w(1)*w(2)/I(3,3) + u(3)/I(3,3) + L(3)/I(3,3)];                              % dw-3                                            % dw-3

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

function [X_p] = rk4_MRPwz(X,dt,I,u,L,K)
%%%
%%% Inputs:
%           1) X  - State, [9x1]
%           2) dt - step size, [1x1]
%           3) I  - Inertia tensor, [3x3]
%           4) u  - control input, [3x1]
%           5) L  - External torque, [3x1]
%           6) K  - Attitude error feedback gain, [1x1]
%%% Outputs:
%           1) X_p - State at time t+dt, [9x1]
% ========================================================================
%%% EOM
dXdt = @(s,w,z,I,u) [0.25 * [1-norm(s)^2+2*s(1)^2, 2*(s(1)*s(2)-s(3)), 2*(s(1)*s(3)+s(2))] * w;... % dMRP-1
                   0.25 * [2*(s(2)*s(1)+s(3)), 1-norm(s)^2+2*s(2)^2, 2*(s(2)*s(3)-s(1))] * w;...   % dMRP-2
                   0.25 * [2*(s(3)*s(1)-s(2)), 2*(s(3)*s(2)+s(1)), 1-norm(s)^2+2*s(3)^2] * w;...   % dMRP-3
                   -(I(3,3)-I(2,2))*w(2)*w(3)/I(1,1) + u(1)/I(1,1) + L(1)/I(1,1);...               % dw-1
                   -(I(1,1)-I(3,3))*w(1)*w(3)/I(2,2) + u(2)/I(2,2) + L(2)/I(2,2);...               % dw-2
                   -(I(2,2)-I(1,1))*w(1)*w(2)/I(3,3) + u(3)/I(3,3) + L(3)/I(3,3);...               % dw-3
                   K*s(1) - (I(3,3)-I(2,2))*w(2)*w(3) + u(1) + L(1);...                            % dz-1
                   K*s(2) - (I(1,1)-I(3,3))*w(1)*w(3) + u(2) + L(2);...                            % dz-2
                   K*s(3) - (I(2,2)-I(1,1))*w(1)*w(2) + u(3) + L(3)];                              % dz-3

%%% Runge Kutta Time!
k1 = dXdt(X(1:3),X(4:6),X(7:9),I,u);
X1 = X + k1.*(dt/2);
k2 = dXdt(X1(1:3),X1(4:6),X(7:9),I,u);
X2 = X + k2.*(dt/2);
k3 = dXdt(X2(1:3),X2(4:6),X(7:9),I,u);
X3 = X + k3.*dt;
k4 = dXdt(X3(1:3),X3(4:6),X(7:9),I,u);

X_p = X + (k1 + 2*k2 + 2*k3 + k4)*(dt/6);

end










