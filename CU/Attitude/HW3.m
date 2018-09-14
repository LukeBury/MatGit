clear
clc
close all
addpath(genpath('/Users/lukebury/Documents/MATLAB/mbin'))
% ------------------------------------------------------------------------
%%% 3.23
% ------------------------------------------------------------------------
% q = [0.5; -0.2; 0.8];
% Q = [0 -q(3) q(2); q(3) 0 -q(1); -q(2) q(1) 0]
% I = eye(3);
% 
% %Cayley
% C1 = inv(I + Q)*(I - Q)
% 
% %Eqn 3.120
% C2 = (1/(1+q'*q))*((1-q'*q)*I + 2*q*q' - 2*Q)
% 
% %Compare
% round(C1,14)==round(C2,14)



% ------------------------------------------------------------------------
%%% 3.28
% ------------------------------------------------------------------------
% clear
% clc
% close all
% t0 = 0; dt = .0001; tf = 5;
% time = dt:dt:tf;
% sigs = zeros(3,length(time)+1);
% 
% w = [1; 0.5; -0.7]; % rad/s
% 
% switchTime = [];
% switchIndex = [];
% 
% for ti = 1:length(time)
%     % Assigning variables for swiftness!
%     s1 = sigs(1,ti); 
%     s2 = sigs(2,ti);
%     s3 = sigs(3,ti); % so swift!
%     s  = norm(sigs(:,ti));
%     
%     % Creating A matrix and updating MRP values
%     A = [1-s^2+2*s1^2, 2*(s1*s2-s3), 2*(s1*s3+s2);...
%                2*(s2*s1+s3), 1-s^2+2*s2^2, 2*(s2*s3-s1);...
%                2*(s3*s1-s2), 2*(s3*s2+s1), 1-s^2+s*s3^2];
%     rate = .25*A*w;
%     
%     % Updating Sigma & Shadow values
%     sigs(:,ti+1) = sigs(:,ti) + rate.*dt;
%     
%     % Check if need to switch to shadow set
%     if norm(sigs(:,ti+1)) > 1
%         switchIndex = [switchIndex,ti];
%         switchTime = [switchTime, time(ti)];
%         sigs(:,ti+1) = -sigs(:,ti+1)./(norm(sigs(:,ti+1))^2);
%     end
%     
% end
% sigs1 = sigs(:,1:switchIndex);
% time1 = time(1:switchIndex);
% sigs2 = sigs(:,(switchIndex+1):end-1);
% time2 = time((switchIndex+1):end);
% switchTime
% 
% figure; hold all
% plot([time1,nan,time2],[sigs1(1,:),nan,sigs2(1,:)],'linewidth',3)
% plot([time1,nan,time2],[sigs1(2,:),nan,sigs2(2,:)],'linewidth',3)
% plot([time1,nan,time2],[sigs1(3,:),nan,sigs2(3,:)],'linewidth',3)
% PlotBoi2('Time, sec','MRP Value',14)
% legend('\sigma1', '\sigma2', '\sigma3')


% ------------------------------------------------------------------------
%%% 2.12
% ------------------------------------------------------------------------
% clear
% clc
% % --------------------
% %%% Setting up system
% % --------------------
% m1 = 1; m2 = 1; m3 = 2; m4 = 2; 
% M = m1+m2+m3+m4;
% 
% R1 = [1;-1;2];
% R2 = [-1;-3;2];
% R3 = [2;-1;-1];
% R4 = [3;-1;-2];
% 
% Rd1 = [2;1;1];
% Rd2 = [0;-1;1];
% Rd3 = [3;2;-1];
% Rd4 = [0;0;1];
% 
% % --------------------
% %%% Position & velocity of COM
% % --------------------
% Rc = (1/M)*(m1*R1 + m2*R2 + m3*R3 + m4*R4);
% Rdc = (1/M)*(m1*Rd1 + m2*Rd2 + m3*Rd3 + m4*Rd4);
% 
% % --------------------
% %%% Relative positions and velocities
% % --------------------
% r1 = R1-Rc;
% r2 = R2-Rc;
% r3 = R3-Rc;
% r4 = R4-Rc;
% 
% rd1 = Rd1-Rdc;
% rd2 = Rd2-Rdc;
% rd3 = Rd3-Rdc;
% rd4 = Rd4-Rdc;
% % --------------------
% %%% Calculations
% % --------------------
% %%% Kinetic Energy
% KE_trans = .5*M*dot(Rdc,Rdc)
% KE_rot = .5*(m1*dot(rd1,rd1) + m2*dot(rd2,rd2) + m3*dot(rd3,rd3) + m4*dot(rd4,rd4))
% 
% %%% Angular Momentum - Origin
% H0 = m1*cross(R1,Rd1) + m2*cross(R2,Rd2) + m3*cross(R3,Rd3) + m4*cross(R4,Rd4)
% 
% %%% Angular Momentum - COM
% Hcom = m1*cross(r1,rd1) + m2*cross(r2,rd2) + m3*cross(r3,rd3) + m4*cross(r4,rd4)

% ------------------------------------------------------------------------
%%% 4.1
% ------------------------------------------------------------------------
% syms r1 r2 r3 real
% rt = [0 -r3 r2; r3 0 -r1; -r2 r1 0];
% -rt*rt

% ------------------------------------------------------------------------
%%% 4.3
% ------------------------------------------------------------------------
% I_B = [15 0 0; 0 11 5; 0 5 16]
% [v,e] = eig(I_B);
% 
% I = e
% 
% v'*I_B*v
% v1_B = v(:,1)
% v2_B = v(:,2)
% v3_B = v(:,3)
% 
% b1_B = [0;1;0];
% b2_B = [0;0;1];
% b3_B = [1;0;0];
% 
% NB = [b1_B,b2_B,b3_B];
% 
% v1_N = NB*v1_B
% v2_N = NB*v2_B
% v3_N = NB*v3_B

% ------------------------------------------------------------------------
%%% 4.7
% ------------------------------------------------------------------------
syms m L th g real

u = m*(.5*L*cos(th)*(3*g*sin(th))/(2*L) - .5*L*sin(th)*(3*g*(1-cos(th))/L))/(m*g - .5*m*L*sin(th)*(3*g*sin(th))/(2*L) - .5*m*L*cos(th)*(3*g*(1-cos(th))/L))
simplify(u)













