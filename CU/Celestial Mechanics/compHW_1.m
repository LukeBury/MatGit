clear
clc
close all
addpath('/Users/lukebury/Documents/MATLAB/mbin')
addpath('/Users/lukebury/Documents/MATLAB/mbin/plot')
addpath('/Users/lukebury/Documents/MATLAB/mbin/matrix')

% Run switch
run_1a = 0;
run_1b = 0;
run_1c = 1;

% ========================================================================
%%% Problem 1
% ========================================================================
%%% Integrator options
% Setting integrator options
tol = 1e-6;
options = odeset('RelTol',tol,'AbsTol',tol);

% -------------------------------------------------------------
%%% 1a - Energy and Angular Momentum Analysis
% -------------------------------------------------------------
if run_1a == 1
% Setting timing options
t0 = 0;
dt = .01;
tf = 14;
tVec = t0:dt:tf;

%%% Setup
r1 = [5,7,10];
r2 = [10,13,2];
r3 = [0, 8, 10];
r12 = r2-r1;
r23 = r3-r2;
r13 = r3-r1;
m = [6; 3; 2]; m = m/sum(m);
G = 1;

%%% Converting to Jacobi Coordinates
R0 = r12;
P0 = (m(1)*r13 + m(2)*r23)/(m(1)+m(2));
Rdot0 = [0,0,0];
Pdot0 = [0,0,0];

%%% Setting intial state vector
X0 = [R0 P0 Rdot0 Pdot0];

%%% Propagating the State
[times,states] = ode45(@jacobiCoord_3Bp_Int,tVec,X0,options,m,G);

%%% Unpacking Results
R = states(:,1:3);    % Position
P = states(:,4:6);    % Position
Rdot = states(:,7:9); % Velocity
Pdot = states(:,10:12); % Velocity

%%% Calculating Energy and Angular Momentum
[H] = jacobiCoord_3Bp_AngMom(R,P,Rdot,Pdot,m);
[E] = jacobiCoord_3Bp_Energy(R,P,Rdot,Pdot,m,G);

figure;
subplot(2,1,1)
plot(times,rownorm(H),'linewidth',2); PlotBoi2('','Angular Momentum',14)
subplot(2,1,2)
plot(times,E-E(1),'linewidth',2); PlotBoi2('Time','Energy Deviation',14)

end % if run_1a == 1

% -------------------------------------------------------------
%%% 1b - Equilateral triangle, m = 1/2, 1/3, 1/6
% -------------------------------------------------------------
if run_1b == 1
% Setting timing options
t0 = 0;
dt = 1;
tf = 120;
tVec = t0:dt:tf;

%%% Setup (equillateral)
% Creating triangle
r1_0 = [0,0,0];
r2_0 = [.5, sqrt(3)/2, 0];
r3_0 = [1,0,0];
% Centering at Origin
center = (r1_0 + r2_0 + r3_0)./3;
r1_0 = r1_0 - center;
r2_0 = r2_0 - center;
r3_0 = r3_0 - center;
% Relative positions
r12_0 = r2_0 - r1_0;
r23_0 = r3_0 - r2_0;
r13_0 = r3_0 - r1_0;

m = [1/2; 1/3; 1/6];
G = 1;
w = [0, 0, 1]; % angular velocity

% Velocities
v1_0 = cross(w,r1_0);
v2_0 = cross(w,r2_0);
v3_0 = cross(w,r3_0);
v12_0 = v2_0 - v1_0;
v23_0 = v3_0 - v2_0;
v13_0 = v3_0 - v1_0;

%%% Converting to Jacobi Coordinates
R0 = r12_0;
P0 = (m(1)*r13_0 + m(2)*r23_0)/(m(1)+m(2));
Rdot0 = v12_0;
Pdot0 = (m(1)*v13_0 + m(2)*v23_0)/(m(1)+m(2));

%%% Setting intial state vector
X0 = [R0 P0 Rdot0 Pdot0];

%%% Propagating the State
[times,states] = ode45(@jacobiCoord_3Bp_Int,tVec,X0,options,m,G);

%%% Unpacking Results
R = states(:,1:3);    % Position
P = states(:,4:6);    % Position
Rdot = states(:,7:9); % Velocity
Pdot = states(:,10:12); % Velocity

%%% Calculating Energy and Angular Momentum
[H] = jacobiCoord_3Bp_AngMom(R,P,Rdot,Pdot,m);
[E] = jacobiCoord_3Bp_Energy(R,P,Rdot,Pdot,m,G);
% 
% %%% Taking magnitude of Angular Momentum
% Hmag = rownorm(H);
% 
% %%% Plotting Hmag, E, and R/P Magnitudes
% figure; subplot(3,1,1)
% plot(times,Hmag-Hmag(1),'linewidth',2); 
% PlotBoi2('','\DeltaH',14)
% title('Relative Equilibrium Lagrange Solution')
% 
% subplot(3,1,2)
% plot(times,E - E(1),'linewidth',2); 
% PlotBoi2('','\DeltaE',14)

% subplot(3,1,3); hold all
figure; hold all
title('Relative Equilibrium Lagrange Solution')
p1 = plot(times,rownorm(R),'linewidth',2); 
p2 = plot(times,rownorm(P),'linewidth',2);
legend([p1 p2],'R magnitude','P magnitude')
PlotBoi2('Time','R, P Magnitude',14)

%%% ----------------------
%%% Repeating 1b w/ nudge
%%% ----------------------
%%% Nudging r3
nudge = 0.001;
r3_0 = r3_0 + [nudge,0,0];
r13_0 = r3_0 - r1_0;

%%% Converting to Jacobi Coordinates
P0 = (m(1)*r13_0 + m(2)*r23_0)/(m(1)+m(2));

%%% Setting intial state vector
X0 = [R0 P0 Rdot0 Pdot0];

%%% Propagating the State
[times,states] = ode45(@jacobiCoord_3Bp_Int,tVec,X0,options,m,G);

%%% Unpacking Results
R = states(:,1:3);    % Position
P = states(:,4:6);    % Position
Rdot = states(:,7:9); % Velocity
Pdot = states(:,10:12); % Velocity

%%% Calculating Energy and Angular Momentum
[H] = jacobiCoord_3Bp_AngMom(R,P,Rdot,Pdot,m);
[E] = jacobiCoord_3Bp_Energy(R,P,Rdot,Pdot,m,G);

%%% Taking magnitude of Angular Momentum
Hmag = rownorm(H);

%%% Plotting Hmag, E, and R/P Magnitudes
% figure; subplot(3,1,1)
% plot(times,Hmag-Hmag(1),'linewidth',2); 
% PlotBoi2('','\DeltaH',14)
% title('Nudged System')
% 
% subplot(3,1,2)
% plot(times,E-E(1),'linewidth',2); 
% PlotBoi2('','\DeltaE',14)

% subplot(3,1,3); hold all
figure; hold all
title('Nudged System')
p1 = plot(times,rownorm(R),'linewidth',2); 
p2 = plot(times,rownorm(P),'linewidth',2);
legend([p1 p2],'R magnitude','P magnitude')
PlotBoi2('Time','R, P magnitude',14)

end % if run_1b == 1

% -------------------------------------------------------------
%%% 1c - 
% -------------------------------------------------------------
if run_1c == 1
% Setting timing options
t0 = 0;
dt = 1;
tf = 200;
tVec = t0:dt:tf;

% Setting integrator options
tol = 1e-7;
options = odeset('RelTol',tol,'AbsTol',tol);

%%% Setup (equillateral)
% Creating triangle
r1_0 = [0,0,0];
r2_0 = [.5, sqrt(3)/2, 0];
r3_0 = [1,0,0];
% Centering at Origin
center = (r1_0 + r2_0 + r3_0)./3;
r1_0 = r1_0 - center;
r2_0 = r2_0 - center;
r3_0 = r3_0 - center;
% Relative positions
r12_0 = r2_0 - r1_0;
r23_0 = r3_0 - r2_0;
r13_0 = r3_0 - r1_0;

%%% Setting masses
m = [0.97; 0.02; .01];
if (m(1)*m(2) + m(2)*m(3) + m(1)*m(3)) >= 1/27
    warning('Routh Stability Criterion not met!')
end
G = 1;
w = [0, 0, 1]; % angular velocity

% Velocities
v1_0 = cross(w,r1_0);
v2_0 = cross(w,r2_0);
v3_0 = cross(w,r3_0);
v12_0 = v2_0 - v1_0;
v23_0 = v3_0 - v2_0;
v13_0 = v3_0 - v1_0;

%%% Converting to Jacobi Coordinates
R0 = r12_0;
P0 = (m(1)*r13_0 + m(2)*r23_0)/(m(1)+m(2));
Rdot0 = v12_0;
Pdot0 = (m(1)*v13_0 + m(2)*v23_0)/(m(1)+m(2));

%%% Setting intial state vector
X0 = [R0 P0 Rdot0 Pdot0];

%%% Propagating the State
[times,states] = ode45(@jacobiCoord_3Bp_Int,tVec,X0,options,m,G);

%%% Unpacking Results
R = states(:,1:3);    % Position
P = states(:,4:6);    % Position
Rdot = states(:,7:9); % Velocity
Pdot = states(:,10:12); % Velocity

%%% Calculating Energy and Angular Momentum
[H] = jacobiCoord_3Bp_AngMom(R,P,Rdot,Pdot,m);
[E] = jacobiCoord_3Bp_Energy(R,P,Rdot,Pdot,m,G);

%%% Taking magnitude of Angular Momentum
Hmag = rownorm(H);

%%% Plotting Hmag, E, and R/P Magnitudes
% figure; subplot(3,1,1)
% plot(times,Hmag-Hmag(1),'linewidth',2); 
% PlotBoi2('','\DeltaH',14)
% title('System satisfying Routh Stability Criterion')
% 
% subplot(3,1,2)
% plot(times,E-E(1),'linewidth',2); 
% PlotBoi2('','\DeltaE',14)

% subplot(3,1,3); hold all
figure; hold all
title('System satisfying Routh Stability Criterion')
p1 = plot(times,rownorm(R),'linewidth',2); 
p2 = plot(times,rownorm(P),'linewidth',2);
legend([p1 p2],'R magnitude','P magnitude')
PlotBoi2('Time','R, P magnitude',14)

%%% ----------------------
%%% Repeating 1c w/ nudge
%%% ----------------------
%%% Nudging r3
nudge = 0.002;
r3_0 = r3_0 + [nudge,0,0];
r13_0 = r3_0 - r1_0;

%%% Converting to Jacobi Coordinates
P0 = (m(1)*r13_0 + m(2)*r23_0)/(m(1)+m(2));

%%% Setting intial state vector
X0 = [R0 P0 Rdot0 Pdot0];

%%% Propagating the State
[times,states] = ode45(@jacobiCoord_3Bp_Int,tVec,X0,options,m,G);

%%% Unpacking Results
R = states(:,1:3);    % Position
P = states(:,4:6);    % Position
Rdot = states(:,7:9); % Velocity
Pdot = states(:,10:12); % Velocity

%%% Calculating Energy and Angular Momentum
[H] = jacobiCoord_3Bp_AngMom(R,P,Rdot,Pdot,m);
[E] = jacobiCoord_3Bp_Energy(R,P,Rdot,Pdot,m,G);

%%% Taking magnitude of Angular Momentum
Hmag = rownorm(H);

%%% Plotting Hmag, E, and R/P Magnitudes
% figure; subplot(3,1,1)
% plot(times,Hmag-Hmag(1),'linewidth',2); 
% PlotBoi2('','\DeltaH',14)
% title(sprintf('Nudged system (r13 + %0.4f)',nudge))
% 
% subplot(3,1,2)
% plot(times,E-E(1),'linewidth',2); 
% PlotBoi2('','\DeltaE',14)

% subplot(3,1,3); hold all
figure; hold all
title(sprintf('Nudged system (r13 + %0.4f)',nudge))
p1 = plot(times,rownorm(R),'linewidth',2); 
p2 = plot(times,rownorm(P),'linewidth',2);
legend([p1 p2],'R magnitude','P magnitude')
PlotBoi2('Time','R, P magnitude',14)

end % if run_1c == 1








% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% ++++++++++++++++++++++++++++ Functions ++++++++++++++++++++++++++++++++++
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


function [dX] = jacobiCoord_3Bp_Int(t,X,m,G)
    % *******************************************************
    % Equations of motion relative to an Inertial Frame
    % *******************************************************
 
    %%% Unpack the state vector
    R = X(1:3);    % Position
    P = X(4:6);    % Position
    Rdot = X(7:9); % Velocity
    Pdot = X(10:12); % Velocity
    m1 = m(1); m2 = m(2); m3 = m(3);

    %%% Relate r13 and r23
    r13 = P + R*m2/(m1+m2);
    r23 = P - R*m1/(m1+m2);

    %%% Equations of Motion
    ddR = -G*(m1+m2)*R/(norm(R)^3) - G*m3*(r13/(norm(r13)^3) - r23/(norm(r23)^3));
    ddP = -G*(m1+m2+m3)*(r23*m2/((m1+m2)*(norm(r23)^3)) + r13*m1/((m1+m2)*(norm(r13)^3)));

    %%% Storing EOMs
    ddX = [ddR; ddP];
    
    %%% Output the derivative of the state
    dX = zeros(12,1);
    dX(1:6) = [Rdot; Pdot];
    dX(7:12) = ddX;
    
end

function [E] = jacobiCoord_3Bp_Energy(R,P,Rdot,Pdot,m,G)
m1 = m(1); m2 = m(2); m3 = m(3);
E = zeros(length(R),1);
for kk = 1:length(R)
    %%% Calculating kinetic energy
    T = dot(Rdot(kk,:),Rdot(kk,:))*m1*m2/(2*(m1+m2)) + dot(Pdot(kk,:),Pdot(kk,:))*m3*(m1+m2)/(2*(m1+m2+m3));

    %%% Converting to relative coordinates
    r12 = norm(R(kk,:));
    r23 = norm(P(kk,:) - R(kk,:)*m1/(m1+m2));
    r31 = norm(-(P(kk,:) + R(kk,:)*m2/(m1+m2)));

    % Calculating Potential Energy
    U = -G*(m1*m2/r12 + m2*m3/r23 + m3*m1/r31);
    
    % Total energy
    E(kk) = T + U;
end
end

function [H] = jacobiCoord_3Bp_AngMom(R,P,Rdot,Pdot,m)
m1 = m(1); m2 = m(2); m3 = m(3);
H = zeros(length(R),3);
for kk = 1:length(R)
    H(kk,:) = cross(R(kk,:),Rdot(kk,:))*m1*m2/(m1+m2) + cross(P(kk,:),Pdot(kk,:))*m3*(m1+m2)/(m1+m2+m3);
end
end

