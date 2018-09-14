clear
clc
close all
addpath('/Users/lukebury/Documents/MATLAB/mbin')
addpath('/Users/lukebury/Documents/MATLAB/mbin/plot')
addpath('/Users/lukebury/Documents/MATLAB/mbin/matrix')

% Run switch
run_2a = 0;
run_2b = 1;
run_2c = 0;
run_2d = 0;
run_sym = 0;

% ========================================================================
%%% Problem 2
% ========================================================================

% -------------------------------------------------------------
%%% 2a - Energy and Angular Momentum Analysis
% -------------------------------------------------------------
if run_2a == 1
%%% Integrator options
% Setting integrator options
tol = 1e-9;
options = odeset('RelTol',tol,'AbsTol',tol);

% Setting timing options
t0 = 0;
dt = .1;
tf = 13;
tVec = t0:dt:tf;

%%% Setup
r1_0 = [5,-7,1]; r1_0 = r1_0./sum(r1_0);
r2_0 = [19,13,-2]; r2_0 = r2_0./sum(r2_0);
r3_0 = [0, 8, 10]; r3_0 = r3_0./sum(r3_0);
r4_0 = [-3, -9, 5]; r4_0 = r4_0./sum(r4_0);
r5_0 = [2, 7, -7]; r5_0 = r5_0./sum(r5_0);
r12_0 = r2_0 - r1_0;
r23_0 = r3_0 - r2_0;
r34_0 = r4_0 - r3_0;
r45_0 = r5_0 - r4_0;
r15_0 = r5_0 - r1_0;
m = [9; 6; 4; 3; 2]; m = m/sum(m);
m1 = m(1); m2 = m(2); m3 = m(3); m4 = m(4); m5 = m(5);
G = 1;

%%% Converting to Jacobi Coordinates
R1_0 = r12_0;
R2_0 = r23_0 + R1_0*m(1)/(m(1)+m(2));
R3_0 = r34_0 + R2_0*(m(1)+m(2))/(m(1)+m(2)+m(3));
R4_0 = r45_0 + R3_0*(m(1)+m(2)+m(3))/(m(1)+m(2)+m(3)+m(4));
R1dot_0 = [0,0,0];
R2dot_0 = [0,0,0];
R3dot_0 = [0,0,0];
R4dot_0 = [0,0,0];

%%% Setting intial state vector
X0 = [R1_0 R2_0 R3_0 R4_0 R1dot_0 R2dot_0 R3dot_0 R4dot_0];

%%% Propagating the State
[times,states] = ode45(@jacobiCoord_5Bp_Int,tVec,X0,options,m,G);

%%% Unpacking Results
R1 = states(:,1:3);      % Position
R2 = states(:,4:6);      % Position
R3 = states(:,7:9);      % Position
R4 = states(:,10:12);    % Position
R1dot = states(:,13:15); % Velocity
R2dot = states(:,16:18); % Velocity
R3dot = states(:,19:21); % Velocity
R4dot = states(:,22:24); % Velocity

%%% Calculating Energy and Angular Momentum
[H] = jacobiCoord_5Bp_AngMom(R1,R2,R3,R4,R1dot,R2dot,R3dot,R4dot,m);
[E] = jacobiCoord_5Bp_Energy(R1,R2,R3,R4,R1dot,R2dot,R3dot,R4dot,m,G);

%%% Taking magnitude of Angular Momentum
Hmag = rownorm(H);

%%% Plotting Hmag and E percent changes and R/P magnitudes
% figure; subplot(3,1,1)
% plot(times(2:end),percentchange(Hmag(2:end,1)),'linewidth',2); 
% PlotBoi2('','Percent Change, %',14)
% title('Angular Momentum Percent Change')
% 
% subplot(3,1,2)
% plot(times,percentchange(E),'linewidth',2); 
% PlotBoi2('','Percent Change, %',14)
% title('Energy Percent Change')
figure; subplot(3,1,1)
plot(times,Hmag-Hmag(1),'linewidth',2); 
PlotBoi2('','\DeltaH',14)

subplot(3,1,2)
plot(times,E-E(1),'linewidth',2); 
PlotBoi2('','\DeltaE',14)

subplot(3,1,3); hold all
p1 = plot(times,rownorm(R1),'linewidth',2); 
p2 = plot(times,rownorm(R2),'linewidth',2);
p3 = plot(times,rownorm(R3),'linewidth',2); 
p4 = plot(times,rownorm(R4),'linewidth',2);
legend([p1 p2 p3 p4],'R1','R2','R3','R4')
PlotBoi2('Time','R magnitudes',14)

% figure; subplot(2,1,1)
% plot(times,Hmag,'linewidth',2); 
% PlotBoi2('','Ang. Mom. Magnitude',14)
% title('Angular Momentum Percent Change')
% 
% subplot(2,1,2)
% plot(times,E,'linewidth',2); 
% PlotBoi2('','Energy',14)
% title('Energy Percent Change')
% 
% %%% Plotting Jacobi Coordinate Values
% figure; hold all
% plot3(R1(:,1),R1(:,2),R1(:,3),'linewidth',2)
% plot3(R2(:,1),R2(:,2),R2(:,3),'linewidth',2)
% plot3(R3(:,1),R3(:,2),R3(:,3),'linewidth',2)
% plot3(R4(:,1),R4(:,2),R4(:,3),'linewidth',2)
% PlotBoi3('X','Y','Z',14)
% 
% %%% Converting to relative positions
% r12 = R1;
% r13 = R2 + m2*R1/(m1+m2);
% r14 = R3 + m3*R2/(m1+m2+m3) + m2*R1/(m1+m2);
% r15 = R4 + m4*R3/(m1+m2+m3+m4) + m3*R2/(m1+m2+m3) + m2*R1/(m1+m2);
% r23 = R2 - m1*R1/(m1+m2);
% r24 = R3 + m3*R2/(m1+m2+m3) - m1*R1/(m1+m2);
% r25 = R4 + m4*R3/(m1+m2+m3+m4) + m3*R2/(m1+m2+m3) - m1*R1/(m1+m2);
% r34 = R3 - (m1+m2)*R2/(m1+m2+m3);
% r35 = R4 + m4*R3/(m1+m2+m3+m4) - (m1+m2)*R2/(m1+m2+m3);
% r45 = R4 - (m1+m2+m3)*R3/(m1+m2+m3+m4);
% 
% %%% Plotting Relative Positions
% figure; hold all
% plot3(r12(:,1),r12(:,2),r12(:,3),'linewidth',2)
% plot3(r13(:,1),r13(:,2),r13(:,3),'linewidth',2)
% plot3(r14(:,1),r14(:,2),r14(:,3),'linewidth',2)
% plot3(r15(:,1),r15(:,2),r15(:,3),'linewidth',2)
% plot3(r23(:,1),r23(:,2),r23(:,3),'linewidth',2)
% plot3(r24(:,1),r24(:,2),r24(:,3),'linewidth',2)
% plot3(r25(:,1),r25(:,2),r25(:,3),'linewidth',2)
% plot3(r34(:,1),r34(:,2),r34(:,3),'linewidth',2)
% plot3(r35(:,1),r35(:,2),r35(:,3),'linewidth',2)
% plot3(r45(:,1),r45(:,2),r45(:,3),'linewidth',2)
% PlotBoi3('X','Y','Z',14)

end % if run_2a == 1

% -------------------------------------------------------------
%%% 2b - pentagon
% -------------------------------------------------------------
if run_2b == 1
%%% Integrator options
% Setting integrator options
tol = 1e-6;
options = odeset('RelTol',tol,'AbsTol',tol);

% Setting timing options
t0 = 0;
dt = 1;
tf = 140;
tVec = t0:dt:tf;

%%% Setup (pentagon)
c1 = cos(2*pi/5);
c2 = cos(pi/5);
s1 = sin(2*pi/5);
s2 = sin(4*pi/5);
r1_0 = [1,0,0];
r2_0 = [c1,s1,0];
r3_0 = [-c2,s2,0];
r4_0 = [-c2,-s2,0];
r5_0 = [c1,-s1,0];
r12_0 = r2_0 - r1_0;
r23_0 = r3_0 - r2_0;
r34_0 = r4_0 - r3_0;
r45_0 = r5_0 - r4_0;
r15_0 = r5_0 - r1_0;
m = [.2; .2; .2; .2; .2]; m = m/sum(m);
m1 = m(1); m2 = m(2); m3 = m(3); m4 = m(4); m5 = m(5);
G = 1;

%%% Determining circular orbit prms
rp = 1;
e = 0;
ra = -rp*(e+1)/(e-1);
a = (ra+rp)/2;
u = G*sum(m);
v = sqrt(u*(2/rp-1/a));
% w = [0, 0, v];
w = [0, 0, .5247];
v1_0 = cross(w,r1_0);
v2_0 = cross(w,r2_0);
v3_0 = cross(w,r3_0);
v4_0 = cross(w,r4_0);
v5_0 = cross(w,r5_0);
v12_0 = v2_0 - v1_0;
v23_0 = v3_0 - v2_0;
v34_0 = v4_0 - v3_0;
v45_0 = v5_0 - v4_0;
v15_0 = v5_0 - v1_0;


%%% Converting to Jacobi Coordinates
R1_0 = r12_0;
R2_0 = r23_0 + R1_0*m(1)/(m(1)+m(2));
R3_0 = r34_0 + R2_0*(m(1)+m(2))/(m(1)+m(2)+m(3));
R4_0 = r45_0 + R3_0*(m(1)+m(2)+m(3))/(m(1)+m(2)+m(3)+m(4));

R1dot_0 = v12_0;
R2dot_0 = v23_0 + R1dot_0*m(1)/(m(1)+m(2));
R3dot_0 = v34_0 + R2dot_0*(m(1)+m(2))/(m(1)+m(2)+m(3));
R4dot_0 = v45_0 + R3dot_0*(m(1)+m(2)+m(3))/(m(1)+m(2)+m(3)+m(4));

%%% Setting intial state vector
X0 = [R1_0 R2_0 R3_0 R4_0 R1dot_0 R2dot_0 R3dot_0 R4dot_0];

%%% Propagating the State
[times,states] = ode45(@jacobiCoord_5Bp_Int,tVec,X0,options,m,G);

%%% Unpacking Results
R1 = states(:,1:3);      % Position
R2 = states(:,4:6);      % Position
R3 = states(:,7:9);      % Position
R4 = states(:,10:12);    % Position
R1dot = states(:,13:15); % Velocity
R2dot = states(:,16:18); % Velocity
R3dot = states(:,19:21); % Velocity
R4dot = states(:,22:24); % Velocity

%%% Calculating Energy and Angular Momentum
[H] = jacobiCoord_5Bp_AngMom(R1,R2,R3,R4,R1dot,R2dot,R3dot,R4dot,m);
[E] = jacobiCoord_5Bp_Energy(R1,R2,R3,R4,R1dot,R2dot,R3dot,R4dot,m,G);

%%% Taking magnitude of Angular Momentum
Hmag = rownorm(H);

figure; hold all
p1 = plot(times,rownorm(R1),'linewidth',2); 
p2 = plot(times,rownorm(R2),'linewidth',2);
p3 = plot(times,rownorm(R3),'linewidth',2); 
p4 = plot(times,rownorm(R4),'linewidth',2);
legend([p1 p2 p3 p4],'R1','R2','R3','R4')
PlotBoi2('Time','R magnitudes',14)

end % if run_2b == 1


% -------------------------------------------------------------
%%% 2c - eccentric orbits
% -------------------------------------------------------------
if run_2c == 1
%%% Integrator options
% Setting integrator options
tol = 1e-6;
options = odeset('RelTol',tol,'AbsTol',tol);

% Setting timing options
t0 = 0;
dt = 1;
tf = 140;
tVec = t0:dt:tf;

%%% Setup (pentagon)
c1 = cos(2*pi/5);
c2 = cos(pi/5);
s1 = sin(2*pi/5);
s2 = sin(4*pi/5);
r1_0 = [1,0,0];
r2_0 = [c1,s1,0];
r3_0 = [-c2,s2,0];
r4_0 = [-c2,-s2,0];
r5_0 = [c1,-s1,0];
r12_0 = r2_0 - r1_0;
r23_0 = r3_0 - r2_0;
r34_0 = r4_0 - r3_0;
r45_0 = r5_0 - r4_0;
r15_0 = r5_0 - r1_0;
m = [.2; .2; .2; .2; .2]; m = m/sum(m);
m1 = m(1); m2 = m(2); m3 = m(3); m4 = m(4); m5 = m(5);
G = 1;

%%% New stuff
Ip0 = m1*dot(r1_0,r1_0) + m2*dot(r2_0,r2_0) + m3*dot(r3_0,r3_0) + m4*dot(r4_0,r4_0) + m5*dot(r5_0,r5_0);
r13_0 = r12_0 + r23_0;
r14_0 = r13_0 + r34_0;
r24_0 = r23_0 + r34_0;
r25_0 = r24_0 + r45_0;
r35_0 = r34_0 + r45_0;
U0 = -G*(m1*m2/norm(r12_0) + m1*m3/norm(r13_0) + m1*m4/norm(r14_0) + m1*m5/norm(r15_0) + m2*m3/norm(r23_0) + m2*m4/norm(r24_0) + m2*m5/norm(r25_0) + m3*m4/norm(r34_0) + m3*m5/norm(r35_0) + m4*m5/norm(r45_0));
wcirc = sqrt(-U0/Ip0);
e = 0.1;
wecc  = sqrt(wcirc*wcirc*(e+1));

% w = [0, 0, 0.5];
% w = [0, 0, wecc];
w = [0, 0, 0.5503];
v1_0 = cross(w,r1_0);
v2_0 = cross(w,r2_0);
v3_0 = cross(w,r3_0);
v4_0 = cross(w,r4_0);
v5_0 = cross(w,r5_0);
v12_0 = v2_0 - v1_0;
v23_0 = v3_0 - v2_0;
v34_0 = v4_0 - v3_0;
v45_0 = v5_0 - v4_0;
v15_0 = v5_0 - v1_0;

%%% Converting to Jacobi Coordinates
R1_0 = r12_0;
R2_0 = r23_0 + R1_0*m(1)/(m(1)+m(2));
R3_0 = r34_0 + R2_0*(m(1)+m(2))/(m(1)+m(2)+m(3));
R4_0 = r45_0 + R3_0*(m(1)+m(2)+m(3))/(m(1)+m(2)+m(3)+m(4));
R1dot_0 = v12_0;
R2dot_0 = v23_0 + R1dot_0*m(1)/(m(1)+m(2));
R3dot_0 = v34_0 + R2dot_0*(m(1)+m(2))/(m(1)+m(2)+m(3));
R4dot_0 = v45_0 + R3dot_0*(m(1)+m(2)+m(3))/(m(1)+m(2)+m(3)+m(4));

%%% Setting intial state vector
X0 = [R1_0 R2_0 R3_0 R4_0 R1dot_0 R2dot_0 R3dot_0 R4dot_0];

%%% Propagating the State
[times,states] = ode45(@jacobiCoord_5Bp_Int,tVec,X0,options,m,G);

%%% Unpacking Results
R1 = states(:,1:3);      % Position
R2 = states(:,4:6);      % Position
R3 = states(:,7:9);      % Position
R4 = states(:,10:12);    % Position
R1dot = states(:,13:15); % Velocity
R2dot = states(:,16:18); % Velocity
R3dot = states(:,19:21); % Velocity
R4dot = states(:,22:24); % Velocity

%%% Calculating Energy and Angular Momentum
[H] = jacobiCoord_5Bp_AngMom(R1,R2,R3,R4,R1dot,R2dot,R3dot,R4dot,m);
[E] = jacobiCoord_5Bp_Energy(R1,R2,R3,R4,R1dot,R2dot,R3dot,R4dot,m,G);

%%% Taking magnitude of Angular Momentum
Hmag = rownorm(H);

figure; hold all
p1 = plot(times,rownorm(R1),'linewidth',2); 
p2 = plot(times,rownorm(R2),'linewidth',2);
p3 = plot(times,rownorm(R3),'linewidth',2); 
p4 = plot(times,rownorm(R4),'linewidth',2);
legend([p1 p2 p3 p4],'R1','R2','R3','R4')
PlotBoi2('Time','R magnitudes',14)
ylim([0 10])

end % if run_2c == 1


% -------------------------------------------------------------
%%% 2d - 2B & 3B Orbiting Each Other
% -------------------------------------------------------------
if run_2d == 1
%%% Integrator options
% Setting integrator options
tol = 1e-6;
options = odeset('RelTol',tol,'AbsTol',tol);

% Setting timing options
t0 = 0;
dt = .1;
tf = 30;
tVec = t0:dt:tf;

for trial = 1:3
%%% Setup (2B)
r1_0 = [.1,3,0];
r2_0 = [-.1,3,0];
if trial == 1
    v1_0 = [0,-1.4142/2,0] + [0,-.1,0];
    v2_0 = [0,1.4142/2,0]  + [0,-.1,0];
elseif trial == 2
    v1_0 = [0,-1.4142/2,0] + [0,-.4,0];
    v2_0 = [0,1.4142/2,0]  + [0,-.4,0];
elseif trial == 3
    v1_0 = [0,-1.4142/2,0] + [.4,0,0];
    v2_0 = [0,1.4142/2,0]  + [.4,0,0];
end

%%% Setup (3B)
r3_0 = [0,sqrt(3)/4,0];
r4_0 = [-.5,-sqrt(3)/4,0];
r5_0 = [.5,-sqrt(3)/4,0];
w3 = [0,0,.7746];
v3_0 = cross(w3,r3_0);
v4_0 = cross(w3,r4_0);
v5_0 = cross(w3,r5_0);

%%% General Setup
m = [.2; .2; .2; .2; .2]; m = m/sum(m);
m1 = m(1); m2 = m(2); m3 = m(3); m4 = m(4); m5 = m(5);
G = 1;

r12_0 = r2_0 - r1_0;
r23_0 = r3_0 - r2_0;
r34_0 = r4_0 - r3_0;
r45_0 = r5_0 - r4_0;
r15_0 = r5_0 - r1_0;

v12_0 = v2_0 - v1_0;
v23_0 = v3_0 - v2_0;
v34_0 = v4_0 - v3_0;
v45_0 = v5_0 - v4_0;
v15_0 = v5_0 - v1_0;

%%% Converting to Jacobi Coordinates
R1_0 = r12_0;
R2_0 = r23_0 + R1_0*m(1)/(m(1)+m(2));
R3_0 = r34_0 + R2_0*(m(1)+m(2))/(m(1)+m(2)+m(3));
R4_0 = r45_0 + R3_0*(m(1)+m(2)+m(3))/(m(1)+m(2)+m(3)+m(4));
R1dot_0 = v12_0;
R2dot_0 = v23_0 + R1dot_0*m(1)/(m(1)+m(2));
R3dot_0 = v34_0 + R2dot_0*(m(1)+m(2))/(m(1)+m(2)+m(3));
R4dot_0 = v45_0 + R3dot_0*(m(1)+m(2)+m(3))/(m(1)+m(2)+m(3)+m(4));

%%% Setting intial state vector
X0 = [R1_0 R2_0 R3_0 R4_0 R1dot_0 R2dot_0 R3dot_0 R4dot_0];

%%% Propagating the State
[times,states] = ode45(@jacobiCoord_5Bp_Int,tVec,X0,options,m,G);

%%% Unpacking Results
R1 = states(:,1:3);      % Position
R2 = states(:,4:6);      % Position
R3 = states(:,7:9);      % Position
R4 = states(:,10:12);    % Position
R1dot = states(:,13:15); % Velocity
R2dot = states(:,16:18); % Velocity
R3dot = states(:,19:21); % Velocity
R4dot = states(:,22:24); % Velocity

%%% Calculating Energy and Angular Momentum
[H] = jacobiCoord_5Bp_AngMom(R1,R2,R3,R4,R1dot,R2dot,R3dot,R4dot,m);
[E] = jacobiCoord_5Bp_Energy(R1,R2,R3,R4,R1dot,R2dot,R3dot,R4dot,m,G);

%%% Relate little r's
r12 = R1;
r23 = R2 - R1*m1/(m1+m2);
r34 = R3 - R2*(m1+m2)/(m1+m2+m3);
r45 = R4 - R3*(m1+m2+m3)/(m1+m2+m3+m4);
r15 = r12 + r23 + r34 + r45;

figure; hold all
p1 = plot(times,rownorm(r12),'linewidth',2);
p2 = plot(times,rownorm(r23),'linewidth',2);
p3 = plot(times,rownorm(r34),'linewidth',2);
p4 = plot(times,rownorm(r45),'linewidth',2);
p5 = plot(times,rownorm(r15),'linewidth',2);
legend([p1 p2 p3 p4 p5],'r12','r23','r34','r45','r15')
PlotBoi2('Time','r magnitudes',14)
title(sprintf('Case %1d',trial))
end % trials
end % if run_2d == 1












if run_sym == 1
tic %
syms m1 m2 m3 m4 m5 real
syms R1 R2 R3 R4 real
syms r12 r23 r34 r45
syms G real

toc %%
U = -G*(m1*m2/R1 + m1*m3/(R2 + m2*R1/(m1+m2)) + m1*m4/(R3 + m3*R2/(m1+m2+m3) + m2*R1/(m1+m2)) + m1*m5/(R4 + m4*R3/(m1+m2+m3+m4) + m3*R2/(m1+m2+m3) + m2*R1/(m1+m2)) + m2*m3/(R2 - m1*R1/(m1+m2)) + m2*m4/(R3 + m3*R2/(m1+m2+m3) - m1*R1/(m1+m2)) + m2*m5/(R4 + m4*R3/(m1+m2+m3+m4) + m3*R2/(m1+m2+m3) - m1*R1/(m1+m2)) + m3*m4/(R3 - (m1+m2)*R2/(m1+m2+m3)) + m3*m5/(R4 + m4*R3/(m1+m2+m3+m4) - (m1+m2)*R2/(m1+m2+m3)) + m4*m5/(R4 - (m1+m2+m3)*R3/(m1+m2+m3+m4)));

toc %%
dR1 = -diff(U,R1);
dR2 = -diff(U,R2);
dR3 = -diff(U,R3);
dR4 = -diff(U,R4);

toc %%
R1 = r12;
R2 = r23 + R1*m1/(m1+m2);
R3 = r34 + (m1+m2)*R2/(m1+m2+m3);
R4 = r45 + (m1+m2+m3)*R3/(m1+m2+m3+m4);

toc %%
dR1_sub = subs(dR1);
dR2_sub = subs(dR2);
dR3_sub = subs(dR3);
dR4_sub = subs(dR4);

toc %%
sdR1 = simplify(dR1_sub)
sdR2 = simplify(dR2_sub)
sdR3 = simplify(dR3_sub)
sdR4 = dR4_sub
% sdR4 = simplify(dR4_sub)    this gets pretty nasty

toc %%
ddR1 = sdR1*(m1+m2)/(m2*m1)
%-(G*(m1 + m2)*((m1*m2)/r12^2 - (m1*m2*m3)/(r23^2*(m1 + m2)) + (m1*m2*m4)/((m1 + m2)*(r12 + r23 + r34)^2) - (m1*m2*m5)/((m1 + m2)*(r23 + r34 + r45)^2) + (m1*m2*m5)/((m1 + m2)*(r12 + r23 + r34 + r45)^2) + (m1*m2*m3)/((m1 + m2)*(r12 + r23)^2) - (m1*m2*m4)/((m1 + m2)*(r23 + r34)^2)))/(m1*m2)
ddR2 = sdR2*(m1+m2+m3)/(m3*(m1+m2))
%-(G*(m1 + m2 + m3)*((m2*m3)/r23^2 + (m1*m3)/(r12 + r23)^2 - (m3*m5*(m1 + m2))/((r34 + r45)^2*(m1 + m2 + m3)) + (m1*m3*m5)/((m1 + m2 + m3)*(r12 + r23 + r34 + r45)^2) - (m3*m4*(m1 + m2))/(r34^2*(m1 + m2 + m3)) + (m2*m3*m4)/((r23 + r34)^2*(m1 + m2 + m3)) + (m1*m3*m4)/((m1 + m2 + m3)*(r12 + r23 + r34)^2) + (m2*m3*m5)/((m1 + m2 + m3)*(r23 + r34 + r45)^2)))/(m3*(m1 + m2))
ddR3 = sdR3*(m1+m2+m3+m4)/(m4*(m1+m2+m3))
%-(G*(m1 + m2 + m3 + m4)*((m3*m4)/r34^2 + (m1*m4)/(r12 + r23 + r34)^2 + (m2*m4)/(r23 + r34)^2 - (m4*m5*(m1 + m2 + m3))/(r45^2*(m1 + m2 + m3 + m4)) + (m2*m4*m5)/((r23 + r34 + r45)^2*(m1 + m2 + m3 + m4)) + (m1*m4*m5)/((m1 + m2 + m3 + m4)*(r12 + r23 + r34 + r45)^2) + (m3*m4*m5)/((r34 + r45)^2*(m1 + m2 + m3 + m4))))/(m4*(m1 + m2 + m3))
ddR4 = sdR4*(m1+m2+m3+m4+m5)/(m5*(m1+m2+m3+m4))
%-(G*((m4*m5)/r45^2 + (m3*m5)/(r45 - ((m1 + m2)*(r23 + (m1*r12)/(m1 + m2)))/(m1 + m2 + m3) + (m4*(r34 + ((m1 + m2)*(r23 + (m1*r12)/(m1 + m2)))/(m1 + m2 + m3)))/(m1 + m2 + m3 + m4) + ((r34 + ((m1 + m2)*(r23 + (m1*r12)/(m1 + m2)))/(m1 + m2 + m3))*(m1 + m2 + m3))/(m1 + m2 + m3 + m4))^2 + (m1*m5)/(r45 + (m4*(r34 + ((m1 + m2)*(r23 + (m1*r12)/(m1 + m2)))/(m1 + m2 + m3)))/(m1 + m2 + m3 + m4) + (m3*(r23 + (m1*r12)/(m1 + m2)))/(m1 + m2 + m3) + ((r34 + ((m1 + m2)*(r23 + (m1*r12)/(m1 + m2)))/(m1 + m2 + m3))*(m1 + m2 + m3))/(m1 + m2 + m3 + m4) + (m2*r12)/(m1 + m2))^2 + (m2*m5)/(r45 + (m4*(r34 + ((m1 + m2)*(r23 + (m1*r12)/(m1 + m2)))/(m1 + m2 + m3)))/(m1 + m2 + m3 + m4) + (m3*(r23 + (m1*r12)/(m1 + m2)))/(m1 + m2 + m3) + ((r34 + ((m1 + m2)*(r23 + (m1*r12)/(m1 + m2)))/(m1 + m2 + m3))*(m1 + m2 + m3))/(m1 + m2 + m3 + m4) - (m1*r12)/(m1 + m2))^2)*(m1 + m2 + m3 + m4 + m5))/(m5*(m1 + m2 + m3 + m4))
r13 = r12 + r23;
r14 = r13 + r34;
r15 = r14 + r45;
r24 = r23 + r34;
r25 = r24 + r45;
r35 = r34 + r45;
tt = subs(ddR4);
tt
tt_s = simplify(tt)



% r12 = R1;
% r13 = R2 + m2*R1/(m1+m2);
% r14 = R3 + m3*R2/(m1+m2+m3) + m2*R1/(m1+m2);
% r15 = R4 + m4*R3/(m1+m2+m3+m4) + m3*R2/(m1+m2+m3) + m2*R1/(m1+m2);

% r23 = R2 - m1*R1/(m1+m2);
% r24 = R3 + m3*R2/(m1+m2+m3) - m1*R1/(m1+m2);
% r25 = R4 + m4*R3/(m1+m2+m3+m4) + m3*R2/(m1+m2+m3) - m1*R1/(m1+m2);

% r34 = R3 - (m1+m2)*R2/(m1+m2+m3);
% r35 = R4 + m4*R3/(m1+m2+m3+m4) - (m1+m2)*R2/(m1+m2+m3);

% r45 = R4 - (m1+m2+m3)*R3/(m1+m2+m3+m4);
end






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


function [dX] = jacobiCoord_5Bp_Int(t,X,m,G)
    % *******************************************************
    % Equations of motion relative to an Inertial Frame
    % *******************************************************
 
    %%% Unpack the state vector
    R1 = X(1:3);      % Position
    R2 = X(4:6);      % Position
    R3 = X(7:9);      % Position
    R4 = X(10:12);    % Position
    R1dot = X(13:15); % Velocity
    R2dot = X(16:18); % Velocity
    R3dot = X(19:21); % Velocity
    R4dot = X(22:24); % Velocity
    m1 = m(1); m2 = m(2); m3 = m(3); m4 = m(4); m5 = m(5);

    %%% Relate little r's
    r12 = R1;
    r23 = R2 - R1*m1/(m1+m2);
    r34 = R3 - R2*(m1+m2)/(m1+m2+m3);
    r45 = R4 - R3*(m1+m2+m3)/(m1+m2+m3+m4);
    r15 = r12 + r23 + r34 + r45;
%     r15 = R4 + R1*m2/(m1+m2) + R2*m3/(m1+m2+m3) + R3*m4/(m1+m2+m3+m4);

    %%% Equations of Motion
% % %     ddR1 = (G*(m1 + m2)*((m1*m2)/r12^2 - (m1*m2*m3)/(r23^2*(m1 + m2)) + (m1*m2*m4)/((m1 + m2)*(r12 + r23 + r34)^2) - (m1*m2*m5)/((m1 + m2)*(r23 + r34 + r45)^2) + (m1*m2*m5)/((m1 + m2)*(r12 + r23 + r34 + r45)^2) + (m1*m2*m3)/((m1 + m2)*(r12 + r23)^2) - (m1*m2*m4)/((m1 + m2)*(r23 + r34)^2)))/(m1*m2);
% % %     ddR2 = (G*(m1 + m2 + m3)*((m2*m3)/r23^2 + (m1*m3)/(r12 + r23)^2 - (m3*m5*(m1 + m2))/((r34 + r45)^2*(m1 + m2 + m3)) + (m1*m3*m5)/((m1 + m2 + m3)*(r12 + r23 + r34 + r45)^2) - (m3*m4*(m1 + m2))/(r34^2*(m1 + m2 + m3)) + (m2*m3*m4)/((r23 + r34)^2*(m1 + m2 + m3)) + (m1*m3*m4)/((m1 + m2 + m3)*(r12 + r23 + r34)^2) + (m2*m3*m5)/((m1 + m2 + m3)*(r23 + r34 + r45)^2)))/(m3*(m1 + m2));
% % %     ddR3 = (G*(m1 + m2 + m3 + m4)*((m3*m4)/r34^2 + (m1*m4)/(r12 + r23 + r34)^2 + (m2*m4)/(r23 + r34)^2 - (m4*m5*(m1 + m2 + m3))/(r45^2*(m1 + m2 + m3 + m4)) + (m2*m4*m5)/((r23 + r34 + r45)^2*(m1 + m2 + m3 + m4)) + (m1*m4*m5)/((m1 + m2 + m3 + m4)*(r12 + r23 + r34 + r45)^2) + (m3*m4*m5)/((r34 + r45)^2*(m1 + m2 + m3 + m4))))/(m4*(m1 + m2 + m3));
% % %     ddR4 = (G*((m4*m5)/r45^2 + (m3*m5)/(r45 - ((m1 + m2)*(r23 + (m1*r12)/(m1 + m2)))/(m1 + m2 + m3) + (m4*(r34 + ((m1 + m2)*(r23 + (m1*r12)/(m1 + m2)))/(m1 + m2 + m3)))/(m1 + m2 + m3 + m4) + ((r34 + ((m1 + m2)*(r23 + (m1*r12)/(m1 + m2)))/(m1 + m2 + m3))*(m1 + m2 + m3))/(m1 + m2 + m3 + m4))^2 + (m1*m5)/(r45 + (m4*(r34 + ((m1 + m2)*(r23 + (m1*r12)/(m1 + m2)))/(m1 + m2 + m3)))/(m1 + m2 + m3 + m4) + (m3*(r23 + (m1*r12)/(m1 + m2)))/(m1 + m2 + m3) + ((r34 + ((m1 + m2)*(r23 + (m1*r12)/(m1 + m2)))/(m1 + m2 + m3))*(m1 + m2 + m3))/(m1 + m2 + m3 + m4) + (m2*r12)/(m1 + m2))^2 + (m2*m5)/(r45 + (m4*(r34 + ((m1 + m2)*(r23 + (m1*r12)/(m1 + m2)))/(m1 + m2 + m3)))/(m1 + m2 + m3 + m4) + (m3*(r23 + (m1*r12)/(m1 + m2)))/(m1 + m2 + m3) + ((r34 + ((m1 + m2)*(r23 + (m1*r12)/(m1 + m2)))/(m1 + m2 + m3))*(m1 + m2 + m3))/(m1 + m2 + m3 + m4) - (m1*r12)/(m1 + m2))^2)*(m1 + m2 + m3 + m4 + m5))/(m5*(m1 + m2 + m3 + m4));
    ddR1 = -(G*(m1 + m2)*((m1*m2)*r12/norm(r12)^3 - (m1*m2*m3)*r23/(norm(r23)^3*(m1 + m2)) + (m1*m2*m4)*(r12 + r23 + r34)/((m1 + m2)*norm(r12 + r23 + r34)^3) - (m1*m2*m5)*(r23 + r34 + r45)/((m1 + m2)*norm(r23 + r34 + r45)^3) + (m1*m2*m5)*(r12 + r23 + r34 + r45)/((m1 + m2)*norm(r12 + r23 + r34 + r45)^3) + (m1*m2*m3)*(r12 + r23)/((m1 + m2)*norm(r12 + r23)^3) - (m1*m2*m4)*(r23 + r34)/((m1 + m2)*norm(r23 + r34)^3)))/(m1*m2);
    ddR2 = -(G*(m1 + m2 + m3)*((m2*m3)*r23/norm(r23)^3 + (m1*m3)*(r12 + r23)/norm(r12 + r23)^3 - (m3*m5*(m1 + m2))*(r34 + r45)/(norm(r34 + r45)^3*(m1 + m2 + m3)) + (m1*m3*m5)*(r12 + r23 + r34 + r45)/((m1 + m2 + m3)*norm(r12 + r23 + r34 + r45)^3) - (m3*m4*(m1 + m2))*r34/(norm(r34)^3*(m1 + m2 + m3)) + (m2*m3*m4)*(r23 + r34)/(norm(r23 + r34)^3*(m1 + m2 + m3)) + (m1*m3*m4)*(r12 + r23 + r34)/((m1 + m2 + m3)*norm(r12 + r23 + r34)^3) + (m2*m3*m5)*(r23 + r34 + r45)/((m1 + m2 + m3)*norm(r23 + r34 + r45)^3)))/(m3*(m1 + m2));
    ddR3 = -(G*(m1 + m2 + m3 + m4)*((m3*m4)*r34/norm(r34)^3 + (m1*m4)*(r12 + r23 + r34)/norm(r12 + r23 + r34)^3 + (m2*m4)*(r23 + r34)/norm(r23 + r34)^3 - (m4*m5*(m1 + m2 + m3))*r45/(norm(r45)^3*(m1 + m2 + m3 + m4)) + (m2*m4*m5)*(r23 + r34 + r45)/(norm(r23 + r34 + r45)^3*(m1 + m2 + m3 + m4)) + (m1*m4*m5)*(r12 + r23 + r34 + r45)/((m1 + m2 + m3 + m4)*norm(r12 + r23 + r34 + r45)^3) + (m3*m4*m5)*(r34 + r45)/(norm(r34 + r45)^3*(m1 + m2 + m3 + m4))))/(m4*(m1 + m2 + m3));
    val1 = r45 - ((m1 + m2)*(r23 + (m1*r12)/(m1 + m2)))/(m1 + m2 + m3) + (m4*(r34 + ((m1 + m2)*(r23 + (m1*r12)/(m1 + m2)))/(m1 + m2 + m3)))/(m1 + m2 + m3 + m4) + ((r34 + ((m1 + m2)*(r23 + (m1*r12)/(m1 + m2)))/(m1 + m2 + m3))*(m1 + m2 + m3))/(m1 + m2 + m3 + m4);
    val2 = r45 + (m4*(r34 + ((m1 + m2)*(r23 + (m1*r12)/(m1 + m2)))/(m1 + m2 + m3)))/(m1 + m2 + m3 + m4) + (m3*(r23 + (m1*r12)/(m1 + m2)))/(m1 + m2 + m3) + ((r34 + ((m1 + m2)*(r23 + (m1*r12)/(m1 + m2)))/(m1 + m2 + m3))*(m1 + m2 + m3))/(m1 + m2 + m3 + m4) + (m2*r12)/(m1 + m2);
    val3 = r45 + (m4*(r34 + ((m1 + m2)*(r23 + (m1*r12)/(m1 + m2)))/(m1 + m2 + m3)))/(m1 + m2 + m3 + m4) + (m3*(r23 + (m1*r12)/(m1 + m2)))/(m1 + m2 + m3) + ((r34 + ((m1 + m2)*(r23 + (m1*r12)/(m1 + m2)))/(m1 + m2 + m3))*(m1 + m2 + m3))/(m1 + m2 + m3 + m4) - (m1*r12)/(m1 + m2);
    ddR4 = -(G*((m4*m5)*r45/norm(r45)^3 + (m3*m5)*val1/norm(val1)^3 + (m1*m5)*val2/norm(val2)^3 + (m2*m5)*val3/norm(val3)^3)*(m1 + m2 + m3 + m4 + m5))/(m5*(m1 + m2 + m3 + m4));
    
%     %%% Test
%     r12_mag3 = norm(r12)^3;
%     r23_mag3 = norm(r23)^3;
%     r34_mag3 = norm(r34)^3;
%     r45_mag3 = norm(r45)^3;
%     r15_mag3 = norm(r15)^3;
%     r13 = r12 + r23; r13_mag3 = norm(r13)^3;
%     r14 = r13 + r34; r14_mag3 = norm(r14)^3;
%     r15 = r14 + r45; r15_mag3 = norm(r15)^3;
%     r24 = r23 + r34; r24_mag3 = norm(r24)^3;
%     r25 = r24 + r45; r25_mag3 = norm(r25)^3;
%     r35 = r34 + r45; r35_mag3 = norm(r35)^3;
%     
%     ddR1 = (-G*m1*m2/r12_mag3*r12 - G*m1*m2/(m1 + m2)*(m3/r13_mag3*r13 + m4/r14_mag3*r14 + m5/r15_mag3*r15 - m3/r23_mag3*r23 - m4/r24_mag3*r24 - m5/r25_mag3*r25))*((m1 + m2)/m1/m2);
%     ddR2 = (-G*m1*m3/r13_mag3*r13 - G*m2*m3/r23_mag3*r23 - G*m3/(m1 + m2 + m3)*(m1*m4/r14_mag3*r14 + m1*m5/r15_mag3*r15 + m2*m4/r24_mag3*r24 + m2*m5/r25_mag3*r25 - (m1 + m2)*m4/r34_mag3*r34 - (m1 + m2)*m5/r35_mag3*r35))*((m1 + m2 + m3)/m3/(m1 + m2));
%     ddR3 = (-G*m1*m4/r14_mag3*r14 - G*m2*m4/r24_mag3*r24 - G*m3*m4/r34_mag3*r34 - G*m4*m5/(m1 + m2 + m3 + m4)*(m1/r15_mag3*r15 + m2/r25_mag3*r25 + m3/r35_mag3*r35 - (m1 + m2 + m3)/r45_mag3*r45))*((m1 + m2 + m3 + m4)/m4/(m1 + m2 + m3));
%     ddR4 = -G*m5*(m1/r15_mag3*r15 + m2/r25_mag3*r25 + m3/r35_mag3*r35 + m4/r45_mag3*r45)*((m1 + m2 + m3 + m4 + m5)/m5/(m1 + m2 + m3 + m4));

    
    %%% Storing EOMs
    ddX = [ddR1; ddR2; ddR3; ddR4];
    
    %%% Output the derivative of the state
    dX = zeros(24,1);
%     warning('check these dimensions ... should it be ; or , ?')
    dX(1:12) = [R1dot; R2dot; R3dot; R4dot];
    dX(13:24) = ddX;
    
end

function [E] = jacobiCoord_5Bp_Energy(R1,R2,R3,R4,R1dot,R2dot,R3dot,R4dot,m,G)
m1 = m(1); m2 = m(2); m3 = m(3); m4 = m(4); m5 = m(5);
E = zeros(length(R1),1);
for kk = 1:size(R1,1)
    %%% Calculating kinetic energy
    T = dot(R1dot(kk,:),R1dot(kk,:))*m1*m2/(2*(m1+m2)) + dot(R2dot(kk,:),R2dot(kk,:))*m3*(m1+m2)/(2*(m1+m2+m3)) + dot(R3dot(kk,:),R3dot(kk,:))*m4*(m1+m2+m3)/(2*(m1+m2+m3+m4)) + dot(R4dot(kk,:),R4dot(kk,:))*m5*(m1+m2+m3+m4)/(2*(m1+m2+m3+m4+m5));

    %%% Converting to relative coordinates
    r12 = R1;
    r13 = R2 + m2*R1/(m1+m2);
    r14 = R3 + m3*R2/(m1+m2+m3) + m2*R1/(m1+m2);
    r15 = R4 + m4*R3/(m1+m2+m3+m4) + m3*R2/(m1+m2+m3) + m2*R1/(m1+m2);
    r23 = R2 - m1*R1/(m1+m2);
    r24 = R3 + m3*R2/(m1+m2+m3) - m1*R1/(m1+m2);
    r25 = R4 + m4*R3/(m1+m2+m3+m4) + m3*R2/(m1+m2+m3) - m1*R1/(m1+m2);
    r34 = R3 - (m1+m2)*R2/(m1+m2+m3);
    r35 = R4 + m4*R3/(m1+m2+m3+m4) - (m1+m2)*R2/(m1+m2+m3);
    r45 = R4 - (m1+m2+m3)*R3/(m1+m2+m3+m4);

    % Calculating Potential Energy
%     U = -G*(m1*m2/norm(R1(kk,:)) + m1*m3/norm(R2(kk,:) + m2*R1(kk,:)/(m1+m2)) + m1*m4/norm(R3(kk,:) + m3*R2(kk,:)/(m1+m2+m3) + m2*R1(kk,:)/(m1+m2)) + m1*m5/norm(R4(kk,:) + m4*R3(kk,:)/(m1+m2+m3+m4) + m3*R2(kk,:)/(m1+m2+m3) + m2*R1(kk,:)/(m1+m2)) + m2*m3/norm(R2(kk,:) - m1*R1(kk,:)/(m1+m2)) + m2*m4/norm(R3(kk,:) + m3*R2(kk,:)/(m1+m2+m3) - m1*R1(kk,:)/(m1+m2)) + m2*m5/norm(R4(kk,:) + m4*R3(kk,:)/(m1+m2+m3+m4) + m3*R2(kk,:)/(m1+m2+m3) - m1*R1(kk,:)/(m1+m2)) + m3*m4/norm(R3(kk,:) - (m1+m2)*R2(kk,:)/(m1+m2+m3)) + m3*m5/norm(R4(kk,:) + m4*R3(kk,:)/(m1+m2+m3+m4) - (m1+m2)*R2(kk,:)/(m1+m2+m3)) + m4*m5/norm(R4(kk,:) - (m1+m2+m3)*R3(kk,:)/(m1+m2+m3+m4)));
    U = -G*(m1*m2/norm(r12(kk,:)) + m1*m3/norm(r13(kk,:)) + m1*m4/norm(r14(kk,:)) + m1*m5/norm(r15(kk,:)) + m2*m3/norm(r23(kk,:)) + m2*m4/norm(r24(kk,:)) + m2*m5/norm(r25(kk,:)) + m3*m4/norm(r34(kk,:)) + m3*m5/norm(r35(kk,:)) + m4*m5/norm(r45(kk,:)));
    
    % Total energy
    E(kk) = T + U;

end
end

function [H] = jacobiCoord_5Bp_AngMom(R1,R2,R3,R4,R1dot,R2dot,R3dot,R4dot,m)
m1 = m(1); m2 = m(2); m3 = m(3); m4 = m(4); m5 = m(5);
H = zeros(length(R1),3);
for kk = 1:size(R1,1)
    H(kk,:) = cross(R1(kk,:),R1dot(kk,:))*m1*m2/(m1+m2) + cross(R2(kk,:),R2dot(kk,:))*m3*(m1+m2)/(m1+m2+m3) + cross(R3(kk,:),R3dot(kk,:))*m4*(m1+m2+m3)/(m1+m2+m3+m4) + cross(R4(kk,:),R4dot(kk,:))*m5*(m1+m2+m3+m4)/(m1+m2+m3+m4+m5);
end
end















