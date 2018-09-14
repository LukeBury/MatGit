clear
clc
close all
addpath(genpath('/Users/lukebury/Documents/MATLAB/mbin'))

% ------------------------------------------------------------------------
%%% P1
% ------------------------------------------------------------------------
% syms I1 I2 I3 I12 I13 I23 wd1 wd2 wd3 w1 w2 w3 dwd1 dwd2 dwd3 dw1 dw2 dw3 L1 L2 L3 real
% syms wr1 wr2 wr3 wdr1 wdr2 wdr3 real
% 
% %%% (a)
% I = [I1 I12 I13; I12 I2 I23; I13 I23 I3];
% wd = [wd1;wd2;wd3];
% wt = [0 -w3 w2; w3 0 -w1; -w2 w1 0];
% w = [w1;w2;w3];
% L = [L1;L2;L3];
% 
% I*wd
% 
% -wt*I*w + L
% 
% %%% (b)
% wd = [wdr1 + dwd1; wdr2 + dwd2; wdr3 + dwd3];
% w = [wr1 + dw1; wr2 + dw2; wr3 + dw3];
% 
% wt = [0 -w(3) w(2); w(3) 0 -w(1); -w(2) w(1) 0];
% 
% eq = -wt*I*w + L

% ------------------------------------------------------------------------
%%% P3
% ------------------------------------------------------------------------
% % syms wn k p11 p12 p21 p22 R1 R2 real
% % A = [0 1; -wn^2 -2*k*wn];
% % P = [p11 p12; p21 p22];
% % 
% % R = A'*P + P*A
% % 
% % clear P
% % syms R11 R12 R21 R22 real
% % 
% % % P = zeros(2,2);
% % % P(1,1) = R12 - (R22+R11/wn^2)*wn/(4*k)-R11*k/wn;
% % % P(1,2) = -R11/(2*wn^2);
% % % P(2,1) = P(1,2);
% % % P(2,2) = -(R22 + R11/wn^2)/(4*k*wn);
% % P = [(R22+R11/wn^2)*wn/(4*k)+R11*k/wn, R11/(2*wn^2);...
% %     R11/(2*wn^2), (R22 + R11/wn^2)/(4*k*wn)];
% % assume(R22 > 0); assume(R11 > 0); assume(wn > 0); assume(1>k>0);
% % % assume(wn > 1)
% % [v,e] = eig(P)
% % isAlways(e(1,1)>0)
% % % syms x1 x2 real
% % % x = [x1;x2];
% % % assume(x1>0); assume(x2>0);
% % % eqn = x'*P*x;
% % % isAlways(eqn(1,1)>0)
% ------------------------------------------------------------------------
%%% P6 - 8.8
% ------------------------------------------------------------------------
clear; clc
m = 1; % kg
c = 0.1; % kg/s
k1 = 1; % kg/s^2
k3 = 0.1; % kg*s^-2*m^-2

x0 = 2; % m
dx0 = 0; % m/s
X0 = [x0; dx0];

tol = 1e-9;
options = odeset('RelTol',tol,'AbsTol',tol);

%%% Setting time vector and normalizing 
ti = 0; % sec
dt = .1; % sec
tf = 40; % sec
times = [ti:dt:tf];

%%% Creating Plot titles
titles = {'x_r(t) = 0','x_r(t) = sin(0.5t)'};

for problem = 1:2
    %%% Create reference
    if problem == 1
        xRef = zeros(length(times),1);
    elseif problem == 2
        xRef = sin(0.5*times);
    end
    
    %%% Propagating the States without J2
    [times, states] = ode45(@Int_5010hw5P6, times, X0, options,m,c,k1,k3,problem);
    
    %%% Plottin Results
    figure
    subplot(2,1,1); hold all
    plot(times,states(:,1),'linewidth',1.5)
    plot(times,xRef,'linewidth',1.5)
    PlotBoi2('','x, m',12)
    title(titles{problem})
    legend('State','Reference')
    
    subplot(2,1,2)
    plot(times,states(:,1)-xRef,'linewidth',1.5);
    title(sprintf('Final dx value = %1.2e m',states(end,1)-xRef(end)))
    PlotBoi2('Time, s','dx, m',12)
end

function [dX] = Int_5010hw5P6(t,X,m,c,k1,k3,problem)
%%% Preallocate state output
dX = zeros(2,1);

%%% Create reference
if problem == 1
    xr = 0;
    xrd = 0;
    xrdd = 0;
elseif problem == 2
    xr = sin(0.5*t);
    xrd = 0.5*cos(0.5*t);
    xrdd = -0.25*sin(0.5*t);
end

%%% Unpack the state vector
x = X(1);
xd = X(2);

dx = x-xr;
dxd = xd - xrd;

%%% create control
u = -dxd + m*xrdd + c*(xrd + dxd) + k1*xr + k3*(3*dx^2*xr + 3*dx*xr^2 + xr^3);

%%% Equations of Motion
xdd = (u - c*xd - k1*x - k3*(x^3))/m;

%%% Output the derivative of the state
dX(1) = xd;
dX(2) = xdd;
end

%%%%%% other symbolic work
syms xrdd c m xrd dxd k1 k3 dx xr real
u = -dxd + m*xrdd + c*(xrd + dxd) + k1*xr + k3*(3*dx^2*xr + 3*dx*xr^2 + xr^3);

dxdd = -xrdd + u/m - c*(xrd+dxd)/m - k1*(xr+dx)/m - k3*((xr+dx)^3)/m;

vddd = -2*dxdd*dxdd



