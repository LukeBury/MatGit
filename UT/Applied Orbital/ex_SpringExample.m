function SpringExample
clc; clear all; close all

% Parameters of the system
m = 1;      % Mass [kg]
c = 0.5;    % Damping coefficient [Ns/m]
k = 1;      % Spring constant [N/m]

% Set up the initial conditions
y0 = -1;    % Initial height [m]
dy0 = 0;    % Initial speed [m/s]
% Put them together in one vector
Y0 = [y0; dy0];

% Set the time range
t0 = 0;     % Initial time [s]
tf = 10;    % Final time [s]

% Set up the options
options = odeset('AbsTol',1E-6,'RelTol',1E-6);
% Can also include other options!

% Call the integrator
[t,Y] = ode45(@springIntegrator,[t0,tf],Y0,options,m,c,k);

% Do something with the output
y = Y(:,1);
dy = Y(:,2);
figure; hold on
plot(t,y,'-b')  % Plot the position over time
plot(t,dy,'-r') % Plot the speed over time
title('Oscillating spring','FontSize',14)
xlabel('Time [s]','FontSize',12)
ylabel('Parameter','FontSize',12)
legend('y [m]','dy [m/s]','Location','Best')








function [dY] = springIntegrator(t,Y,m,c,k)
dY = zeros(2,1);

% Unpack the state vector
y = Y(1);
dy = Y(2);

% Dynamics of the system
% Here's where the fun happens!
ddy = (-c*dy - k*y)/m;

% Output the DERIVATIVE of the state
dY(1) = dy;
dY(2) = ddy;