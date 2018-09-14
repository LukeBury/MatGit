clear
clc
close all
addpath('../bin')

% ------------------------------------------------------------------------
%%% Givens
% ------------------------------------------------------------------------
%%% Orbital Elements
a = 7000; % km
e = 0.001; 
i = 30*pi/180; % rad
RAAN = 80*pi/180; % rad
w = 40*pi/180; % rad
ta = 0; % rad

%%% State deviation vector
dx = [1; 0; 0; 0; .01; 0];%%%%%%%%%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ 10 m/s? .01 km/s

%%% Other
uE = 398600.4415; % km^3/s^2
J2 = 0.0010826269;
g = 9.81;

% ------------------------------------------------------------------------
%%% Finding Initial Position and Setting State
% ------------------------------------------------------------------------
%%% Converting OE to ECI state
[r0, v0] = OE2ECI(a, e, i, RAAN, w, ta, uE); % km, km/s

%%% Setting initial state vector (6x1)
X0 = [r0; v0];

% ------------------------------------------------------------------------
%%% Propagating the State with Numerical Integration
% ------------------------------------------------------------------------
%%% Setting time frame
ti = 0; % sec
tf = 24*3600; % sec
dt = 1; % sec
time = ti:dt:tf; % sec

%%% Setting integrator accuracy
tol = 1E-10;
options = odeset('RelTol',tol,'AbsTol',tol);

%%% Propagating the State
[Times,States] = ode45(@StatOD_Hw1_Int,time,X0,options,uE);


plot3(States(:,1),States(:,2),States(:,3))






% ------------------------------------------------------------------------
%%% Numerical Integrator
% ------------------------------------------------------------------------
function [ dY ] = StatOD_Hw1_Int(t,Y,uE)
dY = zeros(6,1);

%%% Unpack the state vector (ECI)
y = Y(1:3); % Hopper Position, km
dy = Y(4:6); % Hopper Velocity, km/s

%%% System Dynamics
aG = (-uE/(norm(y)^3))*y; % km/s^2
% aJ2 = ;
% ddy = aG + aJ2;
ddy = aG;

%%% Output the derivative of the state
dY(1:3) = dy; % km/s
dY(4:6) = ddy; % km/s^2
end




