clc

% ------------------------------------------------------------------------
%%% Givens
% ------------------------------------------------------------------------
r_ECIi = [5492; 3984.001; 2.955]; % km
v_ECIi = [-3.931; 5.498; 3.665]; % km
uE = 398600.4415; % km^3/s^2

global t2_100 r2_100
% ------------------------------------------------------------------------
%%% Preparing for and Calling the Numerical Integrator
% ------------------------------------------------------------------------
ti = 0; % sec
tf = 1000000; % sec
t = [ti:1:tf]; % sec
%%% Setting integrator accuracy
tol = 1E-06;
options = odeset('RelTol',tol,'AbsTol',tol);

%%% Numerically propagating Cygnus orbit
[Time,State] = ode45(@CygnusPropagator,t,[r_ECIi;v_ECIi],options,uE);

% ------------------------------------------------------------------------
%%% Printing Results
% ------------------------------------------------------------------------
%%% 100 sec results
fprintf('At t = %f sec:\n',t2_100)
r2_100 = State(101,1:3)';
table(r1_100, r2_100)

%%% 1,000,000 sec results
fprintf('At t = 1,000,000 sec:\n')
r2_1000000 = State(end,1:3)'; % km
table(r1_1000000, r2_1000000) % km

% ------------------------------------------------------------------------
%%% Calculating & Printing Difference in Position
% ------------------------------------------------------------------------
if r1_1000000
    fprintf('tolerance = %1.1d\n', tol)
    diff = sqrt((r1_1000000(1)-r2_1000000(1))^2 + (r1_1000000(2)-...
        r2_1000000(2))^2 + (r1_1000000(3)-r2_1000000(3))^2); % km
    fprintf('t_1000000 difference = %1.4d km\n',diff)
end

% ------------------------------------------------------------------------
%%% Numerical Integration Function
% ------------------------------------------------------------------------
function [ dY ] = CygnusPropagator(t,Y,uE)
%%% Bring in global variable for state(100 sec)
global t2_100 r2_100

%%% Unpack the Cygnus state vector
r = Y(1:3); % km
v = Y(4:6); % km/s

%%% Catching the 100 sec state
if t > 99.5 && t < 100.05
    t2_100 = t; 
    r2_100 = r;
end

%%% Dynamics of the system (acceleration cause by Earth)
v_dot = (-uE/(norm(r)^3))*r; % m/s^2

%%% Initializing output state vector
dY = zeros(6,1);

%%% Output the derivative of the state
dY(1:3) = v;
dY(4:6) = v_dot;
end