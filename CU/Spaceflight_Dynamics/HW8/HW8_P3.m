clear
clc
close all
addpath('../../bin')
% ------------------------------------------------------------------------
%%% Given
% ------------------------------------------------------------------------
Re = 6378.1363; % km
uE = 398600.4415; % km^3/s^2
rp = 400 + Re; % km
ra = 1000 + Re; % km
% ra = 40000 + Re; % km
i = 51.5*pi/180; % rad
raan = 0; % rad
w = 60*pi/180; % rad
ta = 0; % rad
Force = 1E-9; % km/s^2

% ------------------------------------------------------------------------
%%% Preparing to Propagate
% ------------------------------------------------------------------------
a = (rp+ra)/2; % km
e = (ra-rp)/(ra+rp); 
[rECIi, vECIi] = OE2ECI(a, e, i, raan, w, ta, uE);
% ------------------------------------------------------------------------
%%% Propogating Orbit
% ------------------------------------------------------------------------
ti = 0; % sec
tf = 3600*4; % sec
t = [ti:1:tf]; % sec
%%% Setting integrator accuracy
tol = 1E-12;
options = odeset('RelTol',tol,'AbsTol',tol);

%%% Numerically propagating Cygnus orbit
[Times,States] = ode45(@orbitPropagator,t,[rECIi;vECIi],options,uE);

% plot3(States(:,1),States(:,2),States(:,3))


% ------------------------------------------------------------------------
%%% Analyzing Orbital Elements Throughout Propagation
% ------------------------------------------------------------------------
%%% Rate of change equations
didt = zeros(size(States,1),1);
dOmdt = zeros(size(States,1),1);
dwdt = zeros(size(States,1),1);
ta_prop = zeros(size(States,1),1);

%%% Actual change
a_prop = zeros(size(States,1),1);
e_prop = zeros(size(States,1),1);
i_prop = zeros(size(States,1),1);
raan_prop = zeros(size(States,1),1);
w_prop = zeros(size(States,1),1);

for k = 1:size(States,1)
    %%% Determining OE at state
    [a,e,i,raan,w,ta] = ECI2OE(States(k,1:3),States(k,4:6),uE);
    %%% Rate vs true anomaly
    r = norm(States(k,1:3));
    n = sqrt(uE/(a^3));
    h = norm(cross(States(k,1:3),States(k,4:6)));
    didt(k) = r*cos(w+ta)*Force/(n*(a^2)*sqrt(1-e^2))*(180/pi)*(86400); % deg/day
    dOmdt(k) = r*sin(w+ta)*Force/(n*(a^2)*sqrt(1-e^2)*sin(i))*(180/pi)*(86400); % deg/day
    dwdt(k) = (-r*cot(i)*sin(w+ta)*Force/h)*(180/pi)*(86400); % deg/day
    
    
    a_prop(k) = a;
    e_prop(k) = e;
    i_prop(k) = i;
    raan_prop(k) = raan;
    if abs(raan_prop(k)-2*pi) < 0.001
        raan_prop(k) = 0;
    end
    w_prop(k) = w;
    ta_prop(k) = ta;
end

%%% Functions of True Anomaly
figure
subplot(3,1,1)
plot(ta_prop*180/pi,didt,'.');
PlotBoi2('','di/dt, [°/day]',18); xlim([0 360])
subplot(3,1,2)
plot(ta_prop*180/pi,dOmdt,'.');
PlotBoi2('','d\Omega/dt, [°/day]',18); xlim([0 360])
subplot(3,1,3)
plot(ta_prop*180/pi,dwdt,'.');
PlotBoi2('True Anomaly, [°]','d\omega/dt, [°/day]',18); xlim([0 360])


% %%% Functions of Time
% figure
% subplot(3,1,1)
% plot((ti:tf)/3600,didt,'.');
% PlotBoi2('','di/dt, [°/day]',18)
% subplot(3,1,2)
% plot((ti:tf)/3600,dOmdt,'.');
% PlotBoi2('','d\Omega/dt, [°/day]',18)
% subplot(3,1,3)
% plot((ti:tf)/3600,dwdt,'.');
% PlotBoi2('Time, [hour]','d\omega/dt, [°/day]',18)

% %%% OEs over time
% figure
% subplot(3,2,1)
% plot((ti:tf)/3600,ta_prop*180/pi,'.');
% PlotBoi2('','True Anomaly [°]',16);
% subplot(3,2,2)
% plot((ti:tf)/3600,a_prop);
% PlotBoi2('','Semimajor Axis [km]',16);
% subplot(3,2,3)
% plot((ti:tf)/3600,e_prop);
% PlotBoi2('','Eccentricity',16);
% subplot(3,2,4)
% plot((ti:tf)/3600,i_prop*180/pi,'.');
% PlotBoi2('','Inclination [°]',16);
% subplot(3,2,5)
% plot((ti:tf)/3600,raan_prop*180/pi,'.');
% PlotBoi2('Time [hour]','\Omega, [°]',16);
% subplot(3,2,6)
% plot((ti:tf)/3600,w_prop*180/pi,'.');
% PlotBoi2('Time [hour]','\omega, [°]',16);

% %%% Isolating Semimajor Axis and Eccentricity
% figure
% subplot(2,1,1)
% plot((ti:tf)/3600,a_prop);
% PlotBoi2('','Semimajor Axis [km]',16);
% subplot(2,1,2)
% plot((ti:tf)/3600,e_prop);
% PlotBoi2('Time [hour]','Eccentricity',16);
% ------------------------------------------------------------------------
%%% Plotting Book Equations
% ------------------------------------------------------------------------



% ------------------------------------------------------------------------
%%% Numerical Integration Function
% ------------------------------------------------------------------------
function [ dY ] = orbitPropagator(t,Y,uE)
%%% Unpack the Cygnus state vector
r = Y(1:3); % km
v = Y(4:6); % km/s

%%% Direction and magnitude of unknown forcing
h = cross(r,v);
Fdir = h/norm(h);
Force = 1E-9; % km/s^2
%%% Dynamics of the system (acceleration cause by Earth)
v_dot = (-uE/(norm(r)^3))*r + Fdir*Force; % km/s^2

%%% Initializing output state vector
dY = zeros(6,1);

%%% Output the derivative of the state
dY(1:3) = v;
dY(4:6) = v_dot;
end





