clear
clc
close all
addpath('/Users/lukebury/Documents/MATLAB/mbin/plot')

% Run switch
run_sym = 0;
run_3a = 0;
run_3b = 1;
run_3c = 0;

% ========================================================================
%%% Problem 3
% ========================================================================

% -------------------------------------------------------------
%%% Partials Calculations
% -------------------------------------------------------------
if run_sym == 1
   syms m r phi real
   V = -(m/r)*(1 + (1/(20*r*r))*(1 + 3*cos(2*phi)));
   dVdr = diff(V,r);
   dVdphi = diff(V,phi);
   
   dVdr = simplify(dVdr);
   dVdphi = simplify(dVdphi);
   
%    dVdr = (m*(20*r^2 + 9*cos(2*phi) + 3))/(20*r^4);
%    dVdphi = (3*m*sin(2*phi))/(10*r^3);

    
end


% -------------------------------------------------------------
%%% 3a
% -------------------------------------------------------------
if run_3a == 1
%%% Integrator options
tol = 1e-7;
options = odeset('RelTol',tol,'AbsTol',tol);

% Setting timing options
t0 = 0;
dt = .1;
tf = 50;
tVec = t0:dt:tf;

%%% Setting constants
Iz = 0.2;
v = 0.75;
m = v*(1-v);

%%% Setting initial state conditions
r_0        = 15;
phi_0      = .1;
theta_0    = pi/5;
rdot_0     = .001;
phidot_0   = .02;
thetadot_0 = .004;

%%% Setting intial state vector
X0 = [r_0 phi_0 theta_0 rdot_0 phidot_0 thetadot_0];

%%% Propagating the state
[times,states] = ode45(@CH3_EOM_Int,tVec,X0,options,Iz,m);

% %%% Unpacking results
% r = states(:,1);
% phi = states(:,2);
% theta = states(:,3);
% rdot = states(:,4);
% phidot = states(:,5);
% thetadot = states(:,6);

%%% Calculating Energy and Angular Momentum
[E, H] = CH3_EH_comp(states,Iz,m);

%%% Plotting Energy and AngMom results
figure
subplot(2,1,1); plot(times,E-E(1),'linewidth',1.5)
PlotBoi2('','\DeltaE',14)
subplot(2,1,2); plot(times,H-H(1),'linewidth',1.5)
PlotBoi2('Time','\DeltaH',14)

end


% -------------------------------------------------------------
%%% 3b
% -------------------------------------------------------------
if run_3b == 1
%%% Integrator options
tol = 1e-7;
options = odeset('RelTol',tol,'AbsTol',tol);

% Setting timing options
t0 = 0;
dt = .1;
tf = 180;
tVec = t0:dt:tf;

for trial = 1:5
    if trial == 1
        v = 0.75;
        r_0 = 2;
        phi_0 = 0;
    elseif trial == 2
        v = 0.5;
        r_0 = 1;
        phi_0 = 1*pi/180;
    elseif trial == 3
        v = 0.5;
        r_0 = 1.5;
        phi_0 = 1*pi/180;
    elseif trial == 4
        v = 0.5;
        r_0 = 2;
        phi_0 = 1*pi/180;
    elseif trial == 5
        v = 0.5;
        r_0 = 3;
        phi_0 = 1*pi/180;
    end
%%% Setting constants
Iz = 0.2;
m = v*(1-v);

%%% Setting initial state conditions
theta_0    = pi/5;
rdot_0     = 0;
phidot_0   = 0;
dVdr_phi0 = (m*(20*r_0^2 + 9*cos(2*0) + 3))/(20*r_0^4);
thetadot_0 = sqrt(dVdr_phi0/(m*r_0));

%%% Setting intial state vector
X0 = [r_0 phi_0 theta_0 rdot_0 phidot_0 thetadot_0];

%%% Propagating the state
[times,states] = ode45(@CH3_EOM_Int,tVec,X0,options,Iz,m);

% %%% Unpacking results
r = states(:,1);
phi = states(:,2);
theta = states(:,3);
rdot = states(:,4);
phidot = states(:,5);
thetadot = states(:,6);

%%% Plotting state results
figure
subplot(3,1,1); plot(times,r,'linewidth',1.5)
title(strcat('\nu',sprintf(' = %0.2f  |  r_0 = %0.1f',v,r_0),'  |  \phi',sprintf('_0 = %0.0f°',phi_0*180/pi)))
PlotBoi2('','r',14)
subplot(3,1,2); plot(times,phi,'linewidth',1.5)
PlotBoi2('','\phi (rad)',14)
subplot(3,1,3); plot(times,theta,'linewidth',1.5)
PlotBoi2('Time','\theta (rad)',14)

end % 3b trials

end % 3b



% -------------------------------------------------------------
%%% 3c
% -------------------------------------------------------------
if run_3c == 1
%%% Integrator options
tol = 1e-9;
options = odeset('RelTol',tol,'AbsTol',tol);

% Setting timing options
t0 = 0;
dt = 1;
tf = 500;
tVec = t0:dt:tf;

for trial = 1:5
    if trial == 1
        v = 0.001;
    elseif trial == 2
        v = 0.01;
    elseif trial == 3
        v = 0.1;
    elseif trial == 4
        v = 0.25;
    elseif trial == 5
        v = 0.5;
    end
%%% Setting constants
Iz = 0.2;
m = v*(1-v);

%%% Setting initial state conditions
r_0        = 1;
phi_0      = 1*pi/180;
theta_0    = pi/5;
rdot_0     = 0;
phidot_0   = 0;
dVdr_phi0 = (m*(20*r_0^2 + 9*cos(2*0) + 3))/(20*r_0^4);
thetadot_0 = sqrt(dVdr_phi0/(m*r_0));

%%% Setting intial state vector
X0 = [r_0 phi_0 theta_0 rdot_0 phidot_0 thetadot_0];

%%% Propagating the state
[times,states] = ode45(@CH3_EOM_Int,tVec,X0,options,Iz,m);

% %%% Unpacking results
r = states(:,1);
phi = states(:,2);
theta = states(:,3);
rdot = states(:,4);
phidot = states(:,5);
thetadot = states(:,6);

%%% Plotting state results
figure
subplot(3,1,1); plot(times,r,'linewidth',1.5)
title(strcat('\nu',sprintf(' = %0.3f  |  r_0 = %0.1f',v,r_0),'  |  \phi',sprintf('_0 = %0.0f°',phi_0*180/pi)))
PlotBoi2('','r',14)
subplot(3,1,2); plot(times,phi,'linewidth',1.5)
PlotBoi2('','\phi (rad)',14)
subplot(3,1,3); plot(times,theta,'linewidth',1.5)
PlotBoi2('Time','\theta (rad)',14)

end % 3c trials

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


function [dX] = CH3_EOM_Int(t,X,Iz,m)
    %%% Unpack the state vector
    r        = X(1);
    phi      = X(2);
    theta    = X(3);
    rdot     = X(4);
    phidot   = X(5);
    thetadot = X(6);

    %%% Partials of Potential
    dVdr   = (m*(20*r^2 + 9*cos(2*phi) + 3))/(20*r^4);
    dVdphi = (3*m*sin(2*phi))/(10*r^3);
    
    %%% Equations of Motion
    ddr     = thetadot*thetadot*r - dVdr/m;
    ddphi   = -(Iz + m*r*r)*dVdphi/(Iz*m*r*r) + 2*rdot*thetadot/r;
    ddtheta = dVdphi/(m*r*r) - 2*rdot*thetadot/r;
    
    %%% Storing EOMs
    ddX = [ddr; ddphi; ddtheta];
    
    %%% Output the derivative of the state
    dX = zeros(6,1);
    dX(1:3) = [rdot; phidot; thetadot];
    dX(4:6) = ddX;
    
end


function [E, H] = CH3_EH_comp(X,Iz,m)
    E = zeros(size(X,1),1);
    H = zeros(size(X,1),1);

    for kk = 1:size(X,1)
        %%% Unpack the state vector
        r        = X(kk,1);
        phi      = X(kk,2);
        theta    = X(kk,3);
        rdot     = X(kk,4);
        phidot   = X(kk,5);
        thetadot = X(kk,6);
        
        %%% Calculate Potential Energy
        V = -m*(1 + (1 + 3*cos(2*phi))/(20*r*r))/r;
        
        %%% Calculate Energy
        E(kk) = .5*Iz*((thetadot+phidot)^2) + .5*m*rdot*rdot + .5*m*r*r*thetadot*thetadot + V;
        
        %%% Calculate Angular Momentum
        H(kk) = (Iz + m*r*r)*thetadot + Iz*phidot;
    end
end














