clear
clc
% close all
mbinPath = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/mbin';
addpath(genpath(mbinPath))

% ========================================================================
%%% Importing Data
% ========================================================================
%%% General data on solar system bodies
bodies = getBodyData(mbinPath);

%%% Color options/schemes
colors = get_colors();

% ========================================================================
%%% Run Switches
% ========================================================================


% ========================================================================
%%% HW3
% ========================================================================
% ------------------------------------
%%% Given
% ------------------------------------
JD_departure = 2454085.5; % days
transferTime_days = 830;  % days
transferTime_sec = transferTime_days*86400; % sec

JD_arrival = JD_departure + transferTime_days;

% ------------------------------------
%%% Departure and arrival positions of Mars and Jupiter
% ------------------------------------
%%% Find Mars position at departure time
[L_rad, a, e, i_rad, Omega_rad, Pi_rad, w_rad, M_rad, ta_rad] = ...
    getPlanetElements_Meeus(JD_departure, 'Mars', 'radians');
[rMars_departure, vMars_departure] = OE2ECI(a, e, i_rad, Omega_rad, w_rad, ta_rad, bodies.sun.u);

%%% Find Jupiter position at arrival time
[L_rad, a, e, i_rad, Omega_rad, Pi_rad, w_rad, M_rad, ta_rad] = ...
    getPlanetElements_Meeus(JD_arrival, 'Jupiter', 'radians');
[rJup_arrival, vJup_arrival] = OE2ECI(a, e, i_rad, Omega_rad, w_rad, ta_rad, bodies.sun.u);

% ----------------------
%%% Computing nominal Lambert solution
% ----------------------
transferType = 1;
[v1,v2] = lambertTargeting_psiTest(rMars_departure, rJup_arrival, transferTime_sec, transferType, 1, bodies.sun.u, 0);
% transferType = 1;
% [v1,v2] = lambertTargeting(rMars_departure, rJup_arrival, transferTime_sec, transferType, 1, bodies.sun.u, 0);

% -------------------------------------------------
%
% -------------------------------------------------
%%% Creating time vector
t0 = 0;             % sec
tf = transferTime_sec; % sec
time = linspace(t0,tf,10000);
% tf_sc = transferTime_sec;
% time_sc = linspace(t0,tf_sc,10000);

%%% Choosing ode tolerance
tol = 2.22045e-14;

%%% Setting integrator options
options = odeset('RelTol',tol,'AbsTol',tol);


[time, X_Mars] = ode113(@Int_2BI,time, [rMars_departure; vMars_departure], options, bodies.sun.u);
[time, X_Jup] = ode113(@Int_2BI,time, [rJup_arrival; vJup_arrival], options, bodies.sun.u);

[time, X_sc] = ode113(@Int_2BI,time, [rMars_departure; v1], options, bodies.sun.u);

figure; hold all
plot3(X_Mars(:,1),X_Mars(:,2),X_Mars(:,3),'r','linewidth',2)
plot3(X_Mars(1,1),X_Mars(1,2),X_Mars(1,3),'ro','linewidth',2)
plot3(X_Jup(:,1),X_Jup(:,2),X_Jup(:,3),'g','linewidth',2)
plot3(X_Jup(1,1),X_Jup(1,2),X_Jup(1,3),'go','linewidth',2)
plot3(X_sc(:,1),X_sc(:,2),X_sc(:,3),'c','linewidth',2)
PlotBoi3('x','y','z',15)

% missDist = norm(X_sc(end,1:3)' - rJup_arrival)





% % ----------------------
% %%% Looping through elliptical TOFs
% % ----------------------
% %%% Setting TOFs
% % tofs_elliptical_days = linspace(150,1200,1000);
% tofs_elliptical_days = linspace(400,1200,1000);
% tofs_elliptical_sec = tofs_elliptical_days*86400;
% 
% %%% Lambert options
% nOrbits = 1;
% 
% %%% Preallocating
% psis_elliptical_rad = zeros(length(tofs_elliptical_days),1);
% 
% %%% Looping through TOFs
% for kk = 1:length(tofs_elliptical_days)
%     %%% Compute nominal Lambert solution
%     [V1, V2,exitflag,psi] = lambertSolver_psi(rMars_departure, rJup_arrival, tofs_elliptical_sec(kk), nOrbits, 0, bodies.sun.u);
%     
%     psis_elliptical_rad(kk) = psi;
% end
% 
% 
% figure; hold all
% plot(psis_elliptical_rad, tofs_elliptical_days,'linewidth',2)
% PlotBoi2('$\Psi$, $rad^2$','TOF,$days$',18,'LaTex')
















