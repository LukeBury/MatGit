clear
clc
close all
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
run_prob = 3; % 1, 2, 3


% ========================================================================
%%% HW4 - P1
% ========================================================================
if (run_prob == 1) 
% ------------------------------------
%%% Setup
% ------------------------------------
v_scM_HC        = [-10.8559; -35.9372];            % km/s
V_scM_HC        = norm(v_scM_HC);                   % km/s
v_venus_HC     = [-15.1945; -31.7927];            % km/s
V_venus_HC     = norm(v_venus_HC);                % km/s
r_venus_HC     = [-96948447.3751; 46106976.1901]; % km
R_venus_HC     = norm(r_venus_HC);                % km
bodies.sun.u   = 1.32712440018e11;                % km^2 / s^2
bodies.venus.u = 3.257e5;                         % km^2 / s^2
bodies.venus.R = 6052;                            % km

% ------------------------------------
%%% Energy
% ------------------------------------
energy_sc_HC = (V_scM_HC^2)/2 - bodies.sun.u/R_venus_HC;
fprintf('Specific energy wrt sun: %1.4f km^2/s^2\n',energy_sc_HC)

% ------------------------------------
%%% Turning Angle Study
% ------------------------------------
%%% Calculating v-infinity
vInfM = v_scM_HC - v_venus_HC; % vInf Minus
VInf = norm(vInfM);

%%% Setting possible flyby radii
n_rfb = 10000;
rfb_values = linspace(0,200000,n_rfb); % km

%%% Looping through flyby radiis and calculating turning angles
turningAngles_deg = zeros(n_rfb,1); 

for kk = 1:n_rfb
    turningAngles_deg(kk) = (pi  - 2*acos(1/(rfb_values(kk)*(VInf^2)/bodies.venus.u + 1)))*180/pi; % deg
end

%%% Plotting turning angles
figure; hold all
plot(rfb_values, turningAngles_deg,'b','linewidth',2)
PlotBoi2('Flyby periapsis, $km$','Turning Angle, $^\circ$',18,'LaTex')

% ------------------------------------
%%% Calculating energy after flyby
% ------------------------------------
%%% flight path angle at flyby
phi = acos(dot(v_scM_HC,v_venus_HC)/(V_scM_HC*V_venus_HC)); % rad

energies_postFB_leading = zeros(n_rfb,1);
energies_postFB_trailing = zeros(n_rfb,1);
for kk = 1:n_rfb
    
    psi = turningAngles_deg(kk)*pi/180; % rad
    
    beta = acos((V_scM_HC^2 + VInf^2 - V_venus_HC^2)/(2*V_scM_HC*VInf)); % rad
    
    V_scP_HC_leading  = sqrt(VInf^2 + V_venus_HC^2 - 2*VInf*V_venus_HC*cos(pi - beta - phi - psi));
    V_scP_HC_trailing = sqrt(VInf^2 + V_venus_HC^2 - 2*VInf*V_venus_HC*cos(pi - beta - phi + psi));
    
    energies_postFB_leading(kk)  = (V_scP_HC_leading^2)/2 - bodies.sun.u/R_venus_HC;
    energies_postFB_trailing(kk) = (V_scP_HC_trailing^2)/2 - bodies.sun.u/R_venus_HC;
end

%%% Plotting turning angles
figure; hold all
p1 = plot(rfb_values, energies_postFB_leading,'b','linewidth',2);
p2 = plot([bodies.venus.R bodies.venus.R],[-850 -500],'--r','linewidth',1.5);
p3 = plot([0, 200000],[energy_sc_HC energy_sc_HC],'--g','linewidth',1.5);
PlotBoi2('Flyby periapsis, $km$','Heliocentric Specific Energy, $km^2/s^2$',18,'LaTex')
legend([p1, p2, p3],'Post-Flyby Specific Energy','Venus Radius','Pre-Flyby Specific Energy','location','best')
title('Leading')

figure; hold all
p1 = plot(rfb_values, energies_postFB_trailing,'b','linewidth',2);
p2 = plot([bodies.venus.R bodies.venus.R],[-700 -350],'--r','linewidth',1.5);
p3 = plot([0, 200000],[energy_sc_HC energy_sc_HC],'--g','linewidth',1.5);
PlotBoi2('Flyby periapsis, $km$','Heliocentric Specific Energy, $km^2/s^2$',18,'LaTex')
legend([p1, p2, p3],'Post-Flyby Specific Energy','Venus Radius','Pre-Flyby Specific Energy','location','best')
title('Trailing')


end % run_p1

% ========================================================================
%%% HW4 - P2
% ========================================================================
if (run_prob == 2) 
% ------------------------------------
%%% Setup
% ------------------------------------
vInfM = [-5.19425; 5.19424; -5.19425]; % km/s
vInfP = [-8.58481; 1.17067; -2.42304]; % km/s

VInfM = norm(vInfM); % km/s
VInfP = norm(vInfP); % km/s
VInf = (VInfM + VInfP)/2; % km/s

bodies.earth.u = 3.986004415e5; % km^3/s^2

% ------------------------------------
%%% Calculations
% ------------------------------------
k_hat = [0; 0; 1];
S_hat = vInfM ./ VInfM;
T_hat = cross(S_hat,k_hat)/norm(cross(S_hat,k_hat));
R_hat = cross(S_hat,T_hat)/norm(cross(S_hat,T_hat));

h_hat = cross(vInfM,vInfP)/norm(cross(vInfM,vInfP));

B_hat = cross(S_hat,h_hat)/norm(cross(S_hat,h_hat));

turningAngle_rad = acos(dot(vInfM,vInfP)/(VInfM*VInfP));


rp = (bodies.earth.u / (VInf^2))*((1)/(cos((pi-turningAngle_rad)/2)) - 1);

b = (bodies.earth.u / (VInf^2)) * sqrt((1 + VInf*VInf*rp/bodies.earth.u)^2 - 1);

B = B_hat.*b;

B_T = dot(B,T_hat);
B_R = dot(B,R_hat);

theta_rad = acos(dot(T_hat,B_hat));
if (dot(B_hat,R_hat) < 0)
    theta_rad = 2*pi - theta_rad;
end

fprintf('Flyby radius:\t%1.4f km\n\n',rp)
fprintf('Turning angle:\t%1.4f deg\n\n',turningAngle_rad*180/pi)
fprintf('B_T:\t\t%1.4f km\n\n',B_T)
fprintf('B_R:\t\t%1.4f km\n\n',B_R)
fprintf('|B|:\t\t%1.4f km\n\n',norm(B))
fprintf('Theta:\t\t%1.4f deg\n',theta_rad*180/pi)





end % run_p2

% ========================================================================
%%% HW4 - P3
% ========================================================================
if (run_prob == 3)
% ------------------------------------
%%% Bodies
% ------------------------------------
%%% Body shortcuts
Sun   = bodies.sun;
Earth = bodies.earth;
Venus = bodies.venus;

Earth.u = 3.986004415e5;    % km^3/s^2
Sun.u   = 1.32712440018e11; % km^2 / s^2
Venus.u = 3.257e5;          % km^2 / s^2

Venus.R = 6052;                            % km

% ------------------------------------
%%% Julian Dates
% ------------------------------------
%%% Julian Dates
JD_Launch   = 2447807.5; 
JD_VenusFB  = 2447932.5;
JD_EarthFB1 = 2448235.5;
JD_EarthFB2 = 2448966.0;
JD_Jupiter  = 2450164.0;

%%% Relevant TOFs
tof1_sec = (JD_VenusFB - JD_Launch)*86400;
tof2_sec = (JD_EarthFB1 - JD_VenusFB)*86400;

% ------------------------------------
%%% Positions of planets
% ------------------------------------
%%% Find Earth position at departure time
[~, a, e, i_rad, Omega_rad, ~, w_rad, ~, ta_rad] = ...
    getPlanetElements_Meeus(JD_Launch, 'Earth', 'radians');
[rEarth_L, vEarth_L] = OE2ECI(a, e, i_rad, Omega_rad, w_rad, ta_rad, Sun.u);

%%% Find Venus position at flyby time
[~, a, e, i_rad, Omega_rad, ~, w_rad, ~, ta_rad] = ...
    getPlanetElements_Meeus(JD_VenusFB, 'Venus', 'radians');
[rVen_FB, vVen_FB] = OE2ECI(a, e, i_rad, Omega_rad, w_rad, ta_rad, Sun.u);

%%% Find Earth position at FB1 time
[~, a, e, i_rad, Omega_rad, ~, w_rad, ~, ta_rad] = ...
    getPlanetElements_Meeus(JD_EarthFB1, 'Earth', 'radians');
[rEarth_FB1, vEarth_FB1] = OE2ECI(a, e, i_rad, Omega_rad, w_rad, ta_rad, Sun.u);

% ------------------------------------
%%% Lambert solutions
% ------------------------------------
%%% Leg 1
[v_EarthDeparture, v_VenusArrival,~] = lambertSolver(rEarth_L,rVen_FB,tof1_sec,0,0,Sun.u);

%%% Leg 2
[v_VenusDeparture, v_EarthArrival,~] = lambertSolver(rVen_FB,rEarth_FB1,tof2_sec,0,0,Sun.u);


% ------------------------------------
%%% What goes down at Venus
% ------------------------------------
% V_VenusArrival
% V_VenusDeparture

vInfM = v_VenusArrival - vVen_FB;
vInfP = v_VenusDeparture - vVen_FB;

VInfM = norm(vInfM);
VInfP = norm(vInfP);
VInfM - VInfP


dV_Venus = v_VenusDeparture - v_VenusArrival;

turningAngle_rad = acos(dot(vInfM,vInfP)/(norm(vInfM)*norm(vInfP)));
turningAngle_deg = turningAngle_rad*180/pi;


r_flyby = (Venus.u / (VInfM*VInfP)) * ((1/(cos((pi-turningAngle_rad)/2))) - 1)
flybyAltitude = r_flyby - Venus.R

energy_pre  = (norm(v_VenusArrival)^2)/2 - Sun.u/Venus.R;
energy_post = (norm(v_VenusDeparture)^2)/2 - Sun.u/Venus.R;

energy_post - energy_pre % km^2/s^2
end % run_p3

















