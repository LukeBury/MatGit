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


% ========================================================================
%%% Setup
% ========================================================================
% ------------------------------------
%%% Setup - Julian Dates
% ------------------------------------
JD_Launch = 2447814;        % days
JD_VGA1   = 2447932;        % days
JD_EGA1   = 2448235;        % days
JD_EGA2   = 2448965.484378; % days
JD_JOI    = 2450154;        % days

% ------------------------------------
%%% Setup - Bodies
% ------------------------------------
Sun     = bodies.sun;
Venus   = bodies.venus;
Earth   = bodies.earth;
Jupiter = bodies.jupiter; 

% ------------------------------------
%%% Setup - Other
% ------------------------------------
altitude_Earth = 300; % km

% ========================================================================
%%% Preparations (State and Lambert solutions)
% ========================================================================
% ------------------------------------
%%% Transfer TOFs
% ------------------------------------
tof_Launch_2_VGA1 = (JD_VGA1 - JD_Launch)*86400; % sec
tof_VGA1_2_EGA1   = (JD_EGA1 - JD_VGA1)*86400;   % sec
tof_EGA1_2_EGA2   = (JD_EGA2 - JD_EGA1)*86400;   % sec
tof_EGA2_2_JOI    = (JD_JOI - JD_EGA2)*86400;    % sec

% ------------------------------------
%%% Ephemeris positions of planets
% ------------------------------------
%%% Find Earth position at Launch
[~, a, e, i_rad, Omega_rad, ~, w_rad, ~, ta_rad] = ...
    getPlanetElements_Meeus(JD_Launch, 'Earth', 'radians');
[rEarth_Launch, vEarth_Launch] = OE2ECI(a, e, i_rad, Omega_rad, w_rad, ta_rad, Sun.u);

%%% Find Venus position at VGA1 
[~, a, e, i_rad, Omega_rad, ~, w_rad, ~, ta_rad] = ...
    getPlanetElements_Meeus(JD_VGA1, 'Venus', 'radians');
[rVenus_VGA1, vVenus_VGA1] = OE2ECI(a, e, i_rad, Omega_rad, w_rad, ta_rad, Sun.u);

%%% Find Earth position at EGA1 
[~, a, e, i_rad, Omega_rad, ~, w_rad, ~, ta_rad] = ...
    getPlanetElements_Meeus(JD_EGA1, 'Earth', 'radians');
[rEarth_EGA1, vEarth_EGA1] = OE2ECI(a, e, i_rad, Omega_rad, w_rad, ta_rad, Sun.u);

%%% Find Earth position at EGA2
[~, a, e, i_rad, Omega_rad, ~, w_rad, ~, ta_rad] = ...
    getPlanetElements_Meeus(JD_EGA2, 'Earth', 'radians');
[rEarth_EGA2, vEarth_EGA2] = OE2ECI(a, e, i_rad, Omega_rad, w_rad, ta_rad, Sun.u);

%%% Find Jupiter position at JOI
[~, a, e, i_rad, Omega_rad, ~, w_rad, ~, ta_rad] = ...
    getPlanetElements_Meeus(JD_JOI, 'Jupiter', 'radians');
[rJupiter_JOI, vJupiter_JOI] = OE2ECI(a, e, i_rad, Omega_rad, w_rad, ta_rad, Sun.u);

% ------------------------------------
%%% Lambert solutions
% ------------------------------------
%%% Launch to VGA1
[vSc_Launch_p, vSc_VGA1_m, ~] = lambertSolver(rEarth_Launch, rVenus_VGA1,  tof_Launch_2_VGA1, 0, 0, Sun.u);

%%% VGA1 to EGA1
[vSc_VGA1_p, vSc_EGA1_m, ~]   = lambertSolver(rVenus_VGA1,   rEarth_EGA1,  tof_VGA1_2_EGA1,   0, 0, Sun.u);

%%% EGA1 to EGA2
[vSc_EGA1_p, vSc_EGA2_m, ~]   = lambertSolver(rEarth_EGA1,   rEarth_EGA2,  tof_EGA1_2_EGA2,   0, 0, Sun.u);

%%% EGA2 to JOI
[vSc_EGA2_p, vSc_JOI_m, ~]    = lambertSolver(rEarth_EGA2,   rJupiter_JOI, tof_EGA2_2_JOI,    0, 0, Sun.u);

% ========================================================================
%%% Calculations
% ========================================================================
% ------------------------------------
%%% P1 - Hyperbolic excess velocity upon arrival at EGA1
% ------------------------------------
vInf_EGA1_m = vSc_EGA1_m - vEarth_EGA1
VInf_EGA1_m = norm(vInf_EGA1_m);

% ------------------------------------
%%% P2, P3 - Hyperbolic excess velocity upon departure at EGA2
% ------------------------------------
vInf_EGA2_p = vSc_EGA2_p - vEarth_EGA2

norm(vInf_EGA2_p) - norm(vInf_EGA1_m)

% ------------------------------------
%%% P4 - Plotting rp vs phi
% ------------------------------------
%%% phi values
n_phis = 10000;
phis = linspace(0,2*pi,n_phis);

%%% Compute theta
% First Flyby
VSc_EGA1_p = norm(vSc_EGA1_p);
VEarth_EGA1 = norm(vEarth_EGA1);
theta_EGA1 = acos((-VSc_EGA1_p^2 + VInf_EGA1_m^2 + VEarth_EGA1^2)/(2*VInf_EGA1_m*VEarth_EGA1));
% Second Flyby
VSc_EGA2_p = norm(vSc_EGA2_p);
VEarth_EGA2 = norm(vEarth_EGA2);
% vInf_EGA2_m = vSc_EGA2_m - vEarth_EGA2;
% VInf_EGA2_m = norm(vInf_EGA2_m);
% theta_EGA2 = acos((-VSc_EGA2_p^2 + VInf_EGA2_m^2 + VEarth_EGA2^2)/(2*VInf_EGA2_m*VEarth_EGA2));

%%% Loop through phi values
rps_EGA1 = zeros(n_phis,1);
rps_EGA2 = rps_EGA1;
for kk = 1:n_phis
    %%% Set phi
    phi = phis(kk);
    
    % ---------------------
    %%% EGA 1
    % ---------------------
    %%% Compute v-inf plus
    vInf_EGA1_p_VNC = VInf_EGA1_m * [cos(pi-theta_EGA1); sin(pi-theta_EGA1)*cos(phi); -sin(pi-theta_EGA1)*sin(phi)];

    %%% Rotate from VNC to inertial
    vHat_EGA1 = vEarth_EGA1 / VEarth_EGA1;
    hHat_EGA1 = cross(rEarth_EGA1,vEarth_EGA1)/norm(cross(rEarth_EGA1,vEarth_EGA1));
    nHat_EGA1 = hHat_EGA1;
    cHat_EGA1 = cross(vHat_EGA1,nHat_EGA1);
    T_VNC2Ecliptic_EGA1 = [vHat_EGA1, nHat_EGA1, cHat_EGA1];
    vInf_EGA1_p = T_VNC2Ecliptic_EGA1 * vInf_EGA1_p_VNC;
    VInf_EGA1_p = norm(vInf_EGA1_p);
    
    %%% Calculate turning angle
    turningAngle_EGA1 = acos(dot(vInf_EGA1_m,vInf_EGA1_p)/(VInf_EGA1_m*VInf_EGA1_p));
    
    %%% Calculate rp
    rps_EGA1(kk) = (Earth.u/(VInf_EGA1_m*VInf_EGA1_p))*((1/cos((pi-turningAngle_EGA1)/2)) - 1);
    
    % ---------------------
    %%% EGA 2
    % ---------------------
    %%% Theta
    vInf_EGA2_m = vInf_EGA1_p + vEarth_EGA1 - vEarth_EGA2;
%     vInf_EGA2_m = vSc_EGA2_m - vEarth_EGA2;
    VInf_EGA2_m = norm(vInf_EGA2_m);
%     theta_EGA2 = acos((-VSc_EGA2_p^2 + VInf_EGA2_m^2 + VEarth_EGA2^2)/(2*VInf_EGA2_m*VEarth_EGA2));
%     
%     % Compute v-inf plus
%     vInf_EGA2_p_VNC = VInf_EGA2_m * [cos(pi-theta_EGA2); sin(pi-theta_EGA2)*cos(phi); -sin(pi-theta_EGA2)*sin(phi)];
% 
%     %%% Rotate from VNC to inertial
%     vHat_EGA2 = vEarth_EGA2 / VEarth_EGA2;
%     hHat_EGA2 = cross(rEarth_EGA2,vEarth_EGA2)/norm(cross(rEarth_EGA2,vEarth_EGA2));
%     nHat_EGA2 = hHat_EGA2;
%     cHat_EGA2 = cross(vHat_EGA2,nHat_EGA2);
%     T_VNC2Ecliptic_EGA2 = [vHat_EGA2, nHat_EGA2, cHat_EGA2];
%     vInf_EGA2_p = T_VNC2Ecliptic_EGA2 * vInf_EGA2_p_VNC;
    VInf_EGA2_p = norm(vInf_EGA2_p);
    
    %%% Calculate turning angle
    turningAngle_EGA2 = acos(dot(vInf_EGA2_m,vInf_EGA2_p)/(VInf_EGA2_m*VInf_EGA2_p));
    
    %%% Calculate rp
    rps_EGA2(kk) = (Earth.u/(VInf_EGA2_m*VInf_EGA2_p))*((1/cos((pi-turningAngle_EGA2)/2)) - 1);
end

figure; hold all
p4 = fill([73.09, 73.09, 114.5, 114.5], [0, 12000,12000,0], 'g');
p1 = plot(phis.*(180/pi),rps_EGA1,'linewidth',2,'color',colors.std.blue);
p2 = plot(phis.*(180/pi),rps_EGA2,'linewidth',2,'color',colors.std.ltred);
p3 = plot([0,360],[Earth.R+altitude_Earth, Earth.R+altitude_Earth],'--k','linewidth',1);
% p4 = plot([73.09, 73.09],[0 12000],'-g','linewidth',2);
% plot([114.5, 114.5],[0 12000],'-g','linewidth',2)
PlotBoi2('$\Phi$, $^\circ$','Flyby Radius, $km$',18,'LaTex')
xlim([0 360])
legend([p1 p2 p3 p4],'EGA1','EGA2','Minimum r_p','Range of Acceptable \Phi')











































