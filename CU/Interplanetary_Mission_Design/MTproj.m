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
plotNominalTrajectory = 1;

% ========================================================================
%%% Setup
% ========================================================================
% ------------------------------------
%%% Baseline Julian Dates
% ------------------------------------
JD_Launch = 2447807.5; % days
JD_VGA    = 2447932.5;
JD_EGA1   = 2448235.5;
JD_EGA2   = 2448965.984378;
JD_JOI    = 2450164;

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

% ========================================================================
%%% Pre-Resonant trajectories
% ========================================================================
% ------------------------------------
%%% Transfer TOFs
% ------------------------------------
%%% TOFs
tof_Launch2VGA = (JD_VGA - JD_Launch)*86400; % sec
tof_VGA2EGA1 = (JD_EGA1 - JD_VGA)*86400;
tof_EGA12EGA2 = (JD_EGA2 - JD_EGA1)*86400;
tof_EGA22JOI = (JD_JOI - JD_EGA2)*86400;

% ------------------------------------
%%% Ephemeris positions of planets
% ------------------------------------
%%% Earth state at Launch
[rEarth_Launch, vEarth_Launch] = JulianDate_2_HeliocentricState(JD_Launch, 'Earth'); % km, km/s

%%% Venus state at Launch
[rVenus_VGA, vVenus_VGA] = JulianDate_2_HeliocentricState(JD_VGA, 'Venus'); % km, km/s

%%% Earth state at EGA1
[rEarth_EGA1, vEarth_EGA1] = JulianDate_2_HeliocentricState(JD_EGA1, 'Earth'); % km, km/s

%%% Earth state at EGA2
[rEarth_EGA2, vEarth_EGA2] = JulianDate_2_HeliocentricState(JD_EGA2, 'Earth'); % km, km/s

%%% Jupiter state at JOI
[rJupiter_JOI, vJupiter_JOI] = JulianDate_2_HeliocentricState(JD_JOI, 'Jupiter'); % km, km/s

% ------------------------------------
%%% Lambert solutions
% ------------------------------------
%%% Launch to VGA
[vSc_Launch, vSc_m_VGA, ~] = lambertSolver(rEarth_Launch, rVenus_VGA, tof_Launch2VGA, 0, 0, Sun.u); % km/s

%%% VGA to EGA1
[vSc_p_VGA, vSc_m_EGA1, ~] = lambertSolver(rVenus_VGA, rEarth_EGA1, tof_VGA2EGA1, 0, 0, Sun.u); % km/s

%%% EGA1 to EGA2
[vSc_p_EGA1_nom, vSc_m_EGA2_nom, ~] = lambertSolver(rEarth_EGA1, rEarth_EGA2, tof_EGA12EGA2, 0, 0, Sun.u); % km/s

%%% EGA2 to JOI
[vSc_p_EGA2, vSc_m_JOI, ~] = lambertSolver(rEarth_EGA2, rJupiter_JOI, tof_EGA22JOI, 0, 0, Sun.u); % km/s

% ------------------------------------
%%% Launch conditions
% ------------------------------------
%%% vInf and Launch C3
vInf_p_Launch = vSc_Launch - vEarth_Launch; % km/s
VInf_p_Launch = norm(vInf_p_Launch);        % km/s
C3_Launch = VInf_p_Launch^2;                % km^2/s^2

RLA_deg = atan2(vInf_p_Launch(2),vInf_p_Launch(1))*180/pi; % deg
DLA_deg = asin(vInf_p_Launch(3)/VInf_p_Launch)*180/pi;     % deg

% ------------------------------------
%%% VGA conditions
% ------------------------------------
%%% vInf
vInf_m_VGA = vSc_m_VGA - vVenus_VGA; % km/s
vInf_p_VGA = vSc_p_VGA - vVenus_VGA; % km/s

VInf_m_VGA = norm(vInf_m_VGA); % km/s
VInf_p_VGA = norm(vInf_p_VGA); % km/s
dVInf_VGA = VInf_p_VGA - VInf_m_VGA;  %km/s

%%% Turning angle
turningAngle_VGA_rad = acos(dot(vInf_m_VGA,vInf_p_VGA)./(norm(vInf_m_VGA)*norm(vInf_p_VGA)));
turningAngle_VGA_deg = turningAngle_VGA_rad*180/pi;

%%% B-Plane
[BT_VGA, BR_VGA, B_VGA, theta_VGA_rad] = BPlaneTargeter_vInfs(vInf_m_VGA, vInf_p_VGA, Venus.u, 1);

%%% Flyby radius
[ flybyRadius_VGA ] = calcFlybyRadius( turningAngle_VGA_rad, (VInf_m_VGA + VInf_p_VGA)/2, Venus.u ); % km
flybyAltitude_VGA = flybyRadius_VGA - Venus.R; % km

% ========================================================================
%%% Resonant trajectory (EGA1 to EGA2)
% ========================================================================
% ------------------------------------
%%% EGA1 guess conditions
% ------------------------------------
%%% vInf
vInf_m_EGA1 = vSc_m_EGA1 - vEarth_EGA1; % km/s

VInf_m_EGA1 = norm(vInf_m_EGA1); % km/s

% ------------------------------------
%%% EGA2 guess conditions
% ------------------------------------
vInf_p_EGA2 = vSc_p_EGA2 - vEarth_EGA2;
VInf_p_EGA2 = norm(vInf_p_EGA2);

% ------------------------------------
%%% Plotting rp vs phi for resonant earth orbit
% ------------------------------------
%%% phi values
n_phis = 10000;
phis_rad = linspace(0,2*pi,n_phis);

%%% Compute theta
% First Flyby
VSc_p_EGA1 = norm(vSc_p_EGA1_nom);
VEarth_EGA1 = norm(vEarth_EGA1);
theta_EGA1_rad = acos((-VSc_p_EGA1^2 + VInf_m_EGA1^2 + VEarth_EGA1^2)/(2*VInf_m_EGA1*VEarth_EGA1));

% Second Flyby
VSc_EGA2_p = norm(vSc_p_EGA2);
VEarth_EGA2 = norm(vEarth_EGA2);

%%% Loop through phi values
rps_EGA1 = zeros(n_phis,1);
rps_EGA2 = rps_EGA1;
vInfs_p_EGA1 = zeros(n_phis,3);
vInfs_m_EGA2 = vInfs_p_EGA1;

for kk = 1:n_phis
    %%% Set phi
    phi_rad = phis_rad(kk);
    
    % ---------------------
    %%% EGA 1
    % ---------------------
    %%% Compute v-inf plus
    vInf_EGA1_p_VNC = VInf_m_EGA1 * [cos(pi-theta_EGA1_rad); sin(pi-theta_EGA1_rad)*cos(phi_rad); -sin(pi-theta_EGA1_rad)*sin(phi_rad)];

    %%% Rotate from VNC to inertial
    vHat_EGA1 = vEarth_EGA1 / VEarth_EGA1;
    hHat_EGA1 = cross(rEarth_EGA1,vEarth_EGA1)/norm(cross(rEarth_EGA1,vEarth_EGA1));
    nHat_EGA1 = hHat_EGA1;
    cHat_EGA1 = cross(vHat_EGA1,nHat_EGA1);
    T_VNC2Ecliptic_EGA1 = [vHat_EGA1, nHat_EGA1, cHat_EGA1];
    vInf_p_EGA1_kk = T_VNC2Ecliptic_EGA1 * vInf_EGA1_p_VNC;
    VInf_p_EGA1_kk = norm(vInf_p_EGA1_kk);
    vInfs_p_EGA1(kk,:) = vInf_p_EGA1_kk';

    %%% Calculate turning angle
    turningAngle_EGA1_rad = acos(dot(vInf_m_EGA1,vInf_p_EGA1_kk)/(VInf_m_EGA1*VInf_p_EGA1_kk));
    
    %%% Calculate rp
    rps_EGA1(kk) = (Earth.u/(VInf_m_EGA1*VInf_p_EGA1_kk))*((1/cos((pi-turningAngle_EGA1_rad)/2)) - 1);
    
    % ---------------------
    %%% EGA 2
    % ---------------------
    %%% Theta
    vInf_m_EGA2_kk = vInf_p_EGA1_kk + vEarth_EGA1 - vEarth_EGA2;
    VInf_m_EGA2_kk = norm(vInf_m_EGA2_kk);
    vInfs_m_EGA2(kk,:) = vInf_m_EGA2_kk';
    
    %%% Calculate turning angle
    turningAngle_EGA2 = acos(dot(vInf_m_EGA2_kk,vInf_p_EGA2)/(VInf_m_EGA2_kk*VInf_p_EGA2));
    
    %%% Calculate rp
    rps_EGA2(kk) = (Earth.u/(VInf_m_EGA2_kk*VInf_p_EGA2))*((1/cos((pi-turningAngle_EGA2)/2)) - 1);
end

%%% Grabbing desired results
myPhi_index = find(rps_EGA2 == max(rps_EGA2));
myPhi_rad = phis_rad(myPhi_index);
myPhi_deg = myPhi_rad*180/pi;
vInf_p_EGA1 = vInfs_p_EGA1(myPhi_index,:)';
vInf_m_EGA2 = vInfs_m_EGA2(myPhi_index,:)';

VInf_p_EGA1 = norm(vInf_p_EGA1);
VInf_m_EGA2 = norm(vInf_m_EGA2);

vSc_p_EGA1 = vInf_p_EGA1 + vEarth_EGA1;
vSc_m_EGA2 = vInf_m_EGA2 + vEarth_EGA2;

%%% Plotting
figure; hold all
p4 = fill([72.08, 72.08, 114.52, 114.52], [0, 12000,12000,0], 'g');
p1 = plot(phis_rad.*(180/pi),rps_EGA1,'linewidth',2,'color',colors.std.blue);
p2 = plot(phis_rad.*(180/pi),rps_EGA2,'linewidth',2,'color',colors.std.ltred);
% p3 = plot([0,360],[Earth.R, Earth.R],'--k','linewidth',1);
p3 = plot([0,360],[Earth.R+300, Earth.R+300],'--k','linewidth',1);
p5 = plot([myPhi_deg,myPhi_deg],[0, 12000],'m','linewidth',1.5);
PlotBoi2('$\Phi$, $^\circ$','Flyby Radius, $km$',18,'LaTex')
xlim([0 360])
% legend([p1 p2 p3],'EGA1','EGA2','300 km Altitude')
legend([p1 p2 p3 p4 p5],'EGA1','EGA2','300 km Altitude','Range of Acceptable \Phi', 'Selected \Phi')

% ------------------------------------
%%% Wrapping up EGA1 and EGA2
% ------------------------------------
%%% Ensuring that VInfs are close
DVInf_EGA1 = VInf_m_EGA1 - VInf_p_EGA1;
DVInf_EGA2 = VInf_m_EGA2 - VInf_p_EGA2;

%%% Flyby radius
flybyRadius_EGA1 = rps_EGA1(myPhi_index); % km
flybyRadius_EGA2 = rps_EGA2(myPhi_index); % km

%%% B-Plane
[BT_EGA1, BR_EGA1, B_EGA1, theta_EGA1_rad] = BPlaneTargeter_vInfs(vInf_m_EGA1, vInf_p_EGA1, Earth.u, 1);
[BT_EGA2, BR_EGA2, B_EGA2, theta_EGA2_rad] = BPlaneTargeter_vInfs(vInf_m_EGA2, vInf_p_EGA2, Earth.u, 1);

% ========================================================================
%%% Post resonant (headed to Jupiter)
% ========================================================================
% ------------------------------------
%%% Galileo info
% ------------------------------------
%http://planet4589.org/space/jsr/notes/orbitaljup.html
i_JOI_deg = 5.300; % deg
a_JOI = 9814688.8; % km
e_JOI = 0.971233;  % km
raan_JOI_deg = 60.964; % deg
w_JOI_deg = 87.684; % deg
TP_JOI_sec = 2*pi*sqrt((a_JOI^3)/Jupiter.u); % sec
TP_JOI_day = TP_JOI_sec/86400; % days
rp_JOI = a_JOI*(1-e_JOI);

% https://pds.nasa.gov/ds-view/pds/viewMissionProfile.jsp?MISSION_NAME=GALILEO
dV_JOI_galileo = 0.630; % km/s

% ------------------------------------
%%% JOI Conditions
% ------------------------------------
%%% vInf
vInf_m_JOI = vSc_m_JOI - vJupiter_JOI; % km/s
VInf_m_JOI = norm(vInf_m_JOI); % km/s
energy_JOI = VInf_m_JOI*VInf_m_JOI/2; % km^2/s^2
Vp_m_JOI = sqrt(2*(energy_JOI + Jupiter.u/rp_JOI)); % km/s

[ Vp_p_JOI ] = visviva_v( rp_JOI, a_JOI, Jupiter.u);  %km/s

dV_JOI = Vp_m_JOI - Vp_p_JOI;

% ========================================================================
%%% Plot Nominal Trajectory
% ========================================================================
if plotNominalTrajectory == 1
    % ------------------------------------
    %%% ICs
    % ------------------------------------
    %%% Transfer legs
    X0_leg1   = [rEarth_Launch; vSc_Launch]; % Launch to VGA
    X0_leg2   = [rVenus_VGA; vSc_p_VGA]; % VGA to EGA1
    X0_leg3   = [rEarth_EGA1; vSc_p_EGA1]; % EGA1 to EGA2
    X0_leg4   = [rEarth_EGA2; vSc_p_EGA2]; % EGA2 to JOI

    %%% Post JOI
    [r0_JupOrb, v0_JupOrb] = COE2RV(a_JOI, e_JOI, i_JOI_deg, raan_JOI_deg, w_JOI_deg, 0, Jupiter.u);
    X0_JupOrb = [r0_JupOrb; v0_JupOrb]; % Post JOI-burn
    
    
    % ------------------------------------
    %%% Building time vectors
    % ------------------------------------
    tf_leg1 = tof_Launch2VGA;
    tf_leg2 = tof_VGA2EGA1;
    tf_leg3 = tof_EGA12EGA2;
    tf_leg4 = tof_EGA22JOI;
    

    t0 = 0; % sec
    n_dt = 1000; 
    time0_leg1 = linspace(t0,tf_leg1,n_dt);
    time0_leg2 = linspace(t0,tf_leg2,n_dt);
    time0_leg3 = linspace(t0,tf_leg3,n_dt);
    time0_leg4 = linspace(t0,tf_leg4,n_dt);

    % ------------------------------------
    %%% Integration
    % ------------------------------------
    %%% Choosing ode tolerance
    tol = 1e-13;

    %%% Setting integrator options
    options = odeset('RelTol',tol,'AbsTol',tol);

    %%% Integrating
    [time_leg1, X_2BI_leg1] = ode113(@Int_2BI, time0_leg1, X0_leg1, options, Sun.u);
    [time_leg2, X_2BI_leg2] = ode113(@Int_2BI, time0_leg2, X0_leg2, options, Sun.u);
    [time_leg3, X_2BI_leg3] = ode113(@Int_2BI, time0_leg3, X0_leg3, options, Sun.u);
    [time_leg4, X_2BI_leg4] = ode113(@Int_2BI, time0_leg4, X0_leg4, options, Sun.u);

    % ------------------------------------
    %%% Plotting
    % ------------------------------------
    nominalTrajColors = colorScale([colors.std.ltred; colors.std.ltblue],4);

    lw = 5;

    figure; hold all
    p1 = plot3(X_2BI_leg1(:,1),X_2BI_leg1(:,2),X_2BI_leg1(:,3),'color',nominalTrajColors(1,:),'linewidth',lw);
    p2 = plot3(X_2BI_leg2(:,1),X_2BI_leg2(:,2),X_2BI_leg2(:,3),'color',nominalTrajColors(2,:),'linewidth',lw);
    p3 = plot3(X_2BI_leg3(:,1),X_2BI_leg3(:,2),X_2BI_leg3(:,3),'color',nominalTrajColors(3,:),'linewidth',lw);
    p4 = plot3(X_2BI_leg4(:,1),X_2BI_leg4(:,2),X_2BI_leg4(:,3),'color',nominalTrajColors(4,:),'linewidth',lw);

    PlotBoi3('$X$, km','$Y$, km','$Z$, km',18,'LaTex')
    axis equal
    view(0,90)

    legend([p1 p2 p3 p4],'Launch to VGA','VGA to EGA1','EGA1 to EGA2','EGA2 to JOI')
end % plotNominalTrajectory == 1

% ========================================================================
%%% Pork Chop Plot - Launch to VGA (leg1)
% ========================================================================
% ------------------------------------
%%% Setting up julian dates
% ------------------------------------
nPoints = 150;
JD_Departures_Launch_vec = linspace(JD_Launch-30, JD_Launch+30, nPoints);
JD_Arrivals_VGA_vec     = linspace(JD_VGA-30, JD_VGA+30, nPoints);

[JD_Departures_Launch, JD_Arrivals_VGA] = meshgrid(JD_Departures_Launch_vec, JD_Arrivals_VGA_vec);

% ------------------------------------
%%% Looping through options and computing vInf and c3
% ------------------------------------
X0s_Sc_leg1{size(JD_Departures_Launch,1),size(JD_Departures_Launch,2)} = [];
vInfs_VGA_m{size(JD_Departures_Launch,1),size(JD_Departures_Launch,2)} = [];
VInfs_VGA_m = zeros(size(JD_Departures_Launch));
c3s_Launch_p = zeros(size(JD_Departures_Launch));
TOFs_days = zeros(size(JD_Departures_Launch));

for xIndex = 1:nPoints
    for yIndex = 1:nPoints
        %%% Grab current launch window
        JD_D = JD_Departures_Launch(xIndex, yIndex);
        JD_A = JD_Arrivals_VGA(xIndex, yIndex);
        
        %%% Find Earth position at departure time
        [rEarth_D, vEarth_D] = JulianDate_2_HeliocentricState(JD_D, 'Earth');
        
        %%% Find Venus position at arrival time
        [rVenus_A, vVenus_A] = JulianDate_2_HeliocentricState(JD_A, 'Venus');
        
        %%% Use Lambert solver to find solution
        tof_sec = (JD_A - JD_D)*86400;
        [vSc_EarthDeparture, vSc_VenusArrival,~] = lambertSolver(rEarth_D,rVenus_A,tof_sec,0,0,Sun.u);
        
        %%% C3 at Earth departure
        vInf_Earth = vSc_EarthDeparture - vEarth_D;
        c3s_Launch_p(xIndex, yIndex) = norm(vInf_Earth)^2;
        
        %%% vInf at Venus arrival
        vInf_Venus_m = vSc_VenusArrival - vVenus_A;
        vInfs_VGA_m{xIndex, yIndex} = vInf_Venus_m;
        VInfs_VGA_m(xIndex, yIndex) = norm(vInf_Venus_m);
        
        %%% TOFs
        TOFs_days(xIndex, yIndex) = JD_A - JD_D;
        
        %%% Starting state for this leg of the trajectory [r0; v0]
        X0s_Sc_leg1{xIndex, yIndex} = [rEarth_D; vSc_EarthDeparture];

    end
end



% ------------------------------------
%%% Pork Chop Plot
% ------------------------------------
figure('position',[74 233 762 564]); hold all
vInf_contours = [4.5, 5, 5.5, 6, 6.5, 7, 7.5, 10, 20];
[cs1, h1] = contour(JD_Departures_Launch-JD_Departures_Launch(1),JD_Arrivals_VGA-JD_Arrivals_VGA(1),VInfs_VGA_m,vInf_contours,'ShowText','on','color',colors.std.blue,'linewidth',1.2);

c3_contours = [5, 10, 15, 17.5, 20, 25, 30, 50, 200];
[cs2, h2] = contour(JD_Departures_Launch-JD_Departures_Launch(1),JD_Arrivals_VGA-JD_Arrivals_VGA(1),c3s_Launch_p,c3_contours,'ShowText','on','color',colors.std.ltred,'linewidth',1.2);

TOF_contours = [80, 100, 120, 140, 160, 200];
[cs3,h3] = contour(JD_Departures_Launch-JD_Departures_Launch(1),JD_Arrivals_VGA-JD_Arrivals_VGA(1),TOFs_days,TOF_contours,'ShowText','on','color',colors.std.black,'linewidth',1.2);

p1 = plot(30, 30,'o','markeredgecolor',colors.std.black,'markerfacecolor',colors.std.cyan,'markersize',10);

[yr, mo, d, hr, min, sec] = julian2calendar(JD_Departures_Launch(1));
calendar0_launch = sprintf('Days Past %d/%d/%d',mo,d,yr);

[yr, mo, d, hr, min, sec] = julian2calendar(JD_Arrivals_VGA(1));
calendar0_VGA = sprintf('Days Past %d/%d/%d',mo,d,yr);

PlotBoi2(calendar0_launch,calendar0_VGA,16)
grid off
clabel(cs1 ,h1,'color',colors.std.blue)
clabel(cs2 ,h2,'color',colors.std.ltred)
clabel(cs3 ,h3,'color',colors.std.black)

legend([h1, h2, h3, p1],'V_\infty @ VGA, km/s','C_3 @ Launch, km^2/s^2','TOF, days','Galileo')




% ========================================================================
%%% Pork Chop Plot - VGA to EGA1 (leg2)
% ========================================================================

% ------------------------------------
%%% Setting up julian dates
% ------------------------------------
JD_Departures_VGA_vec = linspace(JD_VGA-30, JD_VGA+30, nPoints);
JD_Arrivals_EGA1_vec  = linspace(JD_EGA1-30, JD_EGA1+30, nPoints);

[JD_Departures_VGA, JD_Arrivals_EGA1] = meshgrid(JD_Departures_VGA_vec, JD_Arrivals_EGA1_vec);

% ------------------------------------
%%% Looping through options and computing vInf and c3
% ------------------------------------
X0_sc_leg2{size(JD_Departures_VGA,1),size(JD_Departures_VGA,2)} = [];
VInfs_EGA1_m = zeros(size(JD_Departures_VGA));
VInfs_VGA_p = zeros(size(JD_Departures_VGA));
% flybyRadii_VGA = zeros(size(JD_Departures_VGA));
TOFs_days = zeros(size(JD_Departures_VGA));

for xIndex = 1:nPoints
    for yIndex = 1:nPoints
        %%% Grab current launch window
        JD_D = JD_Departures_VGA(xIndex, yIndex);
        JD_A = JD_Arrivals_EGA1(xIndex, yIndex);
        
        %%% Find Venus position at departure time
        [rVenus_D, vVenus_D] = JulianDate_2_HeliocentricState(JD_D, 'Venus');
        
        %%% Find Earth position at arrival time
        [rEarth_A, vEarth_A] = JulianDate_2_HeliocentricState(JD_A, 'Earth');
        
        %%% Use Lambert solver to find solution
        tof_sec = (JD_A - JD_D)*86400;
        [vSc_VenusDeparture, vSc_EarthArrival,~] = lambertSolver(rVenus_D,rEarth_A,tof_sec,0,0,Sun.u);
        
        %%% vInf at Venus departure
        vInf_Venus_p = vSc_VenusDeparture - vVenus_D;
        VInfs_VGA_p(xIndex, yIndex) = norm(vInf_Venus_p);
        
%         %%% flyby radius at Venus
%         vInf_Venus_m = vInfs_VGA_m{xIndex, yIndex};
%         VInf = (norm(vInf_Venus_m) + norm(vInf_Venus_p))/2;
%         turningAngle_rad = acos(dot(vInf_Venus_m,vInf_Venus_p)./(norm(vInf_Venus_m)*norm(vInf_Venus_p)));
%         [ rp ] = calcFlybyRadius( turningAngle_rad, VInf, Venus.u );
%         flybyRadii_VGA(xIndex, yIndex) = rp;
        
        %%% vInf at Earth arrival
        vInf_Earth = vSc_EarthArrival - vEarth_A;
        VInfs_EGA1_m(xIndex, yIndex) = norm(vInf_Earth);
        
        %%% TOFs
        TOFs_days(xIndex, yIndex) = JD_A - JD_D;
        
        %%% Starting state for this leg of the trajectory [r0; v0]
        X0_sc_leg2{xIndex, yIndex} = [rVenus_D; vSc_VenusDeparture];

    end
end


% ------------------------------------
%%% Pork Chop Plot
% ------------------------------------
figure('position',[440 234 762 564]); hold all
vInf_VGA_contours = [4.4,4.6, 5, 6, 8, 10, 15, 20, 30];
[cs1, h1] = contour(JD_Departures_VGA-JD_Departures_VGA(1),JD_Arrivals_EGA1-JD_Arrivals_EGA1(1),VInfs_VGA_p,vInf_VGA_contours,'ShowText','on','color',colors.std.ltred,'linewidth',1.2);

vInf_EGA1_contours = [7.7, 8, 8.5, 9, 10, 12, 15, 20, 30];
[cs2, h2] = contour(JD_Departures_VGA-JD_Departures_VGA(1),JD_Arrivals_EGA1-JD_Arrivals_EGA1(1),VInfs_EGA1_m,vInf_EGA1_contours,'ShowText','on','color',colors.std.blue,'linewidth',1.2);

TOF_contours = [260, 280, 300, 320, 340];
[cs3,h3] = contour(JD_Departures_VGA-JD_Departures_VGA(1),JD_Arrivals_EGA1-JD_Arrivals_EGA1(1),TOFs_days,TOF_contours,'ShowText','on','color',colors.std.black,'linewidth',1.2);

p1 = plot(30, 30,'o','markeredgecolor',colors.std.black,'markerfacecolor',colors.std.cyan,'markersize',10);

[yr, mo, d, hr, min, sec] = julian2calendar(JD_Departures_VGA(1));
calendar0_VGA = sprintf('Days Past %d/%d/%d',mo,d,yr);

[yr, mo, d, hr, min, sec] = julian2calendar(JD_Arrivals_EGA1(1));
calendar0_VGA = sprintf('Days Past %d/%d/%d',mo,d,yr);

PlotBoi2(calendar0_VGA,calendar0_VGA,16)
grid off
clabel(cs1 ,h1,'color',colors.std.blue)
clabel(cs2 ,h2,'color',colors.std.ltred)
clabel(cs3 ,h3,'color',colors.std.black)

legend([h1, h2, h3, p1],'V_\infty @ VGA, km/s','V_\infty @ EGA1, km/s','TOF, days','Galileo')






% ========================================================================
%%% Pork Chop Plot - EGA2 to JOI (leg4)
% ========================================================================

% ------------------------------------
%%% Setting up julian dates
% ------------------------------------
JD_Departures_EGA2_vec = linspace(JD_EGA2-30, JD_EGA2+30, nPoints);
JD_Arrivals_JOI_vec  = linspace(JD_JOI-30, JD_JOI+30, nPoints);

[JD_Departures_EGA2, JD_Arrivals_JOI] = meshgrid(JD_Departures_EGA2_vec, JD_Arrivals_JOI_vec);

% ------------------------------------
%%% Looping through options and computing vInf and c3
% ------------------------------------
X0_sc_leg3{size(JD_Departures_EGA2,1),size(JD_Departures_EGA2,2)} = [];
VInfs_JOI_m = zeros(size(JD_Departures_EGA2));
VInfs_EGA2_p = zeros(size(JD_Departures_EGA2));
% flybyRadii_EGA2 = zeros(size(JD_Departures_EGA2));
TOFs_days = zeros(size(JD_Departures_EGA2));

for xIndex = 1:nPoints
    for yIndex = 1:nPoints
        %%% Grab current launch window
        JD_D = JD_Departures_EGA2(xIndex, yIndex);
        JD_A = JD_Arrivals_JOI(xIndex, yIndex);
        
        %%% Find Earth position at departure time
        [rEarth_D, vEarth_D] = JulianDate_2_HeliocentricState(JD_D, 'Earth');
        
        %%% Find Earth position at arrival time
        [rJupiter_A, vJupiter_A] = JulianDate_2_HeliocentricState(JD_A, 'Jupiter');
        
        %%% Use Lambert solver to find solution
        tof_sec = (JD_A - JD_D)*86400;
        [vSc_EarthDeparture, vSc_JupiterArrival,~] = lambertSolver(rEarth_D,rJupiter_A,tof_sec,0,0,Sun.u);
        
        %%% vInf at Earth departure
        vInf_Earth_p = vSc_EarthDeparture - vEarth_D;
        VInfs_EGA2_p(xIndex, yIndex) = norm(vInf_Earth_p);
        
%         %%% flyby radius at Earth
%         vInf_Earth_m = vInfs_EGA2_m{xIndex, yIndex};
%         VInf = (norm(vInf_Earth_m) + norm(vInf_Earth_p))/2;
%         turningAngle_rad = acos(dot(vInf_Earth_m,vInf_Earth_p)./(norm(vInf_Earth_m)*norm(vInf_Earth_p)));
%         [ rp ] = calcFlybyRadius( turningAngle_rad, VInf, Venus.u );
%         flybyRadii_EGA2(xIndex, yIndex) = rp;
%         
        %%% vInf at Jupiter arrival
        vInf_Jupiter = vSc_JupiterArrival - vJupiter_A;
        VInfs_JOI_m(xIndex, yIndex) = norm(vInf_Jupiter);
        
        %%% TOFs
        TOFs_days(xIndex, yIndex) = JD_A - JD_D;
        
        %%% Starting state for this leg of the trajectory [r0; v0]
        X0_sc_leg3{xIndex, yIndex} = [rEarth_D; vSc_EarthDeparture];

    end
end


% ------------------------------------
%%% Pork Chop Plot
% ------------------------------------
figure('position',[440 234 762 564]); hold all
vInf_EGA2_contours = [8.8, 9, 9.5, 10, 11, 12, 20];
[cs1, h1] = contour(JD_Departures_EGA2-JD_Departures_EGA2(1),JD_Arrivals_JOI-JD_Arrivals_JOI(1),VInfs_EGA2_p,vInf_EGA2_contours,'ShowText','on','color',colors.std.ltred,'linewidth',1.2);

vInf_JOI_contours = [5.8, 5.85, 5.9, 6, 7];
[cs2, h2] = contour(JD_Departures_EGA2-JD_Departures_EGA2(1),JD_Arrivals_JOI-JD_Arrivals_JOI(1),VInfs_JOI_m,vInf_JOI_contours,'ShowText','on','color',colors.std.blue,'linewidth',1.2);
% 
TOF_contours = [1160, 1180, 1200, 1220, 1240];
[cs3,h3] = contour(JD_Departures_EGA2-JD_Departures_EGA2(1),JD_Arrivals_JOI-JD_Arrivals_JOI(1),TOFs_days,TOF_contours,'ShowText','on','color',colors.std.black,'linewidth',1.2);

p1 = plot(30, 30,'o','markeredgecolor',colors.std.black,'markerfacecolor',colors.std.cyan,'markersize',10);

[yr, mo, d, hr, min, sec] = julian2calendar(JD_Departures_EGA2(1));
calendar0_EGA2 = sprintf('Days Past %d/%d/%d',mo,d,yr);

[yr, mo, d, hr, min, sec] = julian2calendar(JD_Arrivals_JOI(1));
calendar0_JOI = sprintf('Days Past %d/%d/%d',mo,d,yr);

PlotBoi2(calendar0_EGA2,calendar0_JOI,16)
grid off
clabel(cs1 ,h1,'color',colors.std.blue)
clabel(cs2 ,h2,'color',colors.std.ltred)
clabel(cs3 ,h3,'color',colors.std.black)

legend([h1, h2, h3, p1],'V_\infty @ EGA2, km/s','V_\infty @ JOI, km/s','TOF, days','Galileo')









% [ JD ] = calendar2julian(yr, mo, d, h, min, s)

[ JD_new ] = calendar2julian(1996, 4, 4, 14, 53, 21.339)
















