clear
clc
close all
mbinPath = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/mbin';
addpath(genpath(mbinPath))
tic
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
JD_Launch = 2453755.29167; % days
JD_JGA    = 2454159.73681; % days
JD_Pluto  = 2457217.99931; % days

[ JD_Launch_NH ] = calendar2julian(2006, 1, 19, 19, 0, 0);
[ JD_JGA_NH ] = calendar2julian(2007, 2, 28, 5, 41, 0);
[ JD_Pluto_NH ] = calendar2julian(2015, 07, 14, 11, 59, 0);

% ------------------------------------
%%% Setup - Bodies
% ------------------------------------
Sun     = bodies.sun;
Earth   = bodies.earth;
Jupiter = bodies.jupiter; 

% ------------------------------------
%%% Setup - Other
% ------------------------------------
Jupiter.u = 1.266865361e8; % km^3/s^2
Jupiter.R = 71492; % km

%% =======================================================================
%%% Part 1 - Porkchop Plot Creation
% ========================================================================
% ------------------------------------
%%% Transfer TOFs
% ------------------------------------
%%% TOFs
tof_Launch2JGA = (JD_JGA - JD_Launch)*86400; % sec
tof_JGA2Pluto  = (JD_Pluto - JD_JGA)*86400;  % sec

% ------------------------------------
%%% Ephemeris positions of planets
% ------------------------------------
%%% Earth state at Launch
[rEarth_Launch, vEarth_Launch] = JulianDate_2_HeliocentricState(JD_Launch, 'Earth');

%%% Jupiter state at JGA
[rJupiter_JGA, vJupiter_JGA] = JulianDate_2_HeliocentricState(JD_JGA, 'Jupiter');

%%% Pluto state at flyby
[rPluto_flyby, vPluto_flyby] = JulianDate_2_HeliocentricState(JD_Pluto, 'Pluto');

% ------------------------------------
%%% Lambert solutions
% ------------------------------------
%%% Launch-to-JGA
[vSc_Launch, vSc_JGA_m, ~] = lambertSolver(rEarth_Launch, rJupiter_JGA, tof_Launch2JGA, 0, 0, Sun.u);
% [vSc_Launch,vSc_JGA_m] = lambertTargeting(rEarth_Launch, rJupiter_JGA, tof_Launch2JGA,1,1,Sun.u,0);

%%% Launch-to-JGA
[vSc_JGA_p, vSc_Pluto_m, ~] = lambertSolver(rJupiter_JGA, rPluto_flyby, tof_JGA2Pluto, 0, 0, Sun.u);
% [vSc_JGA_p,vSc_Pluto_m] = lambertTargeting(rJupiter_JGA, rPluto_flyby, tof_JGA2Pluto,1,1,Sun.u,0);

% ------------------------------------
%%% Launch conditions
% ------------------------------------
%%% Launch C3
vInf_Launch = vSc_Launch - vEarth_Launch; % km/s
VInf_Launch = norm(vInf_Launch);          % km/s
C3_Launch = VInf_Launch^2; % km^2/s^2

RLA_deg = atan2(vInf_Launch(2),vInf_Launch(1))*180/pi; % deg
DLA_deg = asin(vInf_Launch(3)/VInf_Launch)*180/pi;     % deg

% ------------------------------------
%%% Jupiter arrival conditions
% ------------------------------------
vInf_JGA_m = vSc_JGA_m - vJupiter_JGA; % km/s
VInf_JGA_m = norm(vInf_JGA_m);         % km/s

% ------------------------------------
%%% Jupiter departure conditions
% ------------------------------------
vInf_JGA_p = vSc_JGA_p - vJupiter_JGA; % km/s
VInf_JGA_p = norm(vInf_JGA_p);         % km/s

diff_vInf = VInf_JGA_p - VInf_JGA_m;

dV_JupiterFlyby = norm(vSc_JGA_p - vSc_JGA_m);

% ------------------------------------
%%% Pluto arrival conditions
% ------------------------------------
vInf_Pluto_m = vSc_Pluto_m - vPluto_flyby; % km/s
VInf_Pluto_m = norm(vInf_Pluto_m);         % km/s

% ------------------------------------
%%% Jupiter B-Plane Parameters
% ------------------------------------
%%% Computing B-Plane Parameters
[BT, BR, B, theta_rad] = BPlaneTargeter_vInfs(vInf_JGA_m, vInf_JGA_p, Jupiter.u, 2);

%%% Turning angle
turningAngle = acos(dot(vInf_JGA_m,vInf_JGA_p)/(VInf_JGA_m*VInf_JGA_p));
turningAngle_deg = turningAngle*180/pi;

%%% Jupiter flyby altitude
[ jupiterFlybyRadius ] = calcFlybyRadius( turningAngle, (VInf_JGA_m + VInf_JGA_p)/2, Jupiter.u );
jupiterFlybyAltitude = jupiterFlybyRadius - Jupiter.R;

% ========================================================================
%%% Porkchop Plot - Earth to Jupiter
% ========================================================================
% ------------------------------------
%%% Setting up julian dates
% ------------------------------------
nPoints = 300;
JD_Dep = linspace(2453714.5, 2453794.5, nPoints);
JD_Arr   = linspace(2454129.5, 2454239.5, nPoints);

[JD_Departures_1, JD_Arrivals_1] = meshgrid(JD_Dep, JD_Arr);

% ------------------------------------
%%% Looping through options and computing vInf and c3
% ------------------------------------
X0_sc_leg1{size(JD_Departures_1,1),size(JD_Departures_1,2)} = [];
vInfs_jup_m{size(JD_Departures_1,1),size(JD_Departures_1,2)} = [];
VInfs_jup_m = zeros(size(JD_Departures_1));
c3s_earth_p = zeros(size(JD_Departures_1));
TOFs_days = zeros(size(JD_Departures_1));

for xIndex = 1:nPoints
    for yIndex = 1:nPoints
        %%% Grab current launch window
        JD_D = JD_Departures_1(xIndex, yIndex);
        JD_A = JD_Arrivals_1(xIndex, yIndex);
        
        %%% Find Earth position at departure time
        [rEarth_D, vEarth_D] = JulianDate_2_HeliocentricState(JD_D, 'Earth');
        
        %%% Find Jupiter position at arrival time
        [rJupiter_A, vJupiter_A] = JulianDate_2_HeliocentricState(JD_A, 'Jupiter');
        
        %%% Use Lambert solver to find solution
        tof_sec = (JD_A - JD_D)*86400;
        [vSc_EarthDeparture, vSc_JupiterArrival,~] = lambertSolver(rEarth_D,rJupiter_A,tof_sec,0,0,Sun.u);
        
        %%% C3 at Earth departure
        vInf_Earth = vSc_EarthDeparture - vEarth_D;
        c3s_earth_p(xIndex, yIndex) = norm(vInf_Earth)^2;
        
        %%% vInf at Jupiter arrival
        vInf_Jupiter_m = vSc_JupiterArrival - vJupiter_A;
        vInfs_jup_m{xIndex, yIndex} = vInf_Jupiter_m;
        VInfs_jup_m(xIndex, yIndex) = norm(vInf_Jupiter_m);
        
        %%% TOFs
        TOFs_days(xIndex, yIndex) = JD_A - JD_D;
        
        %%% Starting state for this leg of the trajectory [r0; v0]
        X0_sc_leg1{xIndex, yIndex} = [rEarth_D; vSc_EarthDeparture];

    end
end

% ------------------------------------
%%% Pork Chop Plot
% ------------------------------------
figure('position',[440 234 762 564]); hold all
vInf_contours = [14, 15, 16, 17, 18, 18.51, 19, 20, 21];
[cs1, h1] = contour(JD_Departures_1-JD_Departures_1(1),JD_Arrivals_1-JD_Arrivals_1(1),VInfs_jup_m,vInf_contours,'ShowText','on','color',colors.std.blue,'linewidth',1.2);

c3_contours = [120, 130, 140, 150, 158.4, 180, 210, 260];
[cs2, h2] = contour(JD_Departures_1-JD_Departures_1(1),JD_Arrivals_1-JD_Arrivals_1(1),c3s_earth_p,c3_contours,'ShowText','on','color',colors.std.ltred,'linewidth',1.2);

TOF_contours = [340, 360, 380, 400, 404.5, 420, 440, 460, 480, 500, 520];
[cs3,h3] = contour(JD_Departures_1-JD_Departures_1(1),JD_Arrivals_1-JD_Arrivals_1(1),TOFs_days,TOF_contours,'ShowText','on','color',colors.std.black,'linewidth',1.2);

p1 = plot(JD_Launch_NH - JD_Departures_1(1), JD_JGA_NH - JD_Arrivals_1(1),'o','markeredgecolor',colors.std.black,'markerfacecolor',colors.std.cyan,'markersize',10);

PlotBoi2('Days Past 10 December 2005','Days Past 29 January 2007',16)
grid off
clabel(cs1 ,h1,'color',colors.std.blue)
clabel(cs2 ,h2,'color',colors.std.ltred)
clabel(cs3 ,h3,'color',colors.std.black)

legend([h1, h2, h3, p1],'V_\infty @ Jupiter, km/s','C_3, km^2/s^2','TOF, days','New Horizons')

% ========================================================================
%%% Porkchop Plot - Jupiter to Pluto
% ========================================================================
% ------------------------------------
%%% Setting up julian dates
% ------------------------------------
% nPoints_2 = 150;
JD_Dep  = linspace(2454129.5, 2454239.5, nPoints);
JD_Arr  = linspace(2456917.5, 2457517.5, nPoints);

[JD_Departures_2, JD_Arrivals_2] = meshgrid(JD_Dep, JD_Arr);

% ------------------------------------
%%% Looping through options and computing vInf and c3
% ------------------------------------
X0_sc_leg2{size(JD_Departures_1,1),size(JD_Departures_1,2)} = [];
VInfs_plu_m = zeros(size(JD_Departures_2));
VInfs_jup_p = zeros(size(JD_Departures_2));
flybyRadii = zeros(size(JD_Departures_2));
TOFs_days = zeros(size(JD_Departures_2));

for xIndex = 1:nPoints
    for yIndex = 1:nPoints
        %%% Grab current launch window
        JD_D = JD_Departures_2(xIndex, yIndex);
        JD_A = JD_Arrivals_2(xIndex, yIndex);
        
        %%% Find Jupiter position at departure time
        [rJupiter_D, vJupiter_D] = JulianDate_2_HeliocentricState(JD_D, 'Jupiter');
        
        %%% Find Jupiter position at arrival time
        [rPluto_A, vPluto_A] = JulianDate_2_HeliocentricState(JD_A, 'Pluto');
        
        %%% Use Lambert solver to find solution
        tof_sec = (JD_A - JD_D)*86400;
        [vSc_JupiterDeparture, vSc_PlutoArrival,~] = lambertSolver(rJupiter_D,rPluto_A,tof_sec,0,0,Sun.u);
        
        %%% vInf at Jupiter departure
        vInf_Jupiter_p = vSc_JupiterDeparture - vJupiter_D;
        VInfs_jup_p(xIndex, yIndex) = norm(vInf_Jupiter_p);
        
        %%% flyby radius at jupiter
        vInf_Jupiter_m = vInfs_jup_m{xIndex, yIndex};
        VInf = (norm(vInf_Jupiter_m) + norm(vInf_Jupiter_p))/2;
        turningAngle_rad = acos(dot(vInf_Jupiter_m,vInf_Jupiter_p)./(norm(vInf_Jupiter_m)*norm(vInf_Jupiter_p)));
        [ rp ] = calcFlybyRadius( turningAngle_rad, VInf, Jupiter.u );
        flybyRadii(xIndex, yIndex) = rp;
        
        %%% vInf at Pluto arrival
        vInf_Pluto = vSc_PlutoArrival - vPluto_A;
        VInfs_plu_m(xIndex, yIndex) = norm(vInf_Pluto);
        
        %%% TOFs
        TOFs_days(xIndex, yIndex) = JD_A - JD_D;
        
        %%% Starting state for this leg of the trajectory [r0; v0]
        X0_sc_leg2{xIndex, yIndex} = [rJupiter_D; vSc_JupiterDeparture];

    end
end


% ------------------------------------
%%% Pork Chop Plot
% ------------------------------------
figure('position',[440 234 762 564]); hold all
vInf_plu_contours = [12.5, 13, 13.5, 13.79, 14, 14.5, 15, 15.5];
[cs1, h1] = contour(JD_Departures_2-JD_Departures_2(1),JD_Arrivals_2-JD_Arrivals_2(1),VInfs_plu_m,vInf_plu_contours,'ShowText','on','color',colors.std.blue,'linewidth',1.2);

vInf_jup_contours = [17.5, 18, 18.465, 19, 19.5, 20, 20.5];
[cs2, h2] = contour(JD_Departures_2-JD_Departures_2(1),JD_Arrivals_2-JD_Arrivals_2(1),VInfs_jup_p,vInf_jup_contours,'ShowText','on','color',colors.std.ltred,'linewidth',1.2);

TOF_contours = [2700, 2800, 2900, 3000, 3058, 3100, 3200, 3300];
[cs3,h3] = contour(JD_Departures_2-JD_Departures_2(1),JD_Arrivals_2-JD_Arrivals_2(1),TOFs_days,TOF_contours,'ShowText','on','color',colors.std.black,'linewidth',1.2);

p1 = plot(JD_JGA_NH - JD_Departures_2(1), JD_Pluto_NH - JD_Arrivals_2(1),'o','markeredgecolor',colors.std.black,'markerfacecolor',colors.std.cyan,'markersize',10);

PlotBoi2('Days Past 29 January 2007','Days Past 17 September 2014',16)
grid off
clabel(cs1 ,h1,'color',colors.std.blue)
clabel(cs2 ,h2,'color',colors.std.ltred)
clabel(cs3 ,h3,'color',colors.std.black)

legend([h1, h2, h3, p1],'V_\infty @ Pluto, km/s','V_\infty @ Jupiter, km/s','TOF, days','New Horizons')

%% =======================================================================
%%% Searching for candidates
% ========================================================================
%%% Criteria
% -Arrival date matches departure date at flyby planet
% -Difference between vInf_m and vInf_p is less than some tolerance
% -Flyby radius is above some alitude tolerance
% -vInf and c3 requirements

%%% Method
% -Check every single possibility
% -Store trajectory if criteria met
% -

% ------------------------------------
%%% Set requirements
% ------------------------------------
%%% Requirements
% x-Must depart Earth on Jan 9, 2006
% x-Must depart Earth with C3 no larger than 180 km^2/s^2
% x-They may pass by Jupiter someime in the window 2454129.5 to 2454239.5
% x-Must depart Jupiter on same day as arrival
% x-|vInf_p| must be within 0.1 km/s of |vInf_m|. The closer the better
% x-Must encounter Pluto sometime in the window 2456917.5 to 2457517.5
% x-Value of |vInf_plu_m| must be no larger than 14.5 km/s
% -Must fly no closer to jupiter than 30 R (2144760 km) from the surface

%%% Lower and upper bounds on earth launch date
[ req_JD_Launch_lb ] = calendar2julian(2006, 1, 9, 0, 0, 0);
[ req_JD_Launch_ub ] = calendar2julian(2006, 1, 10, 0, 0, 0);
req_JD_Launch_ub = req_JD_Launch_ub - 1e-9;

%%% C3 requirement on Earth departure
req_c3_earth_p_ub = 180; % km^2/s^2

%%% Lower and uppder bounds on Jupiter arrival date
req_JD_jup_lb = 2454129.5;
req_JD_jup_ub = 2454239.5;

%%% Tolerance on abs(|vInf_m| - |vInf_p|) at Jupiter
req_tol_vInf_jup = 0.05; % km/s

%%% Lower and uppder bounds on Pluto arrival date
req_JD_plu_lb = 2456917.5;
req_JD_plu_ub = 2457517.5;

%%% |vInf_m| requirement at Pluto arrival
req_vInf_plu_m_ub = 14.5; % km/s

%%% Flyby radius for Jupiter
req_flybyRadius = 2144760 + 2144760/30; % km


% ------------------------------------
%%% Searching for good dates
% ------------------------------------
solutionIndex = 0;
JDs_Launch_JupArr_JupDep_PluArr{100000} = [];
for earthDepartureIndex = 1:nPoints
    for jupiterArrivalIndex = 1:nPoints
        %%% Grab current launch window
        JDi_earthDeparture = JD_Departures_1(earthDepartureIndex, jupiterArrivalIndex);
        JDi_jupiterArrival = JD_Arrivals_1(earthDepartureIndex, jupiterArrivalIndex);
        
        %%% Skip this date if it's not in the proper launch window
        if (JDi_earthDeparture < req_JD_Launch_lb) || (JDi_earthDeparture > req_JD_Launch_ub)
            continue
        end
        
        %%% Skip this date if earth-departure C3 is too high
        if c3s_earth_p(earthDepartureIndex, jupiterArrivalIndex) > req_c3_earth_p_ub
            continue
        end
        
        %%% Skip this date if it's not in the proper launch window
        if (JDi_jupiterArrival < req_JD_jup_lb) || (JDi_jupiterArrival > req_JD_jup_ub)
            continue
        end
                
        %%% check leg-2 of trajectory
        for jupiterDepartureIndex = 1:nPoints
            for plutoArrivalIndex = 1:nPoints
                %%% Grab current launch window
                JDi_jupiterDeparture = JD_Departures_2(jupiterDepartureIndex, plutoArrivalIndex);
                JDi_plutoArrival = JD_Arrivals_2(jupiterDepartureIndex, plutoArrivalIndex);
                
                %%% Make sure arrival and departure at Jupiter occur on the
                %%% same day
                if (JDi_jupiterDeparture ~= JDi_jupiterArrival)
                    continue
                end
                
                %%% Check magnitude of difference in vIns
                if abs(VInfs_jup_p(jupiterDepartureIndex,plutoArrivalIndex) - VInfs_jup_m(earthDepartureIndex,jupiterArrivalIndex)) > req_tol_vInf_jup
                    continue
                end
                
                %%% Skip this date if it's not in the proper launch window
                if (JDi_plutoArrival < req_JD_plu_lb) || (JDi_plutoArrival > req_JD_plu_ub)
                    continue
                end
                
                %%% Skip this date if |vInf_plu_m| is too large
                if VInfs_plu_m(jupiterDepartureIndex,plutoArrivalIndex) > req_vInf_plu_m_ub
                    continue
                end
                
                %%% Skip this date if flyby radius too close to jupiter
                if flybyRadii(jupiterDepartureIndex, plutoArrivalIndex) < req_flybyRadius
                    continue
                end
                
                % ------------------------------------
                %%% Store solution
                % ------------------------------------
                solutionIndex = solutionIndex + 1;
                
                %%% Relevant dates
                JDs_Launch_JupArr_JupDep_PluArr{solutionIndex}.JD_EarthDeparture   = JDi_earthDeparture;
                JDs_Launch_JupArr_JupDep_PluArr{solutionIndex}.JD_JupiterArrival   = JDi_jupiterArrival;
                JDs_Launch_JupArr_JupDep_PluArr{solutionIndex}.JD_JupiterDeparture = JDi_jupiterDeparture;
                JDs_Launch_JupArr_JupDep_PluArr{solutionIndex}.JD_PlutoArrival     = JDi_plutoArrival;
                JDs_Launch_JupArr_JupDep_PluArr{solutionIndex}.tof_days            = JDi_plutoArrival-JDi_earthDeparture;
                
                %%% Relevant values
                JDs_Launch_JupArr_JupDep_PluArr{solutionIndex}.c3_earth_p    = c3s_earth_p(earthDepartureIndex, jupiterArrivalIndex);
                JDs_Launch_JupArr_JupDep_PluArr{solutionIndex}.vInf_jup_m    = VInfs_jup_m(earthDepartureIndex,jupiterArrivalIndex);
                JDs_Launch_JupArr_JupDep_PluArr{solutionIndex}.vInf_jup_p    = VInfs_jup_p(jupiterDepartureIndex,plutoArrivalIndex);
                JDs_Launch_JupArr_JupDep_PluArr{solutionIndex}.vInf_jup_diff = abs(VInfs_jup_p(jupiterDepartureIndex,plutoArrivalIndex) - VInfs_jup_m(earthDepartureIndex,jupiterArrivalIndex));
                JDs_Launch_JupArr_JupDep_PluArr{solutionIndex}.vInf_plu_m    = VInfs_plu_m(jupiterDepartureIndex,plutoArrivalIndex);
                JDs_Launch_JupArr_JupDep_PluArr{solutionIndex}.leg1_X0       = X0_sc_leg1(earthDepartureIndex,jupiterArrivalIndex);
                JDs_Launch_JupArr_JupDep_PluArr{solutionIndex}.leg2_X0       = X0_sc_leg2(jupiterDepartureIndex,plutoArrivalIndex);

                
            end
        end
        
        
    end
end

% ------------------------------------
%%% Reassigning Data
% ------------------------------------
%%% Deleting empty cells
if length(JDs_Launch_JupArr_JupDep_PluArr(~cellfun('isempty',JDs_Launch_JupArr_JupDep_PluArr))) == length(JDs_Launch_JupArr_JupDep_PluArr)
    warning('Need to preallocate more space to JDs_Launch_JupArr_JupDep_PluArr')
end
JDs_Launch_JupArr_JupDep_PluArr = JDs_Launch_JupArr_JupDep_PluArr(~cellfun('isempty',JDs_Launch_JupArr_JupDep_PluArr));

%%% Reassigning results to get them out of cells
candidates = [JDs_Launch_JupArr_JupDep_PluArr{:}];

candidates_JD_EarthDeparture   = [candidates.JD_EarthDeparture];
candidates_JD_JupiterArrival   = [candidates.JD_JupiterArrival];
candidates_JD_JupiterDeparture = [candidates.JD_JupiterDeparture];
candidates_JD_PlutoArrival     = [candidates.JD_PlutoArrival];
candidates_tof_days      = [candidates.tof_days];
candidates_c3_earth_p    = [candidates.c3_earth_p];
candidates_vInf_jup_m    = [candidates.vInf_jup_m];
candidates_vInf_jup_p    = [candidates.vInf_jup_p];
candidates_vInf_jup_diff = [candidates.vInf_jup_diff];
candidates_vInf_plu_m    = [candidates.vInf_plu_m];

%%% My optimal solution .... minimize vInf_plu_m
myOptimalIndices = find(candidates_vInf_plu_m == min(candidates_vInf_plu_m));
candidates(myOptimalIndices(1))
candidates(myOptimalIndices(2))
candidates(myOptimalIndices(3))

myOptimalIndex = myOptimalIndices(1);

JDs_Launch_JupArr_JupDep_PluArr{myOptimalIndex}.leg1_X0
% candidates(myOptimalIndex)

% ------------------------------------
%%% Ephemeris positions of planets
% ------------------------------------
%%% Earth state at Launch
[rEarth0_leg1, vEarth0_leg1] = JulianDate_2_HeliocentricState(candidates_JD_EarthDeparture(myOptimalIndex), 'Earth');
X0_earth_leg1 = [rEarth0_leg1; vEarth0_leg1];

%%% Jupiter state at JGA
[rJupiter0_leg2, vJupiter0_leg2] = JulianDate_2_HeliocentricState(candidates_JD_JupiterArrival(myOptimalIndex), 'Jupiter');
X0_jupiter_leg2 = [rJupiter0_leg2; vJupiter0_leg2];

%%% Pluto state at flyby
[rPluto0_leg2, vPluto0_leg2] = JulianDate_2_HeliocentricState(candidates_JD_PlutoArrival(myOptimalIndex), 'Pluto');
X0_pluto_leg2 = [rPluto0_leg2; vPluto0_leg2];
% ------------------------------------
%%% Plotting chosen trajectory
% ------------------------------------
%%% Initial conditions
X0_sc_leg1 = JDs_Launch_JupArr_JupDep_PluArr{myOptimalIndex}.leg1_X0{:};
X0_sc_leg2 = JDs_Launch_JupArr_JupDep_PluArr{myOptimalIndex}.leg2_X0{:};

%%% Time
t_i_leg1 = (JDs_Launch_JupArr_JupDep_PluArr{myOptimalIndex}.JD_EarthDeparture)*86400; % sec
t_f_leg1 = (JDs_Launch_JupArr_JupDep_PluArr{myOptimalIndex}.JD_JupiterArrival)*86400; 

t_i_leg2 = (JDs_Launch_JupArr_JupDep_PluArr{myOptimalIndex}.JD_JupiterDeparture)*86400; % sec
t_f_leg2 = (JDs_Launch_JupArr_JupDep_PluArr{myOptimalIndex}.JD_PlutoArrival)*86400; 

n_dt = 1000;

time0_leg1 = linspace(t_i_leg1,t_f_leg1,n_dt);
time0_leg2 = linspace(t_i_leg2,t_f_leg2,n_dt);

%%% Choosing ode tolerance
tol = 1e-13;

%%% Setting integrator options
options = odeset('RelTol',tol,'AbsTol',tol);

%%% Propagating
[time_leg1, X_sc_leg1] = ode113(@Int_2BI, time0_leg1, X0_sc_leg1, options, Sun.u);
[time_leg1, X_earth_leg1] = ode113(@Int_2BI, time0_leg1, X0_earth_leg1, options, Sun.u);

[time_leg2, X_sc_leg2] = ode113(@Int_2BI, time0_leg2, X0_sc_leg2, options, Sun.u);
[time_leg2, X_jupiter_leg2] = ode113(@Int_2BI, time0_leg2.*100, X0_jupiter_leg2, options, Sun.u);
[time_leg2, X_pluto_leg2] = ode113(@Int_2BI, time0_leg2.*100, X0_pluto_leg2, options, Sun.u);

figure; hold all
plot3(X_sc_leg1(:,1),X_sc_leg1(:,2),X_sc_leg1(:,3),'k')
plot3(X_earth_leg1(:,1),X_earth_leg1(:,2),X_earth_leg1(:,3),'b')
p1 = plot3(X_earth_leg1(1,1),X_earth_leg1(1,2),X_earth_leg1(1,3),'ko','markerfacecolor',colors.std.blue);

plot3(X_sc_leg2(:,1),X_sc_leg2(:,2),X_sc_leg2(:,3),'k')
plot3(X_jupiter_leg2(:,1),X_jupiter_leg2(:,2),X_jupiter_leg2(:,3),'r')
p2 = plot3(X0_jupiter_leg2(1),X0_jupiter_leg2(2),X0_jupiter_leg2(3),'ko','markerfacecolor',colors.std.red);
plot3(X_pluto_leg2(:,1),X_pluto_leg2(:,2),X_pluto_leg2(:,3),'m')
p3 = plot3(X0_pluto_leg2(1),X0_pluto_leg2(2),X0_pluto_leg2(3),'ko','markerfacecolor',colors.std.mag);

PlotBoi3('$X_n$','$Y_n$','$Z_n$',18,'LaTex')
axis equal
legend([p1 p2 p3], 'Earth','Jupiter','Pluto')









































