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
run_problem = 4;

% ========================================================================
%%% Question 2
% ========================================================================
if run_problem == 2
% ------------------------------------
%%% Setting up julian dates
% ------------------------------------
nPoints = 200;
JD_Dep = linspace(2461295, 2461415, nPoints);
JD_Arr   = linspace(2461530, 2461710, nPoints);

[JD_Departures, JD_Arrivals] = meshgrid(JD_Dep, JD_Arr);

% ------------------------------------
%%% Looping through options and computing vInf and c3
% ------------------------------------
vInfs     = zeros(size(JD_Departures));
c3s       = zeros(size(JD_Departures));
TOFs_days = zeros(size(JD_Departures));

for xIndex = 1:nPoints
    for yIndex = 1:nPoints
        %%% Grab current launch window
        JD_D = JD_Departures(xIndex, yIndex);
        JD_A = JD_Arrivals(xIndex, yIndex);
        
        %%% Find Earth position at departure time
        [~, a, e, i_rad, Omega_rad, ~, w_rad, ~, ta_rad] = ...
            getPlanetElements_Meeus(JD_D, 'Earth', 'radians');
        [rEarth_D, vEarth_D] = OE2ECI(a, e, i_rad, Omega_rad, w_rad, ta_rad, bodies.sun.u);
        
        %%% Find Earth position at departure time
        [~, a, e, i_rad, Omega_rad, ~, w_rad, ~, ta_rad] = ...
            getPlanetElements_Meeus(JD_A, 'Mars', 'radians');
        [rMars_A, vMars_A] = OE2ECI(a, e, i_rad, Omega_rad, w_rad, ta_rad, bodies.sun.u);
        
        %%% Use Lambert solver to find solution
        tof_sec = (JD_A - JD_D)*86400;
        [v_EarthDeparture, v_MarsArrival,~] = lambertSolver(rEarth_D,rMars_A,tof_sec,0,0,bodies.sun.u);
        
        %%% C3 at Earth departure
        vInf_Earth = v_EarthDeparture - vEarth_D;
        c3s(xIndex, yIndex) = norm(vInf_Earth)^2;
        
        %%% vInf at Mars arrival
        vInf_Mars = v_MarsArrival - vMars_A;
        vInfs(xIndex, yIndex) = norm(vInf_Mars);
        
        %%% TOFs
        TOFs_days(xIndex, yIndex) = JD_A - JD_D;

    end
end

% ------------------------------------
%%% Pork Chop Plot
% ------------------------------------
figure('position',[440 234 762 564]); hold all
vInf_contours = [8, 7, 6, 5, 4, 3.5, 3, 2.8, 2.6];
[cs1, h1] = contour(JD_Departures-JD_Departures(1),JD_Arrivals-JD_Arrivals(1),vInfs,vInf_contours,'ShowText','on','color',colors.std.blue,'linewidth',1.2);

% figure
c3_contours = [70, 50, 30, 20, 15, 13, 12.0, 10, 9.3, 9.16];
[cs2, h2] = contour(JD_Departures-JD_Departures(1),JD_Arrivals-JD_Arrivals(1),c3s,c3_contours,'ShowText','on','color',colors.std.ltred,'linewidth',1.2);

% figure
TOF_contours = [150, 200,225, 250, 300, 350, 400];
[cs3,h3] = contour(JD_Departures-JD_Departures(1),JD_Arrivals-JD_Arrivals(1),TOFs_days,TOF_contours,'ShowText','on','color',colors.std.black,'linewidth',1.2);

PlotBoi2('Days Past 11 September 2026','Days Past 09 January 2027',16)
grid off
clabel(cs1 ,h1,'color',colors.std.blue)
clabel(cs2 ,h2,'color',colors.std.ltred)
clabel(cs3 ,h3,'color',colors.std.black)

legend([h1, h2, h3],'V_\infty @ Mars, km/s','C_3, km^2/s^2','TOF, days')

end










% ========================================================================
%%% Question 3
% ========================================================================
if run_problem == 3
% ------------------------------------
%%% Setting up julian dates
% ------------------------------------
nPoints = 200;
JD_Dep = linspace(2462075, 2462175 , nPoints);
JD_Arr   = linspace(2462280, 2462460, nPoints);

[JD_Departures, JD_Arrivals] = meshgrid(JD_Dep, JD_Arr);

% ------------------------------------
%%% Looping through options and computing vInf and c3
% ------------------------------------
vInfs     = zeros(size(JD_Departures));
c3s       = zeros(size(JD_Departures));
TOFs_days = zeros(size(JD_Departures));

for xIndex = 1:nPoints
    for yIndex = 1:nPoints
        %%% Grab current launch window
        JD_D = JD_Departures(xIndex, yIndex);
        JD_A = JD_Arrivals(xIndex, yIndex);
        
        %%% Find Earth position at departure time
        [~, a, e, i_rad, Omega_rad, ~, w_rad, ~, ta_rad] = ...
            getPlanetElements_Meeus(JD_D, 'Earth', 'radians');
        [rEarth_D, vEarth_D] = OE2ECI(a, e, i_rad, Omega_rad, w_rad, ta_rad, bodies.sun.u);
        
        %%% Find Earth position at departure time
        [~, a, e, i_rad, Omega_rad, ~, w_rad, ~, ta_rad] = ...
            getPlanetElements_Meeus(JD_A, 'Mars', 'radians');
        [rMars_A, vMars_A] = OE2ECI(a, e, i_rad, Omega_rad, w_rad, ta_rad, bodies.sun.u);
        
        %%% Use Lambert solver to find solution
        tof_sec = (JD_A - JD_D)*86400;
        [v_EarthDeparture, v_MarsArrival,~] = lambertSolver(rEarth_D,rMars_A,tof_sec,0,0,bodies.sun.u);
        
        %%% C3 at Earth departure
        vInf_Earth = v_EarthDeparture - vEarth_D;
        c3s(xIndex, yIndex) = norm(vInf_Earth)^2;
        
        %%% vInf at Mars arrival
        vInf_Mars = v_MarsArrival - vMars_A;
        vInfs(xIndex, yIndex) = norm(vInf_Mars);
        
        %%% TOFs
        TOFs_days(xIndex, yIndex) = JD_A - JD_D;

    end
end
% ------------------------------------
%%% Pork Chop Plot
% ------------------------------------
figure('position',[440 234 762 564]); hold all
vInf_contours = [8, 7, 6, 5, 4.5, 4, 3.7, 3.6, 3.5, 3, 2.97];
[cs1, h1] = contour(JD_Departures-JD_Departures(1),JD_Arrivals-JD_Arrivals(1),vInfs,vInf_contours,'ShowText','on','color',colors.std.blue,'linewidth',1.2);

% figure
c3_contours = [70, 50, 30, 20, 15, 13, 10, 9.3, 9];
[cs2, h2] = contour(JD_Departures-JD_Departures(1),JD_Arrivals-JD_Arrivals(1),c3s,c3_contours,'ShowText','on','color',colors.std.ltred,'linewidth',1.2);

% figure
TOF_contours = [150, 200, 250, 300, 350, 400];
[cs3,h3] = contour(JD_Departures-JD_Departures(1),JD_Arrivals-JD_Arrivals(1),TOFs_days,TOF_contours,'ShowText','on','color',colors.std.black,'linewidth',1.2);

PlotBoi2('Days Past 30 October 2028','Days Past 23 May 2029',16)
grid off
clabel(cs1 ,h1,'color',colors.std.blue)
clabel(cs2 ,h2,'color',colors.std.ltred)
clabel(cs3 ,h3,'color',colors.std.black)

legend([h1, h2, h3],'V_\infty @ Mars, km/s','C_3, km^2/s^2','TOF, days')


end










% ========================================================================
%%% Question 4
% ========================================================================
if run_problem == 4
% ------------------------------------
%%% Setting up julian dates
% ------------------------------------
nPoints = 500;
JD_2020 = calendar2julian(2020,1,1,0,0,0);
JD_Dep = linspace(JD_2020+185, JD_2020+265 , nPoints);
JD_Arr   = linspace(JD_2020+370, JD_2020+500, nPoints);

[JD_Departures, JD_Arrivals] = meshgrid(JD_Dep, JD_Arr);

% ------------------------------------
%%% Looping through options and computing vInf and c3
% ------------------------------------
vInfs     = zeros(size(JD_Departures));
c3s       = zeros(size(JD_Departures));
TOFs_days = zeros(size(JD_Departures));

for xIndex = 1:nPoints
    for yIndex = 1:nPoints
        %%% Grab current launch window
        JD_D = JD_Departures(xIndex, yIndex);
        JD_A = JD_Arrivals(xIndex, yIndex);
        
        %%% Find Earth position at departure time
        [~, a, e, i_rad, Omega_rad, ~, w_rad, ~, ta_rad] = ...
            getPlanetElements_Meeus(JD_D, 'Earth', 'radians');
        [rEarth_D, vEarth_D] = OE2ECI(a, e, i_rad, Omega_rad, w_rad, ta_rad, bodies.sun.u);
        
        %%% Find Earth position at departure time
        [~, a, e, i_rad, Omega_rad, ~, w_rad, ~, ta_rad] = ...
            getPlanetElements_Meeus(JD_A, 'Mars', 'radians');
        [rMars_A, vMars_A] = OE2ECI(a, e, i_rad, Omega_rad, w_rad, ta_rad, bodies.sun.u);
        
        %%% Use Lambert solver to find solution
        tof_sec = (JD_A - JD_D)*86400;
        [v_EarthDeparture, v_MarsArrival,~] = lambertSolver(rEarth_D,rMars_A,tof_sec,0,0,bodies.sun.u);
        
        %%% C3 at Earth departure
        vInf_Earth = v_EarthDeparture - vEarth_D;
        c3s(xIndex, yIndex) = norm(vInf_Earth)^2;
        
        %%% vInf at Mars arrival
        vInf_Mars = v_MarsArrival - vMars_A;
        vInfs(xIndex, yIndex) = norm(vInf_Mars);
        
        %%% TOFs
        TOFs_days(xIndex, yIndex) = JD_A - JD_D;

    end
end
% ------------------------------------
%%% Pork Chop Plot
% ------------------------------------
figure('position',[440 234 762 564]); hold all
% vInf_contours = linspace(0,2.8,100);
vInf_contours = [2.44 , 2.46, 2.52, 2.58, 2.64, 2.7, 2.75, 2.8];
[cs1, h1] = contour(JD_Departures-JD_Departures(1),JD_Arrivals-JD_Arrivals(1),vInfs,vInf_contours,'ShowText','on','color',colors.std.blue,'linewidth',1.2);

% figure
% c3_contours = linspace(0,13.5,100);
c3_contours = [13.1, 13.2, 13.3, 13.4, 13.5];
[cs2, h2] = contour(JD_Departures-JD_Departures(1),JD_Arrivals-JD_Arrivals(1),c3s,c3_contours,'ShowText','on','color',colors.std.ltred,'linewidth',1.2);

% figure
TOF_contours = [50, 75, 100, 125, 150, 175, 200];
[cs3,h3] = contour(JD_Departures-JD_Departures(1),JD_Arrivals-JD_Arrivals(1),TOFs_days,TOF_contours,'ShowText','on','color',colors.std.black,'linewidth',1.2);

PlotBoi2('Days Past 04 July 2020','Days Past 05 January 2021',16)
grid off
clabel(cs1 ,h1,'color',colors.std.blue)
clabel(cs2 ,h2,'color',colors.std.ltred)
clabel(cs3 ,h3,'color',colors.std.black)

legend([h1, h2, h3],'V_\infty @ Mars, km/s','C_3, km^2/s^2','TOF, days')




figure('position',[440 234 762 564]); hold all
% vInf_contours = linspace(0,2.8,100);
vInf_contours2 = [2.44 , 2.46, 2.52, 2.58, 2.64, 2.7, 2.75, 2.8];
[cs12, h12] = contour(JD_Departures-JD_Departures(1),JD_Arrivals-JD_Arrivals(1),vInfs,vInf_contours2,'ShowText','on','color',colors.std.blue,'linewidth',1.2);

% figure
% c3_contours = linspace(0,13.5,100);
c3_contours2 = [13.1, 13.2, 13.3, 13.4, 13.5];
[cs22, h22] = contour(JD_Departures-JD_Departures(1),JD_Arrivals-JD_Arrivals(1),c3s,c3_contours2,'ShowText','on','color',colors.std.ltred,'linewidth',1.2);

% figure
TOF_contours2 = [50, 75, 100, 125, 150, 175, 200];
[cs32,h32] = contour(JD_Departures-JD_Departures(1),JD_Arrivals-JD_Arrivals(1),TOFs_days,TOF_contours2,'ShowText','on','color',colors.std.black,'linewidth',1.2);

PlotBoi2('Days Past 04 July 2020','Days Past 05 January 2021',16)
grid off
clabel(cs12 ,h12,'color',colors.std.blue)
clabel(cs22 ,h22,'color',colors.std.ltred)
clabel(cs32 ,h32,'color',colors.std.black)

legend([h12, h22, h32],'V_\infty @ Mars, km/s','C_3, km^2/s^2','TOF, days')

xlim([10 25])
ylim([20 40])

end























toc



















