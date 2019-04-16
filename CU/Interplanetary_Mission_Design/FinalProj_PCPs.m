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
run_PCP_Launch_2_VGA = 1;
run_PCP_VGA_2_EGA1   = 1;
run_PCP_EGA2_2_SOI   = 1;
run_resonanceStudy   = 1;

%% =======================================================================
%%% Setup
% ========================================================================
% ------------------------------------
%%% Bodies
% ------------------------------------
Sun    = bodies.sun;
Venus  = bodies.venus;
Earth  = bodies.earth;
Saturn = bodies.saturn;

% ------------------------------------
%%% Parameters
% ------------------------------------
au_km = 1.49597870700e8; % km

Sun.u    = 1.32712440018e11; % km^3/s^2
Venus.u  = 3.24858599e5;     % km^3/s^2
Earth.u  = 3.98600433e5;     % km^3/s^2
Saturn.u = 3.7931208e7;      % km^3/s^2

Venus.R  = 6051.8;  % km
Earth.R  = 6378.14; % km
Saturn.R = 60268;   % km

% ------------------------------------
%%% Julian Date ranges
% ------------------------------------
% [ JD ] = calendar2julian(yr,  mo, d, h, min, s)
% [Year, Mon, Day, hr, min, s] = julian2calendar(JD)
% 

%%% Imporant dates
JD_2025 = calendar2julian(2025,  1, 1, 0, 0, 0);
JD_2041 = calendar2julian(2041,  1, 1, 0, 0, 0);

%%% Resonance time
days_per_year = 365.242189; % Note: Used in creating JD_Arrivals_SOI_vec
seconds_per_year = days_per_year*86400;
resonance_period_days = 3*days_per_year; % 2 years
resonance_period_secs = resonance_period_days*86400; % 2 years


%%% Ranges
nPoints = 75;
JDs_Launch = linspace(JD_2025,JD_2025 + days_per_year, nPoints);
JDs_VGA    = linspace(2460751.5, 2460956.5, nPoints);
JDs_EGA1   = linspace(2460957.5, 2461456.5, nPoints);
JDs_EGA2   = JDs_EGA1 + resonance_period_days;
JDs_SOI    = linspace(2463060.737444906, JD_2041, nPoints);


% ------------------------------------
%%% Julian Date values
% ------------------------------------
JD_Launch = 2.460740664168338e+06;
JD_VGA    = 2.460909405405405e+06;
JD_EGA1   = 2.461227229729730e+06;
JD_EGA2   = 2.462322956296730e+06;
JD_SOI    = 2.463995808405742e+06;

%% =======================================================================
%%% Pork Chop Plot - Launch to VGA (leg1)
% ========================================================================
if run_PCP_Launch_2_VGA == 1
% ------------------------------------
%%% Setting up julian date vectors
% ------------------------------------
JD_Departures_Launch_vec = linspace(JD_2025, JD_2025+125, nPoints);
JD_Arrivals_VGA_vec     = linspace(JD_2025 + 75, JD_2025+280, nPoints);

[JD_Departures_Launch, JD_Arrivals_VGA] = meshgrid(JD_Departures_Launch_vec, JD_Arrivals_VGA_vec);

% ------------------------------------
%%% Looping through options and computing vInf and c3
% ------------------------------------
VInfs_VGA_m = zeros(size(JD_Departures_Launch));
c3s_Launch_p = zeros(size(JD_Departures_Launch));
TOFs_days_leg1 = zeros(size(JD_Departures_Launch));

for depIndex = 1:nPoints
    for arrIndex = 1:nPoints
        %%% Grab current launch window
        JD_Dep = JD_Departures_Launch(depIndex, arrIndex);
        JD_Arr = JD_Arrivals_VGA(depIndex, arrIndex);
        
        %%% Find Earth position at departure time
        [rEarth_D, vEarth_D] = JulianDate_2_HeliocentricState(JD_Dep, 'Earth');
        
        %%% Find Venus position at arrival time
        [rVenus_A, vVenus_A] = JulianDate_2_HeliocentricState(JD_Arr, 'Venus');
        
        %%% Use Lambert solver to find solution
        tof_sec = (JD_Arr - JD_Dep)*86400;
        [vSc_EarthDeparture, vSc_VenusArrival,~] = lambertSolver(rEarth_D,rVenus_A,tof_sec,0,0,Sun.u);
        
        %%% C3 at Earth departure
        vInf_Earth_p = vSc_EarthDeparture - vEarth_D;
        c3s_Launch_p(depIndex, arrIndex) = norm(vInf_Earth_p)^2;
        
        %%% vInf at Venus arrival
        vInf_Venus_m = vSc_VenusArrival - vVenus_A;
        VInfs_VGA_m(depIndex, arrIndex) = norm(vInf_Venus_m);
        
        %%% TOFs
        TOFs_days_leg1(depIndex, arrIndex) = JD_Arr - JD_Dep;

    end
end


% ------------------------------------
%%% Pork Chop Plot
% ------------------------------------
figure('position',[74 233 762 564]); hold all
c3_contours = [5,7.2,8, 10, 15, 20, 25];
[cs1, h1] = contour(JD_Departures_Launch-JD_Departures_Launch(1),JD_Arrivals_VGA-JD_Arrivals_VGA(1),c3s_Launch_p,c3_contours,'ShowText','on','color',colors.std.ltred,'linewidth',1.2);

vInf_contours = [3.8, 4, 5, 8, 12, 15, 17, 20];
[cs2, h2] = contour(JD_Departures_Launch-JD_Departures_Launch(1),JD_Arrivals_VGA-JD_Arrivals_VGA(1),VInfs_VGA_m,vInf_contours,'ShowText','on','color',colors.std.blue,'linewidth',1.2);

TOF_contours_leg1 = [80, 100, 120, 140, 160, 200];
[cs3,h3] = contour(JD_Departures_Launch-JD_Departures_Launch(1),JD_Arrivals_VGA-JD_Arrivals_VGA(1),TOFs_days_leg1,'ShowText','on','color',colors.std.black,'linewidth',1.2);

p1 = plot(JD_Launch - JD_Departures_Launch_vec(1), JD_VGA - JD_Arrivals_VGA_vec(1),'o','markersize',10,'markeredgecolor',colors.std.black,'markerfacecolor',colors.std.cyan);

[yr, mo, d, ~, ~, ~] = julian2calendar(JD_Departures_Launch_vec(1));
calendar0_launch = sprintf('Launch, Days Past %d/%d/%d',mo,d,yr);

[yr, mo, d, ~, ~, ~] = julian2calendar(JD_Arrivals_VGA_vec(1));
calendar0_VGA = sprintf('VGA, Days Past %d/%d/%d',mo,d,yr);

PlotBoi2(calendar0_launch,calendar0_VGA,16)
grid off
clabel(cs1 ,h1,'color',colors.std.ltred)
clabel(cs2 ,h2,'color',colors.std.blue)
clabel(cs3 ,h3,'color',colors.std.black)

legend([h1, h2, h3, p1],'C_3 @ Launch, km^2/s^2','V_\infty @ VGA, km/s','TOF, days','My Dates')

end % run_PCP_Launch_2_VGA == 1


%% =======================================================================
%%% Pork Chop Plot - VGA to EGA1 (leg2)
% ========================================================================
if run_PCP_VGA_2_EGA1 == 1
% ------------------------------------
%%% Setting up julian dates
% ------------------------------------
JD_Departures_VGA_vec = JD_Arrivals_VGA_vec;
JD_Arrivals_EGA1_vec  = linspace(JD_Departures_VGA_vec(end)+1, JD_Departures_VGA_vec(end)+500, nPoints);


[JD_Departures_VGA, JD_Arrivals_EGA1] = meshgrid(JD_Departures_VGA_vec, JD_Arrivals_EGA1_vec);

% ------------------------------------
%%% Looping through options and computing vInf and c3
% ------------------------------------
vInfs_EGA1_m{size(JD_Departures_VGA,1),size(JD_Departures_VGA,2)} = [];
VInfs_EGA1_m = zeros(size(JD_Departures_VGA));
VInfs_VGA_p = zeros(size(JD_Departures_VGA));
TOFs_days_leg2 = zeros(size(JD_Departures_VGA));

for depIndex = 1:nPoints
    for arrIndex = 1:nPoints
        %%% Grab current launch window
        JD_Dep = JD_Departures_VGA(depIndex, arrIndex);
        JD_Arr = JD_Arrivals_EGA1(depIndex, arrIndex);
        
        %%% Find Venus position at departure time
        [rVenus_D, vVenus_D] = JulianDate_2_HeliocentricState(JD_Dep, 'Venus');
        
        %%% Find Earth position at arrival time
        [rEarth_A, vEarth_A] = JulianDate_2_HeliocentricState(JD_Arr, 'Earth');
        
        %%% Use Lambert solver to find solution
        tof_leg2_sec = (JD_Arr - JD_Dep)*86400;
        [vSc_VGA_p, vSc_EGA1_m,~] = lambertSolver(rVenus_D,rEarth_A,tof_leg2_sec,0,0,Sun.u);
        
        %%% vInf at Venus departure
        vInf_Venus_p = vSc_VGA_p - vVenus_D;
        VInfs_VGA_p(depIndex, arrIndex) = norm(vInf_Venus_p);
        
        %%% vInf at Earth arrival
        vInf_EGA1_m = vSc_EGA1_m - vEarth_A;
        vInfs_EGA1_m{depIndex, arrIndex} = vInf_EGA1_m;
        VInfs_EGA1_m(depIndex, arrIndex) = norm(vInf_EGA1_m);
        
        %%% TOFs
        TOFs_days_leg2(depIndex, arrIndex) = JD_Arr - JD_Dep;

    end
end

% ------------------------------------
%%% Pork Chop Plot
% ------------------------------------
figure('position',[440 234 762 564]); hold all
vInf_VGA_contours = [4, 5, 6, 8, 12, 25, 35];
[cs1, h1] = contour(JD_Departures_VGA-JD_Departures_VGA(1),JD_Arrivals_EGA1-JD_Arrivals_EGA1(1),VInfs_VGA_p,vInf_VGA_contours,'ShowText','on','color',colors.std.blue,'linewidth',1.2);

vInf_EGA1_contours = [7.5, 8, 9, 10, 12 , 16, 35, 60];
[cs2, h2] = contour(JD_Departures_VGA-JD_Departures_VGA(1),JD_Arrivals_EGA1-JD_Arrivals_EGA1(1),VInfs_EGA1_m,vInf_EGA1_contours,'ShowText','on','color',colors.std.ltred,'linewidth',1.2);

TOF_contours_leg2 = [260, 280, 300, 320, 340];
[cs3,h3] = contour(JD_Departures_VGA-JD_Departures_VGA(1),JD_Arrivals_EGA1-JD_Arrivals_EGA1(1),TOFs_days_leg2,'ShowText','on','color',colors.std.black,'linewidth',1.2);

p1 = plot(JD_VGA - JD_Departures_VGA_vec(1), JD_EGA1 - JD_Arrivals_EGA1_vec(1),'o','markersize',10,'markeredgecolor',colors.std.black,'markerfacecolor',colors.std.cyan);

[yr, mo, d, ~, ~, ~] = julian2calendar(JD_Departures_VGA_vec(1));
calendar0_VGA = sprintf('VGA, Days Past %d/%d/%d',mo,d,yr);

[yr, mo, d, ~, ~, ~] = julian2calendar(JD_Arrivals_EGA1_vec(1));
calendar0_EGA1 = sprintf('EGA1, Days Past %d/%d/%d',mo,d,yr);

PlotBoi2(calendar0_VGA,calendar0_EGA1,16)
grid off
clabel(cs1 ,h1,'color',colors.std.blue)
clabel(cs2 ,h2,'color',colors.std.ltred)
clabel(cs3 ,h3,'color',colors.std.black)

legend([h1, h2, h3, p1],'V_\infty @ VGA, km/s','V_\infty @ EGA1, km/s','TOF, days', 'My Dates')

end % run_PCP_VGA_2_EGA1 == 1


%% =======================================================================
%%%  Pork Chop Plot - EGA2 to SOI (leg2)
% ========================================================================
if run_PCP_EGA2_2_SOI == 1
% ------------------------------------
%%% Setting up julian dates
% ------------------------------------
JD_Departures_EGA2_vec = JD_Arrivals_EGA1_vec + resonance_period_days;
JD_Arrivals_SOI_vec = linspace(JD_Departures_EGA2_vec(1) + 2.5*days_per_year, JD_2041-2*365, nPoints);

[JD_Departures_EGA2, JD_Arrivals_SOI] = meshgrid(JD_Departures_EGA2_vec, JD_Arrivals_SOI_vec);

% ------------------------------------
%%% Looping through options and computing vInf and c3
% ------------------------------------
VInfs_SOI_m = zeros(size(JD_Departures_EGA2));
VInfs_EGA2_p = zeros(size(JD_Departures_EGA2));
TOFs_days_leg4 = zeros(size(JD_Departures_EGA2));

for depIndex = 1:nPoints
    for arrIndex = 1:nPoints
        %%% Grab current launch window
        JD_Dep = JD_Departures_EGA2(depIndex, arrIndex);
        JD_Arr = JD_Arrivals_SOI(depIndex, arrIndex);
        
        %%% Find Earth position at departure time
        [rEarth_D, vEarth_D] = JulianDate_2_HeliocentricState(JD_Dep, 'Earth');
        
        %%% Find Saturn position at arrival time
        [rSaturn_A, vSaturn_A] = JulianDate_2_HeliocentricState(JD_Arr, 'Saturn');
        
        %%% Use Lambert solver to find solution
        tof_leg4_sec = (JD_Arr - JD_Dep)*86400;
        [vSc_EGA2_p, vSc_SOI_m,~] = lambertSolver(rEarth_D,rSaturn_A,tof_leg4_sec,0,0,Sun.u);
        
        %%% vInf at Earth departure
        vInf_Earth_p = vSc_EGA2_p - vEarth_D;
        VInfs_EGA2_p(depIndex, arrIndex) = norm(vInf_Earth_p);

        %%% vInf at Saturn arrival
        vInf_Saturn = vSc_SOI_m - vSaturn_A;
        VInfs_SOI_m(depIndex, arrIndex) = norm(vInf_Saturn);
        
        %%% TOFs
        TOFs_days_leg4(depIndex, arrIndex) = JD_Arr - JD_Dep;
    end
end

% ------------------------------------
%%% Pork Chop Plot
% ------------------------------------
figure('position',[440 234 762 564]); hold all
vInf_EGA2_contours = [10, 11, 13, 15, 20, 25, 30, 40, 50];
[cs1, h1] = contour(JD_Departures_EGA2-JD_Departures_EGA2(1),JD_Arrivals_SOI-JD_Arrivals_SOI(1),VInfs_EGA2_p,vInf_EGA2_contours,'ShowText','on','color',colors.std.blue,'linewidth',1.2);

vInf_SOI_contours = [1, 2, 3, 4, 5,5.8, 6, 7, 8, 9];
[cs2, h2] = contour(JD_Departures_EGA2-JD_Departures_EGA2(1),JD_Arrivals_SOI-JD_Arrivals_SOI(1),VInfs_SOI_m,vInf_SOI_contours,'ShowText','on','color',colors.std.ltred,'linewidth',1.2);

TOF_contours_leg4 = [260, 280, 300, 320, 340];
[cs3,h3] = contour(JD_Departures_EGA2-JD_Departures_EGA2(1),JD_Arrivals_SOI-JD_Arrivals_SOI(1),TOFs_days_leg4,'ShowText','on','color',colors.std.black,'linewidth',1.2);

p1 = plot(JD_EGA2 - JD_Departures_EGA2_vec(1), JD_SOI - JD_Arrivals_SOI_vec(1),'o','markersize',10,'markeredgecolor',colors.std.black,'markerfacecolor',colors.std.cyan);

[yr, mo, d, ~, ~, ~] = julian2calendar(JD_Departures_EGA2_vec(1));
calendar0_EGA2 = sprintf('EGA2, Days Past %d/%d/%d',mo,d,yr);

[yr, mo, d, ~, ~, ~] = julian2calendar(JD_Arrivals_SOI_vec(1));
calendar0_SOI = sprintf('SOI, Days Past %d/%d/%d',mo,d,yr);

PlotBoi2(calendar0_EGA2,calendar0_SOI,16)
grid off
clabel(cs1 ,h1,'color',colors.std.blue)
clabel(cs2 ,h2,'color',colors.std.ltred)
clabel(cs3 ,h3,'color',colors.std.black)

legend([h1, h2, h3, p1],'V_\infty @ EGA2, km/s','V_\infty @ SOI, km/s','TOF, days','My Dates')

end % run_PCP_EGA2_2_SOI == 1




