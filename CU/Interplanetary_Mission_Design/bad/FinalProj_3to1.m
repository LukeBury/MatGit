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
run_study_Launch_2_VGA                       = 1;
run_study_VGA_2_EGA1                         = 1;
run_candidateSearch_prelim                   = 1;
run_study_EGA2_2_SOI                         = 1;
run_candidateSearch_mid                      = 1;
run_resonanceStudy_and_candidateSearch_final = 1;

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
%%% Julian Date guesses for help in making PCPs
% ------------------------------------
% [ JD ] = calendar2julian(yr,  mo, d, h, min, s)
% [Year, Mon, Day, hr, min, s] = julian2calendar(JD)
% 

JD_2025 = 2460676.5;
JD_2041 = calendar2julian(2041,  1, 1, 0, 0, 0);

JD_guess_VGA = JD_2025 + 125;

JD_guess_EGA1 = JD_guess_VGA + 300;

%%% Resonance time
year_Earth_day = 365.242189; % Note: Used in creating JD_Arrivals_SOI_vec
year_Earth_sec = year_Earth_day*86400;
tof_resonance_day = 3*year_Earth_day; % 2 years
tof_resonance_sec = tof_resonance_day*86400; % 2 years

%% =======================================================================
%%% Pork Chop Plot - Launch to VGA (leg1)
% ========================================================================
if run_study_Launch_2_VGA == 1
% ------------------------------------
%%% Setting up julian dates
% ------------------------------------
% nPoints = 175;
989
% nPoints = 300;
nPoints = 75;
% JD_Departures_Launch_vec = linspace(JD_2025-100, JD_2025+125, nPoints);
JD_Departures_Launch_vec = linspace(JD_2025, JD_2025+125, nPoints);
JD_Arrivals_VGA_vec     = linspace(JD_guess_VGA-50, JD_guess_VGA+155, nPoints);
989
% JD_Departures_Launch_vec = linspace(JD_2025-100, JD_2025+465, nPoints);
% JD_Arrivals_VGA_vec     = linspace(JD_guess_VGA-100, JD_guess_VGA+455, nPoints);

[JD_Departures_Launch, JD_Arrivals_VGA] = meshgrid(JD_Departures_Launch_vec, JD_Arrivals_VGA_vec);

% ------------------------------------
%%% Looping through options and computing vInf and c3
% ------------------------------------
X0s_sc_leg1{size(JD_Departures_Launch,1),size(JD_Departures_Launch,2)} = [];
vInfs_VGA_m{size(JD_Departures_Launch,1),size(JD_Departures_Launch,2)} = [];
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
        vInfs_VGA_m{depIndex, arrIndex} = vInf_Venus_m;
        VInfs_VGA_m(depIndex, arrIndex) = norm(vInf_Venus_m);
        
        %%% TOFs
        TOFs_days_leg1(depIndex, arrIndex) = JD_Arr - JD_Dep;
        
        %%% Starting state for this leg of the trajectory [r0; v0]
        X0s_sc_leg1{depIndex, arrIndex} = [rEarth_D; vSc_EarthDeparture];

    end
end


% ------------------------------------
%%% Pork Chop Plot
% ------------------------------------
figure('position',[74 233 762 564]); hold all
c3_contours = [5,7.2,8, 10, 15, 20, 25];
[cs1, h1] = contour(JD_Departures_Launch-JD_Departures_Launch(1),JD_Arrivals_VGA-JD_Arrivals_VGA(1),c3s_Launch_p,c3_contours,'ShowText','on','color',colors.std.ltred,'linewidth',1.2);
% [cs1, h1] = contour(JD_Departures_Launch-JD_Departures_Launch(1),JD_Arrivals_VGA-JD_Arrivals_VGA(1),c3s_Launch_p,'ShowText','on','color',colors.std.ltred,'linewidth',1.2);

vInf_contours = [3.8, 4, 5, 8, 12, 15, 17, 20];
[cs2, h2] = contour(JD_Departures_Launch-JD_Departures_Launch(1),JD_Arrivals_VGA-JD_Arrivals_VGA(1),VInfs_VGA_m,vInf_contours,'ShowText','on','color',colors.std.blue,'linewidth',1.2);
% [cs2, h2] = contour(JD_Departures_Launch-JD_Departures_Launch(1),JD_Arrivals_VGA-JD_Arrivals_VGA(1),VInfs_VGA_m,'ShowText','on','color',colors.std.blue,'linewidth',1.2);

TOF_contours_leg1 = [80, 100, 120, 140, 160, 200];
% [cs3,h3] = contour(JD_Departures_Launch-JD_Departures_Launch(1),JD_Arrivals_VGA-JD_Arrivals_VGA(1),TOFs_days,TOF_contours,'ShowText','on','color',colors.std.black,'linewidth',1.2);
[cs3,h3] = contour(JD_Departures_Launch-JD_Departures_Launch(1),JD_Arrivals_VGA-JD_Arrivals_VGA(1),TOFs_days_leg1,'ShowText','on','color',colors.std.black,'linewidth',1.2);

[yr, mo, d, ~, ~, ~] = julian2calendar(JD_Departures_Launch_vec(1));
calendar0_launch = sprintf('Launch, Days Past %d/%d/%d',mo,d,yr);

[yr, mo, d, ~, ~, ~] = julian2calendar(JD_Arrivals_VGA_vec(1));
calendar0_VGA = sprintf('VGA, Days Past %d/%d/%d',mo,d,yr);

PlotBoi2(calendar0_launch,calendar0_VGA,16)
grid off
clabel(cs1 ,h1,'color',colors.std.ltred)
clabel(cs2 ,h2,'color',colors.std.blue)
clabel(cs3 ,h3,'color',colors.std.black)

legend([h1, h2, h3],'C_3 @ Launch, km^2/s^2','V_\infty @ VGA, km/s','TOF, days')

end % run_study_Launch_2_VGA == 1


%% =======================================================================
%%% Pork Chop Plot - VGA to EGA1 (leg2)
% ========================================================================
if run_study_VGA_2_EGA1 == 1
% ------------------------------------
%%% Setting up julian dates
% ------------------------------------
% JD_Departures_VGA_vec = JD_Arrivals_VGA_vec;
% JD_Arrivals_EGA1_vec  = linspace(JD_Departures_VGA_vec(end)+1, JD_Departures_VGA_vec(end)+400, nPoints);

JD_Departures_VGA_vec = JD_Arrivals_VGA_vec;
989
JD_Arrivals_EGA1_vec  = linspace(JD_Departures_VGA_vec(end)+1, JD_Departures_VGA_vec(end)+500, nPoints);


[JD_Departures_VGA, JD_Arrivals_EGA1] = meshgrid(JD_Departures_VGA_vec, JD_Arrivals_EGA1_vec);

% ------------------------------------
%%% Looping through options and computing vInf and c3
% ------------------------------------
X0s_sc_leg2{size(JD_Departures_VGA,1),size(JD_Departures_VGA,2)} = [];
vInfs_EGA1_m{size(JD_Departures_VGA,1),size(JD_Departures_VGA,2)} = [];
VInfs_EGA1_m = zeros(size(JD_Departures_VGA));
VInfs_VGA_p = zeros(size(JD_Departures_VGA));
flybyRadii_VGA = zeros(size(JD_Departures_VGA));
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
        
        %%% flyby radius at Venus
        vInf_Venus_m = vInfs_VGA_m{depIndex, arrIndex};
        VInf_Venus_m = (norm(vInf_Venus_m) + norm(vInf_Venus_p))/2;
        turningAngle_rad = acos(dot(vInf_Venus_m,vInf_Venus_p)./(norm(vInf_Venus_m)*norm(vInf_Venus_p)));
        [ rp ] = calcFlybyRadius( turningAngle_rad, VInf_Venus_m, Venus.u );
        flybyRadii_VGA(depIndex, arrIndex) = rp;

        %%% vInf at Earth arrival
        vInf_EGA1_m = vSc_EGA1_m - vEarth_A;
        vInfs_EGA1_m{depIndex, arrIndex} = vInf_EGA1_m;
        VInfs_EGA1_m(depIndex, arrIndex) = norm(vInf_EGA1_m);
        
        %%% TOFs
        TOFs_days_leg2(depIndex, arrIndex) = JD_Arr - JD_Dep;
        
        %%% Starting state for this leg of the trajectory [r0; v0]
        X0s_sc_leg2{depIndex, arrIndex} = [rVenus_D; vSc_VGA_p];

    end
end

% ------------------------------------
%%% Pork Chop Plot
% ------------------------------------
figure('position',[440 234 762 564]); hold all
vInf_VGA_contours = [4, 5,6,7,7.517 8, 12, 15, 17, 20];
[cs1, h1] = contour(JD_Departures_VGA-JD_Departures_VGA(1),JD_Arrivals_EGA1-JD_Arrivals_EGA1(1),VInfs_VGA_p,vInf_VGA_contours,'ShowText','on','color',colors.std.blue,'linewidth',1.2);

vInf_EGA1_contours = [7.5, 8, 8.5,9,10, 12, 15, 17, 20];
[cs2, h2] = contour(JD_Departures_VGA-JD_Departures_VGA(1),JD_Arrivals_EGA1-JD_Arrivals_EGA1(1),VInfs_EGA1_m,vInf_EGA1_contours,'ShowText','on','color',colors.std.ltred,'linewidth',1.2);

TOF_contours_leg2 = [260, 280, 300, 320, 340];
% [cs3,h3] = contour(JD_Departures_VGA-JD_Departures_VGA(1),JD_Arrivals_EGA1-JD_Arrivals_EGA1(1),TOFs_days,TOF_contours,'ShowText','on','color',colors.std.black,'linewidth',1.2);
[cs3,h3] = contour(JD_Departures_VGA-JD_Departures_VGA(1),JD_Arrivals_EGA1-JD_Arrivals_EGA1(1),TOFs_days_leg2,'ShowText','on','color',colors.std.black,'linewidth',1.2);

[yr, mo, d, ~, ~, ~] = julian2calendar(JD_Departures_VGA_vec(1));
calendar0_VGA = sprintf('VGA, Days Past %d/%d/%d',mo,d,yr);

[yr, mo, d, ~, ~, ~] = julian2calendar(JD_Arrivals_EGA1_vec(1));
calendar0_EGA1 = sprintf('EGA1, Days Past %d/%d/%d',mo,d,yr);

PlotBoi2(calendar0_VGA,calendar0_EGA1,16)
grid off
clabel(cs1 ,h1,'color',colors.std.blue)
clabel(cs2 ,h2,'color',colors.std.ltred)
clabel(cs3 ,h3,'color',colors.std.black)

legend([h1, h2, h3,],'V_\infty @ VGA, km/s','V_\infty @ EGA1, km/s','TOF, days')

end % run_study_VGA_2_EGA1 == 1


%% =======================================================================
%%% Look at candidates for Launch-VGA-EGA1 (leg 1 and leg 2)
% ========================================================================
if run_candidateSearch_prelim == 1
%%% Criteria
% -Arrival date matches departure date at VGA
% -Difference between VInf_m and VInf_p is less than some tolerance
% -Flyby radius is above some alitude tolerance
% -VInf and c3 requirements

% ------------------------------------
%%% Set requirements
% ------------------------------------
%%% C3 requirement on Earth Launch
req_c3_Launch_ub = 25; % km^2/s^2

%%% Tolerance on abs(|vInf_m| - |vInf_p|) at VGA
req_tol_VInf_VGA = 0.2; % km/s

%%% Flyby radius for Venus
req_altitude_VGA = 200; % km
req_flybyRadius_VGA = Venus.R + req_altitude_VGA; % km

% ------------------------------------
%%% Searching for good dates
% ------------------------------------
%%% Preallocating space for candidates
Index_solution = 0;
candidateStruct_prelim{10000000} = [];

% testCount1 = 0;
% testCount2 = 0;
% testCount3 = 0;
% testCount4 = 0;
%%% Check leg-1 of trajectory
for Index_Launch = 1:nPoints
    for Index_VGA_arr = 1:nPoints
        %%% Grab current launch window
        JDi_Launch = JD_Departures_Launch(Index_Launch, Index_VGA_arr);
        JDi_VGAArrival = JD_Arrivals_VGA(Index_Launch, Index_VGA_arr);
        
        %%% Skip this date if earth-departure C3 is too high
        if c3s_Launch_p(Index_Launch, Index_VGA_arr) > req_c3_Launch_ub
            continue
        end
        
        %%% Check leg-2 of trajectory
        for Index_VGA_dep = 1:nPoints
            for Index_EGA1_arr = 1:nPoints
                %%% Grab current launch window
                JDi_VGA = JD_Departures_VGA(Index_VGA_dep, Index_EGA1_arr);
                JDi_EGA1 = JD_Arrivals_EGA1(Index_VGA_dep, Index_EGA1_arr);
                
%                 if (JDi_VGA - JD_Arrivals_VGA_vec(1)) > 80
%                     testCount1 = testCount1 + 1;
%                 end
                
                %%% Make sure arrival and departure at VGA occur on the
                %%% same day
                if (JDi_VGA ~= JDi_VGAArrival)
                    continue
                end
                
%                 if (JDi_VGA - JD_Arrivals_VGA_vec(1)) > 80
%                     testCount2 = testCount2 + 1;
%                 end
                
                %%% Check magnitude of difference in vIns
                if abs(VInfs_VGA_p(Index_VGA_dep,Index_EGA1_arr) - VInfs_VGA_m(Index_Launch,Index_VGA_arr)) > req_tol_VInf_VGA
                    continue
                end
%                 
%                 if (JDi_VGA - JD_Arrivals_VGA_vec(1)) > 80
%                     testCount3 = testCount3 + 1;
%                 end

                %%% Skip this date if flyby radius too close to Venus
                if flybyRadii_VGA(Index_VGA_dep, Index_EGA1_arr) < req_flybyRadius_VGA
                    continue
                end
                
%                 if (JDi_VGA - JD_Arrivals_VGA_vec(1)) > 80
%                     testCount4 = testCount4 + 1;
%                 end

                % ------------------------------------
                %%% Store solution
                % ------------------------------------
                Index_solution = Index_solution + 1;
                
                %%% Relevant dates
                candidateStruct_prelim{Index_solution}.JD_Launch       = JDi_Launch;
                candidateStruct_prelim{Index_solution}.JD_VGA          = JDi_VGA;
                candidateStruct_prelim{Index_solution}.JD_EGA1         = JDi_EGA1;
                candidateStruct_prelim{Index_solution}.TOF_leg1_sec    = JDi_VGAArrival-JDi_Launch;
                candidateStruct_prelim{Index_solution}.TOF_leg2_sec    = JDi_EGA1-JDi_VGA;
                
                %%% Relevant values
                candidateStruct_prelim{Index_solution}.c3_earth_p     = c3s_Launch_p(Index_Launch, Index_VGA_arr);
                candidateStruct_prelim{Index_solution}.VInf_VGA_m     = VInfs_VGA_m(Index_Launch,Index_VGA_arr);
                candidateStruct_prelim{Index_solution}.VInf_VGA_p     = VInfs_VGA_p(Index_VGA_dep,Index_EGA1_arr);
                candidateStruct_prelim{Index_solution}.VInf_VGA_diff  = abs(VInfs_VGA_p(Index_VGA_dep,Index_EGA1_arr) - VInfs_VGA_m(Index_Launch,Index_VGA_arr));
                candidateStruct_prelim{Index_solution}.flybyRadii_VGA = flybyRadii_VGA(Index_VGA_dep, Index_EGA1_arr);
                candidateStruct_prelim{Index_solution}.VInf_EGA1_m    = VInfs_EGA1_m(Index_VGA_dep,Index_EGA1_arr);
                candidateStruct_prelim{Index_solution}.vInf_EGA1_m    = vInfs_EGA1_m{Index_VGA_dep,Index_EGA1_arr};
                candidateStruct_prelim{Index_solution}.leg1_X0        = X0s_sc_leg1{Index_Launch,Index_VGA_arr};
                candidateStruct_prelim{Index_solution}.leg2_X0        = X0s_sc_leg2{Index_VGA_dep,Index_EGA1_arr};

                
            end
        end
        
        
    end
end

% ------------------------------------
%%% Reassigning Data
% ------------------------------------
%%% Deleting empty cells
if length(candidateStruct_prelim(~cellfun('isempty',candidateStruct_prelim))) == length(candidateStruct_prelim)
    warning('Need to preallocate more space to candidateStruct_prelim')
end
candidateStruct_prelim = candidateStruct_prelim(~cellfun('isempty',candidateStruct_prelim));

%%% Reassigning results to get them out of cells
clear candidates_prelim
candidates_prelim = [candidateStruct_prelim{:}];
candidates_prelim = orderfields(candidates_prelim);

candidates_prelim_JD_EGA1           = [candidates_prelim.JD_EGA1];

candidates_prelim_JD_Launch         = [candidates_prelim.JD_Launch];
candidates_prelim_JD_VGA            = [candidates_prelim.JD_VGA];
% candidates_prelim_TOF_leg1_sec      = [candidates_prelim.TOF_leg1_sec];
% candidates_prelim_TOF_leg2_sec      = [candidates_prelim.TOF_leg2_sec];
% candidates_prelim_c3_earth_p        = [candidates_prelim.c3_earth_p];
% candidates_prelim_VInf_VGA_m        = [candidates_prelim.VInf_VGA_m];
% candidates_prelim_VInf_VGA_p        = [candidates_prelim.VInf_VGA_p];
% candidates_prelim_VInf_VGA_diff     = [candidates_prelim.VInf_VGA_diff];
% candidates_prelim_flybyRadii        = [candidates_prelim.flybyRadii_VGA];
% candidates_prelim_VInf_EGA1_m       = [candidates_prelim.VInf_EGA1_m];

end % candidates_prelim





%% =======================================================================
%%% Take candidates so far and use them to form a PCP for leg4 to get new
%%% class of candidates
% ========================================================================
if run_study_EGA2_2_SOI == 1
% ------------------------------------
%%% Setting up julian dates
% ------------------------------------
%%% Get range of predicted EGA2 JDs, create range of JD_SOIs, and make a
%%% mesh
JD_Departures_EGA2_vec = sort(unique(candidates_prelim_JD_EGA1)) + tof_resonance_day;
nPoints2 = length(JD_Departures_EGA2_vec);
% JD_Arrivals_SOI_vec = linspace(JD_Departures_EGA2_vec(1) + 4.5*year_Earth_day, JD_2041, nPoints2);
JD_Arrivals_SOI_vec = linspace(JD_Departures_EGA2_vec(1) + 2.5*year_Earth_day, JD_2041, nPoints2);
warning('Maybe populate JD_Arrivals_SOI_vec a lot more')
[JD_Departures_EGA2, JD_Arrivals_SOI] = meshgrid(JD_Departures_EGA2_vec, JD_Arrivals_SOI_vec);

% ------------------------------------
%%% Looping through options and computing vInf and c3
% ------------------------------------
X0s_sc_leg4{size(JD_Departures_EGA2,1),size(JD_Departures_EGA2,2)} = [];
vInfs_EGA2_p{size(JD_Departures_EGA2,1),size(JD_Departures_EGA2,2)} = [];
VInfs_SOI_m = zeros(size(JD_Departures_EGA2));
VInfs_EGA2_p = zeros(size(JD_Departures_EGA2));
TOFs_days_leg4 = zeros(size(JD_Departures_EGA2));

for depIndex = 1:nPoints2
    for arrIndex = 1:nPoints2
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
        vInfs_EGA2_p{depIndex, arrIndex} = vInf_Earth_p;
        VInfs_EGA2_p(depIndex, arrIndex) = norm(vInf_Earth_p);

        %%% vInf at Saturn arrival
        vInf_Saturn = vSc_SOI_m - vSaturn_A;
        VInfs_SOI_m(depIndex, arrIndex) = norm(vInf_Saturn);
        
        %%% TOFs
        TOFs_days_leg4(depIndex, arrIndex) = JD_Arr - JD_Dep;
        
        %%% Starting state for this leg of the trajectory [r0; v0]
        X0s_sc_leg4{depIndex, arrIndex} = [rEarth_D; vSc_EGA2_p];

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
% [cs3,h3] = contour(JD_Departures_EGA2-JD_Departures_EGA2(1),JD_Arrivals_SOI-JD_Arrivals_SOI(1),TOFs_days,TOF_contours_leg4,'ShowText','on','color',colors.std.black,'linewidth',1.2);
[cs3,h3] = contour(JD_Departures_EGA2-JD_Departures_EGA2(1),JD_Arrivals_SOI-JD_Arrivals_SOI(1),TOFs_days_leg4,'ShowText','on','color',colors.std.black,'linewidth',1.2);

[yr, mo, d, ~, ~, ~] = julian2calendar(JD_Departures_EGA2_vec(1));
calendar0_EGA2 = sprintf('EGA2, Days Past %d/%d/%d',mo,d,yr);

[yr, mo, d, ~, ~, ~] = julian2calendar(JD_Arrivals_SOI_vec(1));
calendar0_SOI = sprintf('SOI, Days Past %d/%d/%d',mo,d,yr);

PlotBoi2(calendar0_EGA2,calendar0_SOI,16)
grid off
clabel(cs1 ,h1,'color',colors.std.blue)
clabel(cs2 ,h2,'color',colors.std.ltred)
clabel(cs3 ,h3,'color',colors.std.black)

legend([h1, h2, h3,],'V_\infty @ EGA2, km/s','V_\infty @ SOI, km/s','TOF, days')

end % run_study_EGA2_2_SOI == 1


%% =======================================================================
%%% Look at candidates for leg 4 and combine with preliminary candidates
% ========================================================================
if run_candidateSearch_mid == 1
%%% Criteria
% -vInf at SOI < 9 km/s
% -TOF < 15 years

% ------------------------------------
%%% Set requirements
% ------------------------------------
%%% Upper bound of VInf_SOI
req_vInf_SOI = 9; % km/s

%%% Upper bound of total travel time
req_TOF_yrs = 15;

% ------------------------------------
%%% Searching for candidates
% ------------------------------------
%%% Preallocating space for candidates
Index_solution = 0;
candidateStruct_mid{10000000} = [];

%%% Loop through preliminary candidates
for Index_cand_pre = 1:length(candidates_prelim)
    
    %%% Loop through leg-4 porkchop plot results
    for Index_EGA2_dep = 1:nPoints2
        for Index_SOI_arr = 1:nPoints2
            %%% Check to make sure it's the same EGA schedule as the
            %%% preliminary candidate

            if JD_Departures_EGA2(Index_EGA2_dep, Index_SOI_arr) ~= (candidates_prelim(Index_cand_pre).JD_EGA1 + tof_resonance_day)
                continue
            end
            
            %%% Check VInf_SOI requirement
            if VInfs_SOI_m(Index_EGA2_dep, Index_SOI_arr) > req_vInf_SOI
                continue
            end
            
            %%% Check TOF requirement
            if (candidates_prelim(Index_cand_pre).JD_EGA1 + tof_resonance_day + TOFs_days_leg4(Index_EGA2_dep, Index_SOI_arr) - candidates_prelim(Index_cand_pre).JD_Launch)/year_Earth_day > req_TOF_yrs
                continue
            end
            
            % ------------------------------------
            %%% Store solution
            % ------------------------------------
            Index_solution = Index_solution + 1;
            
            %%% Store all data from preliminary candidate
            candidateStruct_mid{Index_solution} = candidates_prelim(Index_cand_pre);
            
            %%% Add new data relating to leg-4
            candidateStruct_mid{Index_solution}.JD_EGA2          = JD_Departures_EGA2(Index_EGA2_dep, Index_SOI_arr);
            candidateStruct_mid{Index_solution}.JD_SOI           = JD_Departures_EGA2(Index_EGA2_dep, Index_SOI_arr) + TOFs_days_leg4(Index_EGA2_dep, Index_SOI_arr);
            candidateStruct_mid{Index_solution}.vInf_EGA2_p      = vInfs_EGA2_p{Index_EGA2_dep, Index_SOI_arr};
            candidateStruct_mid{Index_solution}.VInf_EGA2_p      = VInfs_EGA2_p(Index_EGA2_dep, Index_SOI_arr);
            candidateStruct_mid{Index_solution}.VInf_SOI_m       = VInfs_SOI_m(Index_EGA2_dep, Index_SOI_arr);
            candidateStruct_mid{Index_solution}.TOF_leg4_sec     = TOFs_days_leg4(Index_EGA2_dep, Index_SOI_arr)*86400;
            candidateStruct_mid{Index_solution}.leg4_X0          = X0s_sc_leg4{Index_EGA2_dep, Index_SOI_arr};
            
        end
    end

end

% ------------------------------------
%%% Reassigning Data
% ------------------------------------
%%% Deleting empty cells
if length(candidateStruct_mid(~cellfun('isempty',candidateStruct_mid))) == length(candidateStruct_mid)
    warning('Need to preallocate more space to candidateStruct_mid')
end
candidateStruct_mid = candidateStruct_mid(~cellfun('isempty',candidateStruct_mid));

%%% Reassigning results to get them out of cells
clear candidates_mid
candidates_mid = [candidateStruct_mid{:}];
candidates_mid = orderfields(candidates_mid);


%%% Seperating some data for convenience of analysys
candidates_mid_JD_EGA1           = [candidates_mid.JD_EGA1];
candidates_mid_JD_Launch         = [candidates_mid.JD_Launch];
candidates_mid_JD_VGA            = [candidates_mid.JD_VGA];
end

%% =======================================================================
%%% Look at candidates that have small discrepancy at VInf_EGA2 and save as
%%% final candidates... Have to do full resonance analysis at each
%%% candidate
% ========================================================================
if run_resonanceStudy_and_candidateSearch_final == 1
%%% Criteria
% -abs(VInf_EGA2_m - VInf_EGA2_p) must be < tolerance
% -EGA flyby alitude

% ------------------------------------
%%% Set requirements
% ------------------------------------
%%% Upper bound of VInf_SOI
req_tol_VInf_EGA2 = 0.2; % km/s

%%% Earth flyby altitude requirement
req_Earth_flyby_altitude = 150; % km

% ------------------------------------
%%% Setup for resonance study
% ------------------------------------
%%% phi values
n_phis = 360;
phis_rad = linspace(0,2*pi,n_phis);

% ------------------------------------
%%% Searching for candidates
% ------------------------------------
%%% Preallocating space for candidates
Index_solution = 0;
candidateStruct_final{10000000} = [];

warning('Could be a mistake the way I grab phi results here')
%%% Loop through mid-stage candidates
tic

problem1_counter  = 0;
problem2_counter  = 0;
problem2_minValue = 1000000000;
length(candidates_mid)
for Index_cand_mid = 1:length(candidates_mid)
    %%% Grabbing data from current candidate for clarity
    vInf_EGA1_m_jj  = candidates_mid(Index_cand_mid).vInf_EGA1_m;
    VInf_EGA1_m_jj  = candidates_mid(Index_cand_mid).VInf_EGA1_m;
    JD_EGA1_jj      = candidates_mid(Index_cand_mid).JD_EGA1;
    JD_EGA2_jj      = candidates_mid(Index_cand_mid).JD_EGA2;
    JD_SOI_jj       = candidates_mid(Index_cand_mid).JD_SOI;
    TOF_leg4_sec_jj = candidates_mid(Index_cand_mid).TOF_leg4_sec;
    vInf_EGA2_p_jj  = candidates_mid(Index_cand_mid).vInf_EGA2_p;
    VInf_EGA2_p_jj  = candidates_mid(Index_cand_mid).VInf_EGA2_p;
    
    %%% Earth's state at current EGA1
    [rEarth_EGA1_jj, vEarth_EGA1_jj]  = JulianDate_2_HeliocentricState(JD_EGA1_jj, 'Earth'); % km, km/s
    
    %%% Earth's state at current EGA2
    [rEarth_EGA2_jj, vEarth_EGA2_jj]  = JulianDate_2_HeliocentricState(JD_EGA2_jj, 'Earth'); % km, km/s
    
    %%% Saturn's state at current SOI
    [rSaturn_SOI_jj, vSaturn_SOI_jj]  = JulianDate_2_HeliocentricState(JD_SOI_jj, 'Saturn'); % km, km/s
    
    %%% Lambert solution from current EGA2 to SOI
    [vSc_EGA2_p_jj, vSc_SOI_m_jj, ~]  = lambertSolver(rEarth_EGA2_jj, rSaturn_SOI_jj, TOF_leg4_sec_jj, 0, 0, Sun.u); % km/s
    
    %%% Lambert solve EGA1 to EGA2
    [vSc_EGA1_p_jj, vSc_EGA2_m_jj, ~] = lambertSolver(rEarth_EGA1_jj, rEarth_EGA2_jj, tof_resonance_sec, 0, 0, Sun.u); % km/s

    %%% EGA1
    VSc_EGA1_p_jj     = norm(vSc_EGA1_p_jj);
    VEarth_EGA1_jj    = norm(vEarth_EGA1_jj);
    theta_EGA1_rad_jj = acos((-VSc_EGA1_p_jj^2 + VInf_EGA1_m_jj^2 + VEarth_EGA1_jj^2)/(2*VInf_EGA1_m_jj*VEarth_EGA1_jj));
    
    %%% EGA2
    VSc_EGA2_p_jj  = norm(vSc_EGA2_p_jj);
    VEarth_EGA2_jj = norm(vEarth_EGA2_jj);
    
    %%% Preallocating for this batch of phis
    rps_EGA1_kk     = zeros(n_phis,1);
    rps_EGA2_kk     = zeros(n_phis,1);
    vInfs_EGA1_p_kk = zeros(n_phis,3);
    vInfs_EGA2_m_kk = zeros(n_phis,3);

    for kk = 1:n_phis
        %%% Set phi
        phi_rad = phis_rad(kk);

        % ---------------------
        %%% EGA 1
        % ---------------------
        %%% Compute v-inf plus
        vInf_EGA1_p_VNC_kk = VInf_EGA1_m_jj * [cos(pi-theta_EGA1_rad_jj); sin(pi-theta_EGA1_rad_jj)*cos(phi_rad); -sin(pi-theta_EGA1_rad_jj)*sin(phi_rad)];

        %%% Rotate from VNC to inertial
        vHat_EGA1 = vEarth_EGA1_jj / VEarth_EGA1_jj;
        hHat_EGA1 = cross(rEarth_EGA1_jj,vEarth_EGA1_jj)/norm(cross(rEarth_EGA1_jj,vEarth_EGA1_jj));
        nHat_EGA1 = hHat_EGA1;
        cHat_EGA1 = cross(vHat_EGA1,nHat_EGA1);
        T_VNC2Ecliptic_EGA1 = [vHat_EGA1, nHat_EGA1, cHat_EGA1];
        vInf_EGA1_p_kk = T_VNC2Ecliptic_EGA1 * vInf_EGA1_p_VNC_kk;
        VInf_EGA1_p_kk = norm(vInf_EGA1_p_kk);
        vInfs_EGA1_p_kk(kk,:) = vInf_EGA1_p_kk';

        %%% Calculate turning angle
        turningAngle_EGA1_rad_kk = acos(dot(vInf_EGA1_m_jj,vInf_EGA1_p_kk)/(VInf_EGA1_m_jj*VInf_EGA1_p_kk));

        %%% Calculate EGA1 flyby radius
        rps_EGA1_kk(kk) = (Earth.u/(VInf_EGA1_m_jj*VInf_EGA1_p_kk))*((1/cos((pi-turningAngle_EGA1_rad_kk)/2)) - 1);

        % ---------------------
        %%% EGA 2
        % ---------------------
        %%% Theta
        vInf_EGA2_m_kk = vInf_EGA1_p_kk + vEarth_EGA1_jj - vEarth_EGA2_jj;
        VInf_EGA2_m_kk = norm(vInf_EGA2_m_kk);
        vInfs_EGA2_m_kk(kk,:) = vInf_EGA2_m_kk';

        %%% Calculate turning angle
        turningAngle_EGA2_kk = acos(dot(vInf_EGA2_m_kk,vInf_EGA2_p_jj)/(VInf_EGA2_m_kk*VInf_EGA2_p_jj));

        %%% Calculate EGA2 flyby radius
        rps_EGA2_kk(kk) = (Earth.u/(VInf_EGA2_m_kk*VInf_EGA2_p_jj))*((1/cos((pi-turningAngle_EGA2_kk)/2)) - 1);

    end
    


    % ---------------------
    %%% Checking if altitude requirement met
    % ---------------------
    %%% Max rps of both flybys
    max_rp_EGA1_index_jj = find(rps_EGA1_kk == max(rps_EGA1_kk));
    max_rp_EGA2_index_jj = find(rps_EGA2_kk == max(rps_EGA2_kk));
    
    [ smallest_rpMax_jj ] = smallerValue( rps_EGA1_kk(max_rp_EGA1_index_jj), rps_EGA2_kk(max_rp_EGA2_index_jj) ); % km
    
    if smallest_rpMax_jj < (req_Earth_flyby_altitude + Earth.R)
        problem1_counter = problem1_counter + 1;
        continue
    end
    
    
%         %%% Plotting
%         figure; hold all
%     %     p4 = fill([72.08, 72.08, 114.52, 114.52], [0, 12000,12000,0], 'g');
%         p1 = plot(phis_rad.*(180/pi),rps_EGA1_kk,'linewidth',2,'color',colors.std.blue);
%         p2 = plot(phis_rad.*(180/pi),rps_EGA2_kk,'linewidth',2,'color',colors.std.ltred);
%         % p3 = plot([0,360],[Earth.R, Earth.R],'--k','linewidth',1);
%         p3 = plot([0,360],[Earth.R+100, Earth.R+100],'--k','linewidth',1);
%     %     p5 = plot([myPhi_deg,myPhi_deg],[0, 12000],'m','linewidth',1.5);
%         PlotBoi2('$\Phi$, $^\circ$','Flyby Radius, $km$',18,'LaTex')
%         xlim([0 360])
%         legend([p1 p2 p3],'EGA1','EGA2','100 km Altitude')
%     %     legend([p1 p2 p3 p4 p5],'EGA1','EGA2','300 km Altitude','Range of Acceptable \Phi', 'Selected \Phi')
    
    % ---------------------
    %%% Grabbing info from chosen phi value (lesser of two max flyby radii)
    % ---------------------
    
    
    %%% Get index of chosen values
    if max_rp_EGA1_index_jj > max_rp_EGA2_index_jj
        myPhi_index_jj = max_rp_EGA2_index_jj;
    elseif max_rp_EGA1_index_jj < max_rp_EGA2_index_jj
        myPhi_index_jj = max_rp_EGA1_index_jj;
    else
       myPhi_index_jj = max_rp_EGA1_index_jj;
    end
    
    
    % ---------------------
    %%% Checking if |VInf_EGA2| requirement meant
    % ---------------------
    VInf_EGA2_m_jj = norm(vInfs_EGA2_m_kk(myPhi_index_jj,:));
    VInf_EGA2_p_jj = candidates_mid(Index_cand_mid).VInf_EGA2_p;
    
    if abs(VInf_EGA2_p_jj - VInf_EGA2_m_jj) > req_tol_VInf_EGA2
        problem2_counter = problem2_counter + 1;
        if abs(VInf_EGA2_p_jj - VInf_EGA2_m_jj) < problem2_minValue
            problem2_minValue = abs(VInf_EGA2_p_jj - VInf_EGA2_m_jj);
        end
        continue
    end

    % ---------------------
    %%% Storing candidate info
    % ---------------------
    Index_solution = Index_solution + 1;
    
    candidateStruct_final{Index_solution} = candidates_mid(Index_cand_mid);
    
    candidateStruct_final{Index_solution}.vInf_EGA1_p = vInfs_EGA1_p_kk(myPhi_index_jj,:);
    candidateStruct_final{Index_solution}.VInf_EGA1_p = norm(vInfs_EGA1_p_kk(myPhi_index_jj,:));
    candidateStruct_final{Index_solution}.vInf_EGA2_m = vInfs_EGA2_m_kk(myPhi_index_jj,:);
    candidateStruct_final{Index_solution}.VInf_EGA2_m = VInf_EGA2_m_jj;
    candidateStruct_final{Index_solution}.phi_rad     = phis_rad(myPhi_index_jj);
    
    candidateStruct_final{Index_solution}.TOF_leg3_sec = tof_resonance_sec;
    candidateStruct_final{Index_solution}.leg3_X0      = [rEarth_EGA1_jj; vInfs_EGA1_p_kk(myPhi_index_jj,:)' + vEarth_EGA1_jj];
    
    candidateStruct_final{Index_solution}.rp_EGA1     = rps_EGA1_kk(myPhi_index_jj);
    candidateStruct_final{Index_solution}.rp_EGA2     = rps_EGA2_kk(myPhi_index_jj);
    
    candidateStruct_final{Index_solution}.VInf_EGA2_diff = abs(VInf_EGA2_p_jj - VInf_EGA2_m_jj);
end
toc


% ------------------------------------
%%% Reassigning Data
% ------------------------------------
%%% Deleting empty cells
if length(candidateStruct_final(~cellfun('isempty',candidateStruct_final))) == length(candidateStruct_final)
    warning('Need to preallocate more space to candidateStruct_final')
end
candidateStruct_final = candidateStruct_final(~cellfun('isempty',candidateStruct_final));

%%% Reassigning results to get them out of cells
clear candidates_final
candidates_final = [candidateStruct_final{:}];
candidates_final = orderfields(candidates_final);


candidates_final_VInf_EGA2_diff = [candidates_final.VInf_EGA2_diff];

candidates_final_JD_Launch = [candidates_final.JD_Launch];
candidates_final_JD_VGA = [candidates_final.JD_VGA];




idx = 60;
% idx = 18037;

%%% Choosing ode tolerance
tol = 1e-13;
%%% Setting integrator options
options = odeset('RelTol',tol,'AbsTol',tol);
n_dt = 1000;

time0_leg1 = linspace(candidates_mid(idx).JD_Launch,candidates_mid(idx).JD_VGA,n_dt).*86400;
time0_leg2 = linspace(candidates_mid(idx).JD_VGA,candidates_mid(idx).JD_EGA1,n_dt).*86400;
time0_leg4 = linspace(candidates_mid(idx).JD_EGA2,candidates_mid(idx).JD_SOI,n_dt).*86400;

X0_leg1 = candidates_mid(idx).leg1_X0;
X0_leg2 = candidates_mid(idx).leg2_X0;
X0_leg4 = candidates_mid(idx).leg4_X0;

% time0_leg1 = linspace(candidates_final(idx).JD_Launch,candidates_final(idx).JD_VGA,n_dt).*86400;
% time0_leg2 = linspace(candidates_final(idx).JD_VGA,candidates_final(idx).JD_EGA1,n_dt).*86400;
% time0_leg3 = linspace(candidates_final(idx).JD_EGA1,candidates_final(idx).JD_EGA2,n_dt).*86400;
% time0_leg4 = linspace(candidates_final(idx).JD_EGA2,candidates_final(idx).JD_SOI,n_dt).*86400;
% 
% X0_leg1 = candidates_final(idx).leg1_X0;
% X0_leg2 = candidates_final(idx).leg2_X0;
% X0_leg3 = candidates_final(idx).leg3_X0;
% X0_leg4 = candidates_final(idx).leg4_X0;

[~, X_sc_leg1] = ode113(@Int_2BI, time0_leg1, X0_leg1, options, Sun.u);
[~, X_sc_leg2] = ode113(@Int_2BI, time0_leg2, X0_leg2, options, Sun.u);
% [~, X_sc_leg3] = ode113(@Int_2BI, time0_leg3, X0_leg3, options, Sun.u);
[~, X_sc_leg4] = ode113(@Int_2BI, time0_leg4, X0_leg4, options, Sun.u);


figure; hold all
plot3(X_sc_leg1(:,1),X_sc_leg1(:,2),X_sc_leg1(:,3),'color',colors.std.ltred)
plot3(X_sc_leg2(:,1),X_sc_leg2(:,2),X_sc_leg2(:,3),'color',colors.std.blue)
% plot3(X_sc_leg3(:,1),X_sc_leg3(:,2),X_sc_leg3(:,3),'color',colors.std.mag)
plot3(X_sc_leg4(:,1),X_sc_leg4(:,2),X_sc_leg4(:,3),'color',colors.std.grn)

axis equal
PlotBoi3('X, km','Y, km', 'Z, km', 18, 'LaTex')


warning('Check scheme for Phi selection')




end % run_candidateSearch_final == 1

















