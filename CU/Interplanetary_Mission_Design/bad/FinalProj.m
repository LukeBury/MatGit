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
run_study_Launch_2_VGA        = 1;
run_study_VGA_2_EGA1          = 1;
run_candidateSearch_leg1_leg2 = 1;
run_study_guess_EGA2_2_SOI    = 1;
run_study_EGA1_2_EGA2         = 1;
run_study_EGA2_2_SOI          = 0;

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
%%% Julian Dates
% ------------------------------------
% [ JD ] = calendar2julian(yr,  mo, d, h, min, s)
% [Year, Mon, Day, hr, min, s] = julian2calendar(JD)
% 
% [ JD ] = calendar2julian(2025, 1, 1, 0, 0, 0)


JD_2025 = 2460676.5;

JD_guess_VGA = JD_2025 + 125;

JD_guess_EGA1 = JD_guess_VGA + 300;

JD_guess_EGA2 = JD_guess_EGA1 + 365.242189;

JD_guess_SOI = JD_guess_EGA2 + 2222.85;


JD_Launch_guess = 2460676.5;
JD_VGA_guess    = 2460774.5;
% JD_EGA1_guess   = 2460986.5721;
JD_EGA1_guess   = 2461080.2021;
% JD_EGA2_guess   = JD_EGA1_guess + Earth.Tp*2/86400;

%% =======================================================================
%%% Pork Chop Plot - Launch to VGA (leg1)
% ========================================================================
if run_study_Launch_2_VGA == 1
% ------------------------------------
%%% Setting up julian dates
% ------------------------------------
nPoints = 250;
JD_Departures_Launch_vec = linspace(JD_2025, JD_2025+100, nPoints);
JD_Arrivals_VGA_vec     = linspace(JD_guess_VGA-50, JD_guess_VGA+130, nPoints);
% JD_Departures_Launch_vec = linspace(JD_2025, JD_2025+365, nPoints);
% JD_Arrivals_VGA_vec     = linspace(JD_guess_VGA-50, JD_guess_VGA+550, nPoints);

[JD_Departures_Launch, JD_Arrivals_VGA] = meshgrid(JD_Departures_Launch_vec, JD_Arrivals_VGA_vec);

% ------------------------------------
%%% Looping through options and computing vInf and c3
% ------------------------------------
X0s_sc_leg1{size(JD_Departures_Launch,1),size(JD_Departures_Launch,2)} = [];
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
        X0s_sc_leg1{xIndex, yIndex} = [rEarth_D; vSc_EarthDeparture];

    end
end



% ------------------------------------
%%% Pork Chop Plot
% ------------------------------------
figure('position',[74 233 762 564]); hold all
vInf_contours = [3.8, 4, 5,7.517, 8, 12, 15, 17, 20];
[cs1, h1] = contour(JD_Departures_Launch-JD_Departures_Launch(1),JD_Arrivals_VGA-JD_Arrivals_VGA(1),VInfs_VGA_m,vInf_contours,'ShowText','on','color',colors.std.blue,'linewidth',1.2);
% [cs1, h1] = contour(JD_Departures_Launch-JD_Departures_Launch(1),JD_Arrivals_VGA-JD_Arrivals_VGA(1),VInfs_VGA_m,'ShowText','on','color',colors.std.blue,'linewidth',1.2);

c3_contours = [5,7.2,8, 10, 15, 20, 25, 30, 35, 40, 45, 50];
[cs2, h2] = contour(JD_Departures_Launch-JD_Departures_Launch(1),JD_Arrivals_VGA-JD_Arrivals_VGA(1),c3s_Launch_p,c3_contours,'ShowText','on','color',colors.std.ltred,'linewidth',1.2);
% [cs2, h2] = contour(JD_Departures_Launch-JD_Departures_Launch(1),JD_Arrivals_VGA-JD_Arrivals_VGA(1),c3s_Launch_p,'ShowText','on','color',colors.std.ltred,'linewidth',1.2);

TOF_contours = [80, 100, 120, 140, 160, 200];
% [cs3,h3] = contour(JD_Departures_Launch-JD_Departures_Launch(1),JD_Arrivals_VGA-JD_Arrivals_VGA(1),TOFs_days,TOF_contours,'ShowText','on','color',colors.std.black,'linewidth',1.2);
[cs3,h3] = contour(JD_Departures_Launch-JD_Departures_Launch(1),JD_Arrivals_VGA-JD_Arrivals_VGA(1),TOFs_days,'ShowText','on','color',colors.std.black,'linewidth',1.2);

[yr, mo, d, ~, ~, ~] = julian2calendar(JD_Departures_Launch_vec(1));
calendar0_launch = sprintf('Launch, Days Past %d/%d/%d',mo,d,yr);

[yr, mo, d, ~, ~, ~] = julian2calendar(JD_Arrivals_VGA_vec(1));
calendar0_VGA = sprintf('VGA, Days Past %d/%d/%d',mo,d,yr);

%%% Plotting current Launch and VGA choice
plot(JD_Launch_guess - JD_Departures_Launch_vec(1),JD_VGA_guess - JD_Arrivals_VGA_vec(1),'m.','markersize',16) 

PlotBoi2(calendar0_launch,calendar0_VGA,16)
grid off
clabel(cs1 ,h1,'color',colors.std.blue)
clabel(cs2 ,h2,'color',colors.std.ltred)
clabel(cs3 ,h3,'color',colors.std.black)

legend([h1, h2, h3],'V_\infty @ VGA, km/s','C_3 @ Launch, km^2/s^2','TOF, days')

end % run_study_Launch_2_VGA == 1

%% =======================================================================
%%% Pork Chop Plot - VGA to EGA1 (leg2)
% ========================================================================
if run_study_VGA_2_EGA1 == 1
% ------------------------------------
%%% Setting up julian dates
% ------------------------------------
% JD_Departures_VGA_vec = linspace(JD_guess_VGA-30, JD_guess_VGA+30, nPoints);
% JD_Arrivals_EGA1_vec  = linspace(JD_guess_EGA1-150, JD_guess_EGA1+30, nPoints);

JD_Departures_VGA_vec = JD_Arrivals_VGA_vec;
JD_Arrivals_EGA1_vec  = linspace(JD_guess_EGA1-169, JD_guess_EGA1+100, nPoints);
[JD_Departures_VGA, JD_Arrivals_EGA1] = meshgrid(JD_Departures_VGA_vec, JD_Arrivals_EGA1_vec);

% ------------------------------------
%%% Looping through options and computing vInf and c3
% ------------------------------------
X0s_sc_leg2{size(JD_Departures_VGA,1),size(JD_Departures_VGA,2)} = [];
VInfs_EGA1_m = zeros(size(JD_Departures_VGA));
VInfs_VGA_p = zeros(size(JD_Departures_VGA));
flybyRadii_VGA = zeros(size(JD_Departures_VGA));
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
        
        %%% flyby radius at Venus
        vInf_Venus_m = vInfs_VGA_m{xIndex, yIndex};
        VInf_Venus_m = (norm(vInf_Venus_m) + norm(vInf_Venus_p))/2;
        turningAngle_rad = acos(dot(vInf_Venus_m,vInf_Venus_p)./(norm(vInf_Venus_m)*norm(vInf_Venus_p)));
        [ rp ] = calcFlybyRadius( turningAngle_rad, VInf_Venus_m, Venus.u );
        flybyRadii_VGA(xIndex, yIndex) = rp;

        %%% vInf at Earth arrival
        vInf_Earth = vSc_EarthArrival - vEarth_A;
        VInfs_EGA1_m(xIndex, yIndex) = norm(vInf_Earth);
        
        %%% TOFs
        TOFs_days(xIndex, yIndex) = JD_A - JD_D;
        
        %%% Starting state for this leg of the trajectory [r0; v0]
        X0s_sc_leg2{xIndex, yIndex} = [rVenus_D; vSc_VenusDeparture];

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

TOF_contours = [260, 280, 300, 320, 340];
% [cs3,h3] = contour(JD_Departures_VGA-JD_Departures_VGA(1),JD_Arrivals_EGA1-JD_Arrivals_EGA1(1),TOFs_days,TOF_contours,'ShowText','on','color',colors.std.black,'linewidth',1.2);
[cs3,h3] = contour(JD_Departures_VGA-JD_Departures_VGA(1),JD_Arrivals_EGA1-JD_Arrivals_EGA1(1),TOFs_days,'ShowText','on','color',colors.std.black,'linewidth',1.2);

[yr, mo, d, ~, ~, ~] = julian2calendar(JD_Departures_VGA_vec(1));
calendar0_VGA = sprintf('VGA, Days Past %d/%d/%d',mo,d,yr);

[yr, mo, d, ~, ~, ~] = julian2calendar(JD_Arrivals_EGA1_vec(1));
calendar0_EGA1 = sprintf('EGA1, Days Past %d/%d/%d',mo,d,yr);

%%% Plotting line resulting from current VGA choice
plot([1,1].*(JD_VGA_guess - JD_Departures_VGA_vec(1)),[JD_Arrivals_EGA1_vec(1), JD_Arrivals_EGA1_vec(end)]-JD_Arrivals_EGA1_vec(1),'m') 

%%% Plotting current EGA choice
plot(JD_VGA_guess - JD_Departures_VGA_vec(1),JD_EGA1_guess-JD_Arrivals_EGA1_vec(1),'m.','markersize',16) 


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
if run_candidateSearch_leg1_leg2 == 1
%%% Criteria
% -Arrival date matches departure date at flyby planet
% -Difference between vInf_m and vInf_p is less than some tolerance
% -Flyby radius is above some alitude tolerance
% -vInf and c3 requirements

% ------------------------------------
%%% Set requirements
% ------------------------------------
%%% Lower and upper bounds on earth launch date
[ req_JD_Launch_lb ] = JD_2025;
[ req_JD_Launch_ub ] = JD_2025 + 365.242189;

%%% C3 requirement on Earth Launch
req_c3_Launch_ub = 25; % km^2/s^2

%%% Lower and uppder bounds on Venus arrival date
% % % % % % % % % % % % % % % % req_JD_ven_lb = JD_Departures_VGA_vec(1);
% % % % % % % % % % % % % % % % req_JD_ven_ub = JD_Departures_VGA_vec(end);
% req_JD_ven_lb = JD_Arrivals_VGA_vec(1);
% req_JD_ven_ub = JD_Arrivals_VGA_vec(end);

%%% Tolerance on abs(|vInf_m| - |vInf_p|) at Jupiter
req_tol_vInf_ven = 0.05; % km/s

% %%% Lower and uppder bounds on Earth arrival date
% req_JD_EGA1_lb = JD_Arrivals_EGA1_vec(1);
% req_JD_EGA1_ub = JD_Arrivals_EGA1_vec(end);

%%% Flyby radius for Jupiter
req_flybyRadius = Venus.R + 250; % km
% % % % % % % JD_Departures_Launch, JD_Arrivals_VGA, JD_Departures_VGA, JD_Arrivals_EGA1
% ------------------------------------
%%% Searching for good dates
% ------------------------------------
solutionIndex_1 = 0;
candidateStruct_leg1_leg2{100000} = [];

for LaunchIndex = 1:nPoints
    for venusArrivalIndex = 1:nPoints
        %%% Grab current launch window
        JDi_Launch = JD_Departures_Launch(LaunchIndex, venusArrivalIndex);
        JDi_VGAArrival = JD_Arrivals_VGA(LaunchIndex, venusArrivalIndex);
% %         JD_Arrivals_VGA_vec     = linspace(JD_guess_VGA-50, JD_guess_VGA+150, nPoints);
        %%% Skip this date if it's not in the proper launch window
        if (JDi_Launch < req_JD_Launch_lb) || (JDi_Launch > req_JD_Launch_ub)
            continue
        end
        
        %%% Skip this date if earth-departure C3 is too high
        if c3s_Launch_p(LaunchIndex, venusArrivalIndex) > req_c3_Launch_ub
            continue
        end
        
%         %%% Skip this date if it's not in the proper launch window
%         if (JDi_venusArrival < req_JD_ven_lb) || (JDi_venusArrival > req_JD_ven_ub)
%             continue
%         end
                
        %%% check leg-2 of trajectory
        for VGADepartureIndex = 1:nPoints
            for EGA1ArrivalIndex = 1:nPoints
                %%% Grab current launch window
                JDi_VGADeparture = JD_Departures_VGA(VGADepartureIndex, EGA1ArrivalIndex);
                JDi_EGA1Arrival = JD_Arrivals_EGA1(VGADepartureIndex, EGA1ArrivalIndex);
                
                %%% Make sure arrival and departure at Jupiter occur on the
                %%% same day
                if (JDi_VGADeparture ~= JDi_VGAArrival)
                    continue
                end
                
                %%% Check magnitude of difference in vIns
                if abs(VInfs_VGA_p(VGADepartureIndex,EGA1ArrivalIndex) - VInfs_VGA_m(LaunchIndex,venusArrivalIndex)) > req_tol_vInf_ven
                    continue
                end
                
%                 %%% Skip this date if it's not in the proper  window
%                 if (JDi_EGA1Arrival < req_JD_EGA1_lb) || (JDi_EGA1Arrival > req_JD_EGA1_ub)
%                     continue
%                 end
                
                
                %%% Skip this date if flyby radius too close to Venus
                if flybyRadii_VGA(VGADepartureIndex, EGA1ArrivalIndex) < req_flybyRadius
                    continue
                end

                % ------------------------------------
                %%% Store solution
                % ------------------------------------
                solutionIndex_1 = solutionIndex_1 + 1;
                
                %%% Relevant dates
                candidateStruct_leg1_leg2{solutionIndex_1}.JD_Launch       = JDi_Launch;
                candidateStruct_leg1_leg2{solutionIndex_1}.JD_VGAArrival   = JDi_VGAArrival;
                candidateStruct_leg1_leg2{solutionIndex_1}.JD_VGADeparture = JDi_VGADeparture;
                candidateStruct_leg1_leg2{solutionIndex_1}.JD_EGA1Arrival  = JDi_EGA1Arrival;
                candidateStruct_leg1_leg2{solutionIndex_1}.tof_days        = JDi_EGA1Arrival-JDi_Launch;
                
                %%% Relevant values
                candidateStruct_leg1_leg2{solutionIndex_1}.c3_earth_p     = c3s_Launch_p(LaunchIndex, venusArrivalIndex);
                candidateStruct_leg1_leg2{solutionIndex_1}.vInf_VGA_m     = VInfs_VGA_m(LaunchIndex,venusArrivalIndex);
                candidateStruct_leg1_leg2{solutionIndex_1}.vInf_VGA_p     = VInfs_VGA_p(VGADepartureIndex,EGA1ArrivalIndex);
                candidateStruct_leg1_leg2{solutionIndex_1}.vInf_VGA_diff  = abs(VInfs_VGA_p(VGADepartureIndex,EGA1ArrivalIndex) - VInfs_VGA_m(LaunchIndex,venusArrivalIndex));
                candidateStruct_leg1_leg2{solutionIndex_1}.flybyRadii_VGA = flybyRadii_VGA(VGADepartureIndex, EGA1ArrivalIndex);
                candidateStruct_leg1_leg2{solutionIndex_1}.vInf_EGA1_m    = VInfs_EGA1_m(VGADepartureIndex,EGA1ArrivalIndex);
                candidateStruct_leg1_leg2{solutionIndex_1}.leg1_X0        = X0s_sc_leg1(LaunchIndex,venusArrivalIndex);
                candidateStruct_leg1_leg2{solutionIndex_1}.leg2_X0        = X0s_sc_leg2(VGADepartureIndex,EGA1ArrivalIndex);

                
            end
        end
        
        
    end
end

% ------------------------------------
%%% Reassigning Data
% ------------------------------------
%%% Deleting empty cells
if length(candidateStruct_leg1_leg2(~cellfun('isempty',candidateStruct_leg1_leg2))) == length(candidateStruct_leg1_leg2)
    warning('Need to preallocate more space to JDs_Launch_JupArr_JupDep_PluArr')
end
candidateStruct_leg1_leg2 = candidateStruct_leg1_leg2(~cellfun('isempty',candidateStruct_leg1_leg2));

%%% Reassigning results to get them out of cells
clear candidates_leg1_leg2
candidates_leg1_leg2 = [candidateStruct_leg1_leg2{:}];

candidates_JD_Launch         = [candidates_leg1_leg2.JD_Launch];
candidates_JD_VenusArrival   = [candidates_leg1_leg2.JD_VGAArrival];
candidates_JD_VenusDeparture = [candidates_leg1_leg2.JD_VGADeparture];
candidates_JD_EGA1Arrival    = [candidates_leg1_leg2.JD_EGA1Arrival];
candidates_tof_days          = [candidates_leg1_leg2.tof_days];
candidates_c3_earth_p        = [candidates_leg1_leg2.c3_earth_p];
candidates_vInf_VGA_m        = [candidates_leg1_leg2.vInf_VGA_m];
candidates_vInf_VGA_p        = [candidates_leg1_leg2.vInf_VGA_p];
candidates_vInf_VGA_diff     = [candidates_leg1_leg2.vInf_VGA_diff];
candidates_flybyRadii        = [candidates_leg1_leg2.flybyRadii_VGA];
candidates_vInf_EGA1_m       = [candidates_leg1_leg2.vInf_EGA1_m];

898

myIndxs1 = find(candidates_JD_EGA1Arrival == max(candidates_JD_EGA1Arrival)); % Getting furthest right solution line on next plot (most favorable vInfs at both ends)
myIndex1 = find(candidates_vInf_VGA_diff == min(candidates_vInf_VGA_diff(myIndxs1))); % Finding min vInf-diff solution on that line

% myIndex = 219; % Picked based on smalled vInfDiff and other prms were fine
JD_Launch_new = candidates_JD_Launch(myIndex1);
JD_VGA_new    = candidates_JD_VenusArrival(myIndex1);
JD_EGA1_new   = candidates_JD_EGA1Arrival(myIndex1);
JD_EGA2_new   = JD_EGA1_new + 365.242189*2;
% ------------------------------------
%%% Transfer TOFs
% ------------------------------------
%%% TOFs
tof_Launch2VGA = (JD_VGA_new - JD_Launch_new)*86400; % sec
tof_VGA2EGA1 = (JD_EGA1_new - JD_VGA_new)*86400;
tof_EGA12EGA2 = (JD_EGA2_new - JD_EGA1_new)*86400;

% ------------------------------------
%%% Ephemeris positions of planets
% ------------------------------------
%%% Earth state at Launch
[rEarth_Launch, vEarth_Launch] = JulianDate_2_HeliocentricState(JD_Launch_new, 'Earth'); % km, km/s

%%% Venus state at VGA
[rVenus_VGA, vVenus_VGA] = JulianDate_2_HeliocentricState(JD_VGA_new, 'Venus'); % km, km/s

%%% Earth state at EGA1
[rEarth_EGA1, vEarth_EGA1] = JulianDate_2_HeliocentricState(JD_EGA1_new, 'Earth'); % km, km/s

%%% Earth state at EGA2
[rEarth_EGA2, vEarth_EGA2] = JulianDate_2_HeliocentricState(JD_EGA2_new, 'Earth'); % km, km/s

%%% Saturn state at EGA2
[rSaturn_EGA2, vSaturn_EGA2] = JulianDate_2_HeliocentricState(JD_EGA2_new, 'Saturn'); % km, km/s
% ------------------------------------
%%% Lambert solutions
% ------------------------------------
%%% Launch to VGA
[vSc_Launch, vSc_m_VGA, ~] = lambertSolver(rEarth_Launch, rVenus_VGA, tof_Launch2VGA, 0, 0, Sun.u); % km/s

%%% VGA to EGA1
[vSc_p_VGA, vSc_m_EGA1, ~] = lambertSolver(rVenus_VGA, rEarth_EGA1, tof_VGA2EGA1, 0, 0, Sun.u); % km/s

%%% EGA1 to EGA2
[vSc_p_EGA1_nom, vSc_m_EGA2_nom, ~] = lambertSolver(rEarth_EGA1, rEarth_EGA2, tof_EGA12EGA2, 0, 0, Sun.u); % km/s

% ------------------------------------
%%% EGA1 guess conditions
% ------------------------------------
%%% vInf
vInf_m_EGA1 = vSc_m_EGA1 - vEarth_EGA1; % km/s

VInf_m_EGA1 = norm(vInf_m_EGA1); % km/s

end % run_candidateSearch_leg1_leg2 == 1


% % % % % % % EGA2 at current max candidate (~350.7 days past 6/24/2027
% % % % % % % SOI at ~674.7 days past 3/11/2032
% % % % % % % VInf from EGA2 around 13 km/s
% % % % % % % VInf at SOI around 5.9


%% =======================================================================
%%% Porkchop plot to figure out desired conditions after EGA2
% ========================================================================
if run_study_guess_EGA2_2_SOI == 1
% ------------------------------------
%%% Setting up julian dates
% ------------------------------------
%%% Updating guesses
JD_SOI_guess = JD_EGA2_new + (1.920542535961106e+08)/86400;

%%% Setting up range of JDs
% JD_Departures_EGA2_vec = linspace(JD_EGA2_new-10, JD_EGA2_new+10, nPoints);
JD_Departures_EGA2_vec = linspace(JD_EGA2_new-230, JD_EGA2_new+180, nPoints);
JD_Arrivals_SOI_vec  = linspace(JD_SOI_guess-2*365.242189, JD_SOI_guess+2*365.242189, nPoints);

%%% Inserting my chosen JD_EGA2_new
JD_Departures_EGA2_vec = sort([JD_EGA2_new, JD_Departures_EGA2_vec(1:end-1)]);
[JD_Departures_EGA2, JD_Arrivals_SOI] = meshgrid(JD_Departures_EGA2_vec, JD_Arrivals_SOI_vec);

% JD_EGA2_find(abs(JD_Departures_EGA2_vec-JD_EGA2_new) == min(abs(JD_Departures_EGA2_vec-JD_EGA2_new)))
% ------------------------------------
%%% Looping through options and computing vInf and c3
% ------------------------------------
X0s_sc_leg4{size(JD_Departures_EGA2,1),size(JD_Departures_EGA2,2)} = [];
VInfs_SOI_m = zeros(size(JD_Departures_EGA2));
VInfs_EGA2_p = zeros(size(JD_Departures_EGA2));
% flybyRadii_EGA2 = zeros(size(JD_Departures_EGA2));
TOFs_days_leg4 = zeros(size(JD_Departures_EGA2));

%%% Preallocating solutions for my EGA2
solutionIndex_2 = 0;
candidateStruct_leg4{100000} = [];

for xIndex = 1:nPoints
    for yIndex = 1:nPoints
        %%% Grab current launch window
        JD_D = JD_Departures_EGA2(xIndex, yIndex);
        JD_A = JD_Arrivals_SOI(xIndex, yIndex);
        
        %%% Find Venus position at departure time
        [rEarth_D, vEarth_D] = JulianDate_2_HeliocentricState(JD_D, 'Earth');
        
        %%% Find Saturn position at arrival time
        [rSaturn_A, vSaturn_A] = JulianDate_2_HeliocentricState(JD_A, 'Saturn');
        
        %%% Use Lambert solver to find solution
        tof_sec = (JD_A - JD_D)*86400;
        [vSc_EGA2Departure, vSc_SOIArrival,~] = lambertSolver(rEarth_D,rSaturn_A,tof_sec,0,0,Sun.u);
        
        %%% vInf at Earth departure
        vInf_EGA2_p = vSc_EGA2Departure - vEarth_D;
        VInfs_EGA2_p(xIndex, yIndex) = norm(vInf_EGA2_p);
        
%         %%% flyby radius at EGA2
%         vInf_Earth_m = vInfs_EGA2_m{xIndex, yIndex};
%         VInf_Earth_m = (norm(vInf_Earth_m) + norm(vInf_Earth_p))/2;
%         turningAngle_rad = acos(dot(vInf_Earth_m,vInf_Earth_p)./(norm(vInf_Earth_m)*norm(vInf_Earth_p)));
%         [ rp ] = calcFlybyRadius( turningAngle_rad, VInf_Earth_m, Earth.u );
%         flybyRadii_EGA2(xIndex, yIndex) = rp;
        
        %%% vInf at Saturn arrival
        vInf_SOI_m = vSc_SOIArrival - vSaturn_A;
        VInfs_SOI_m(xIndex, yIndex) = norm(vInf_SOI_m);
        
        %%% TOFs
        TOFs_days_leg4(xIndex, yIndex) = JD_A - JD_D;
        
        %%% Starting state for this leg of the trajectory [r0; v0]
        X0s_sc_leg4{xIndex, yIndex} = [rEarth_D; vSc_EGA2Departure];
        
        if JD_D == JD_EGA2_new
            solutionIndex_2 = solutionIndex_2 + 1;
            
            candidateStruct_leg4{solutionIndex_2}.cand_JD_SOI            = JD_A;
            candidateStruct_leg4{solutionIndex_2}.cand_tof_sec           = tof_sec;
            candidateStruct_leg4{solutionIndex_2}.cand_vSc_EGA2Departure = vSc_EGA2Departure;
            candidateStruct_leg4{solutionIndex_2}.cand_vSc_SOIArrival    = vSc_SOIArrival;
            candidateStruct_leg4{solutionIndex_2}.cand_vInf_EGA2_p       = vInf_EGA2_p;
            candidateStruct_leg4{solutionIndex_2}.cand_VInf_EGA2_p       = norm(vInf_EGA2_p);
            candidateStruct_leg4{solutionIndex_2}.cand_vInf_SOI_m        = vInf_SOI_m;
            candidateStruct_leg4{solutionIndex_2}.cand_VInf_SOI_m        = norm(vInf_SOI_m);
            candidateStruct_leg4{solutionIndex_2}.cand_X0_sc_leg4        = X0s_sc_leg4{xIndex, yIndex};
            
        end
    end
end

% ------------------------------------
%%% Cleaning up candidates and selecting
% ------------------------------------
candidateStruct_leg4 = candidateStruct_leg4(~cellfun('isempty',candidateStruct_leg4));

%%% Reassigning results to get them out of cells
clear candidates_leg1_leg2
candidates_leg4              = [candidateStruct_leg4{:}];

candidates2_JD_SOI            = [candidates_leg4.cand_JD_SOI];
candidates2_tof_sec           = [candidates_leg4.cand_tof_sec];
candidates2_vSc_EGA2Departure = [candidates_leg4.cand_vSc_EGA2Departure];
candidates2_vSc_SOIArrival    = [candidates_leg4.cand_vSc_SOIArrival];
candidates2_vInf_EGA2_p       = [candidates_leg4.cand_vInf_EGA2_p];
candidates2_VInf_EGA2_p       = [candidates_leg4.cand_VInf_EGA2_p];
candidates2_vInf_SOI_m        = [candidates_leg4.cand_vInf_SOI_m];
candidates2_VInf_SOI_m        = [candidates_leg4.cand_VInf_SOI_m];
candidates2_X0_sc_leg4        = [candidates_leg4.cand_X0_sc_leg4];

minVInf_EGA2_indx = find(candidates2_VInf_EGA2_p == min(candidates2_VInf_EGA2_p));
minVInf_SOI_indx  = find(candidates2_VInf_SOI_m == min(candidates2_VInf_SOI_m));

myIndx2 = minVInf_SOI_indx;

JD_SOI_new = candidates2_JD_SOI(myIndx2);
tof_EGA22SOI_days = JD_SOI_new - JD_EGA2_new;
tof_EGA22SOI_sec = tof_EGA22SOI_days*86400;

% ------------------------------------
%%% Pork Chop Plot
% ------------------------------------
figure('position',[440 234 762 564]); hold all
vInf_EGA2_contours = [10, 11, 13, 15, 20, 25, 30, 40, 50];
[cs1, h1] = contour(JD_Departures_EGA2-JD_Departures_EGA2(1),JD_Arrivals_SOI-JD_Arrivals_SOI(1),VInfs_EGA2_p,vInf_EGA2_contours,'ShowText','on','color',colors.std.blue,'linewidth',1.2);
% [cs1, h1] = contour(JD_Departures_EGA2-JD_Departures_EGA2(1),JD_Arrivals_SOI-JD_Arrivals_SOI(1),VInfs_EGA2_p,'ShowText','on','color',colors.std.blue,'linewidth',1.2);

vInf_SOI_contours = [1, 2, 3, 4, 5,5.8, 6, 7, 8, 9];
[cs2, h2] = contour(JD_Departures_EGA2-JD_Departures_EGA2(1),JD_Arrivals_SOI-JD_Arrivals_SOI(1),VInfs_SOI_m,vInf_SOI_contours,'ShowText','on','color',colors.std.ltred,'linewidth',1.2);
% [cs2, h2] = contour(JD_Departures_EGA2-JD_Departures_EGA2(1),JD_Arrivals_SOI-JD_Arrivals_SOI(1),VInfs_SOI_m,'ShowText','on','color',colors.std.ltred,'linewidth',1.2);

TOF_contours_leg4 = [260, 280, 300, 320, 340];
% [cs3,h3] = contour(JD_Departures_EGA2-JD_Departures_EGA2(1),JD_Arrivals_SOI-JD_Arrivals_SOI(1),TOFs_days_leg4,TOF_contours_leg4,'ShowText','on','color',colors.std.black,'linewidth',1.2);
[cs3,h3] = contour(JD_Departures_EGA2-JD_Departures_EGA2(1),JD_Arrivals_SOI-JD_Arrivals_SOI(1),TOFs_days_leg4,'ShowText','on','color',colors.std.black,'linewidth',1.2);

[yr, mo, d, ~, ~, ~] = julian2calendar(JD_Departures_EGA2_vec(1));
calendar0_EGA2 = sprintf('EGA2, Days Past %d/%d/%d',mo,d,yr);

[yr, mo, d, ~, ~, ~] = julian2calendar(JD_Arrivals_SOI_vec(1));
calendar0_SOI = sprintf('SOI, Days Past %d/%d/%d',mo,d,yr);

%%% Plotting line resulting from current EGA2 choice
plot([1,1].*(JD_EGA2_new - JD_Departures_EGA2_vec(1)),[JD_Arrivals_SOI_vec(1), JD_Arrivals_SOI_vec(end)]-JD_Arrivals_SOI_vec(1),'m') 
% plot([1,1].*(min(candidates_JD_EGA1Arrival)+365.242189*2 - JD_Departures_EGA2_vec(1)),[JD_Arrivals_SOI_vec(1), JD_Arrivals_SOI_vec(end)]-JD_Arrivals_SOI_vec(1),'m','linewidth',2) 
% plot([1,1].*(max(candidates_JD_EGA1Arrival)+365.242189*2 - JD_Departures_EGA2_vec(1)),[JD_Arrivals_SOI_vec(1), JD_Arrivals_SOI_vec(end)]-JD_Arrivals_SOI_vec(1),'m','linewidth',2) 

%%% Plotting current choice
plot(JD_EGA2_new - JD_Departures_EGA2_vec(1),JD_SOI_new-JD_Arrivals_SOI_vec(1),'m.','markersize',16) 


PlotBoi2(calendar0_EGA2,calendar0_SOI,16)
grid off
clabel(cs1 ,h1,'color',colors.std.blue)
clabel(cs2 ,h2,'color',colors.std.ltred)
clabel(cs3 ,h3,'color',colors.std.black)

legend([h1, h2, h3,],'V_\infty @ EGA2, km/s','V_\infty @ SOI, km/s','TOF, days')

end % run_study_guess_EGA2_2_SOI == 1




% ========================================================================
%%% Brief study - looking at current candidates to find a solution that
%%% minimizes the discrepency in VInf_EGA2 plus and minus
% ========================================================================
DVInfs_EGA2 = zeros(length(candidates_leg4),1);
for jj = 1:length(candidates_leg4)
%     candidates2_JD_SOI            = [candidates_leg4.cand_JD_SOI];
%     candidates2_tof_sec           = [candidates_leg4.cand_tof_sec];
%     candidates2_vSc_EGA2Departure = [candidates_leg4.cand_vSc_EGA2Departure];
%     candidates2_vSc_SOIArrival    = [candidates_leg4.cand_vSc_SOIArrival];
%     candidates2_vInf_EGA2_p       = [candidates_leg4.cand_vInf_EGA2_p];
%     candidates2_VInf_EGA2_p       = [candidates_leg4.cand_VInf_EGA2_p];
%     candidates2_vInf_SOI_m        = [candidates_leg4.cand_vInf_SOI_m];
%     candidates2_VInf_SOI_m        = [candidates_leg4.cand_VInf_SOI_m];
%     candidates2_X0_sc_leg4        = [candidates_leg4.cand_X0_sc_leg4];
% 
%     minVInf_EGA2_indx = find(candidates2_VInf_EGA2_p == min(candidates2_VInf_EGA2_p));
%     minVInf_SOI_indx  = find(candidates2_VInf_SOI_m == min(candidates2_VInf_SOI_m));
% 
%     myIndx2 = minVInf_SOI_indx;
% 
%     JD_SOI_new = candidates2_JD_SOI(myIndx2);
%     tof_EGA22SOI_days = JD_SOI_new - JD_EGA2_new;
%     tof_EGA22SOI_sec = tof_EGA22SOI_days*86400;
    JD_SOI_jj = candidates2_JD_SOI(jj);
    tof_EGA22SOI_sec_jj = (JD_SOI_jj - JD_EGA2_new)*86400;
    
    %%% Saturn state at EGA2
    [rSaturn_SOI_jj, vSaturn_SOI_jj] = JulianDate_2_HeliocentricState(JD_SOI_jj, 'Saturn'); % km, km/s

    %%% EGA2 to SOI
    [vSc_p_EGA2_nom_jj, vSc_m_SOI_nom_kk, ~] = lambertSolver(rEarth_EGA2, rSaturn_SOI_jj, tof_EGA22SOI_sec_jj, 0, 0, Sun.u); % km/s

    %%% EGA2 guess conditions
    vInf_p_EGA2_jj = vSc_p_EGA2_nom_jj - vEarth_EGA2;
    VInf_p_EGA2_jj = norm(vInf_p_EGA2_jj);
    
    % First Flyby
    VSc_p_EGA1_jj = norm(vSc_p_EGA1_nom);
    VEarth_EGA1_jj = norm(vEarth_EGA1);
    theta_EGA1_rad_jj = acos((-VSc_p_EGA1_jj^2 + VInf_m_EGA1^2 + VEarth_EGA1_jj^2)/(2*VInf_m_EGA1*VEarth_EGA1_jj));

    phi_rad_jj = 1.523196438104142;
    
    %%% Compute v-inf plus
    vInf_EGA1_p_VNC_jj = VInf_m_EGA1 * [cos(pi-theta_EGA1_rad_jj); sin(pi-theta_EGA1_rad_jj)*cos(phi_rad_jj); -sin(pi-theta_EGA1_rad_jj)*sin(phi_rad_jj)];
    
    
    
    
    
    
    
    %%% Rotate from VNC to inertial
    vHat_EGA1_jj = vEarth_EGA1 / VEarth_EGA1_jj;
    hHat_EGA1_jj = cross(rEarth_EGA1,vEarth_EGA1)/norm(cross(rEarth_EGA1,vEarth_EGA1));
    nHat_EGA1_jj = hHat_EGA1_jj;
    cHat_EGA1_jj = cross(vHat_EGA1_jj,nHat_EGA1_jj);
    T_VNC2Ecliptic_EGA1_jj = [vHat_EGA1_jj, nHat_EGA1_jj, cHat_EGA1_jj];
    vInf_p_EGA1_jj = T_VNC2Ecliptic_EGA1_jj * vInf_EGA1_p_VNC_jj;
    
    vInf_m_EGA2_jj = vInf_p_EGA1_jj + vEarth_EGA1 - vEarth_EGA2;
    
    VInf_m_EGA2 = norm(vInf_m_EGA2_jj)
    
    DVInf_EGA2 = VInf_m_EGA2 - VInf_p_EGA2_jj;
    DVInfs_EGA2(jj) = DVInf_EGA2;
    
end

989
warning('Have to change EGA2 ... EGA1 .... At all these candidates, VInf_m_EGA2 is 10.6, which is not possible for VInf_p_EGA2 on current EGA2-SOI PCP')
return






























%% =======================================================================
%%% Resonant trajectory (EGA1 to EGA2)
% ========================================================================
if run_study_EGA1_2_EGA2 == 1
% ------------------------------------
%%% Nominal Lambert study of leg 4
% ------------------------------------
% %%% EGA1 to EGA2
% [vSc_p_EGA1_nom, vSc_m_EGA2_nom, ~] = lambertSolver(rEarth_EGA1, rEarth_EGA2, tof_EGA12EGA2, 0, 0, Sun.u); % km/s
% %%% Earth state at EGA2
% [rEarth_EGA2, vEarth_EGA2] = JulianDate_2_HeliocentricState(JD_EGA2_new, 'Earth'); % km, km/s

%%% Saturn state at EGA2
[rSaturn_SOI, vSaturn_SOI] = JulianDate_2_HeliocentricState(JD_SOI_new, 'Saturn'); % km, km/s

%%% EGA2 to SOI
[vSc_p_EGA2_nom, vSc_m_SOI_nom, ~] = lambertSolver(rEarth_EGA2, rSaturn_SOI, tof_EGA22SOI_sec, 0, 0, Sun.u); % km/s

%%% EGA2 guess conditions
vInf_p_EGA2 = vSc_p_EGA2_nom - vEarth_EGA2;
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
% 
% % Second Flyby
% VSc_EGA2_p = norm(vSc_p_EGA2);
% VEarth_EGA2 = norm(vEarth_EGA2);

%%% Loop through phi values
rps_EGA1 = zeros(n_phis,1);
rps_EGA2 = zeros(n_phis,1);
vInfs_p_EGA1 = zeros(n_phis,3);
vInfs_m_EGA2 = zeros(n_phis,3);

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
p4 = fill([67.093, 67.093, 117.33, 117.33], [0, 13000,13000,0], 'g');
p1 = plot(phis_rad.*(180/pi),rps_EGA1,'linewidth',2,'color',colors.std.blue);
p2 = plot(phis_rad.*(180/pi),rps_EGA2,'linewidth',2,'color',colors.std.ltred);
% % p3 = plot([0,360],[Earth.R, Earth.R],'--k','linewidth',1);
p3 = plot([0,360],[Earth.R+300, Earth.R+300],'--k','linewidth',1);
p5 = plot([myPhi_deg,myPhi_deg],[0, 13000],'m','linewidth',1.5);
PlotBoi2('$\Phi$, $^\circ$','Flyby Radius, $km$',18,'LaTex')
xlim([0 360])
ylim([0 13000])
% legend([p1 p2 p3],'EGA1','EGA2','300 km Altitude')
legend([p1 p2 p3 p4 p5],'EGA1','EGA2','300 km Altitude','Range of Acceptable \Phi', 'Selected \Phi')
% 
% % ------------------------------------
% %%% Wrapping up EGA1 and EGA2
% % ------------------------------------
% %%% Ensuring that VInfs are close
DVInf_EGA1 = VInf_m_EGA1 - VInf_p_EGA1;
DVInf_EGA2 = VInf_m_EGA2 - VInf_p_EGA2;
% 
% %%% Flyby radius
% flybyRadius_EGA1 = rps_EGA1(myPhi_index); % km
% flybyRadius_EGA2 = rps_EGA2(myPhi_index); % km
% 
% %%% B-Plane
% [BT_EGA1, BR_EGA1, B_EGA1, theta_EGA1_rad] = BPlaneTargeter_vInfs(vInf_m_EGA1, vInf_p_EGA1, Earth.u, 1);
% [BT_EGA2, BR_EGA2, B_EGA2, theta_EGA2_rad] = BPlaneTargeter_vInfs(vInf_m_EGA2, vInf_p_EGA2, Earth.u, 1);

end % run_study_EGA1_2_EGA2 == 1





% %% =======================================================================
% %%% Pork Chop Plot - EGA2 to SOI (leg4)
% % ========================================================================
% if run_study_EGA2_2_SOI == 1
% % ------------------------------------
% %%% Setting up julian dates
% % ------------------------------------
% JD_Departures_EGA2_vec = linspace(JD_guess_EGA2-30, JD_guess_EGA2+30, nPoints);
% JD_Arrivals_SOI_vec  = linspace(JD_guess_SOI-30, JD_guess_SOI+30, nPoints);
% 
% [JD_Departures_EGA2, JD_Arrivals_SOI] = meshgrid(JD_Departures_EGA2_vec, JD_Arrivals_SOI_vec);
% 
% % ------------------------------------
% %%% Looping through options and computing vInf and c3
% % ------------------------------------
% X0_sc_leg3{size(JD_Departures_EGA2,1),size(JD_Departures_EGA2,2)} = [];
% VInfs_SOI_m = zeros(size(JD_Departures_EGA2));
% VInfs_EGA2_p = zeros(size(JD_Departures_EGA2));
% % flybyRadii_EGA2 = zeros(size(JD_Departures_EGA2));
% TOFs_days = zeros(size(JD_Departures_EGA2));
% 
% for xIndex = 1:nPoints
%     for yIndex = 1:nPoints
%         %%% Grab current launch window
%         JD_D = JD_Departures_EGA2(xIndex, yIndex);
%         JD_A = JD_Arrivals_SOI(xIndex, yIndex);
%         
%         %%% Find Earth position at departure time
%         [rEarth_D, vEarth_D] = JulianDate_2_HeliocentricState(JD_D, 'Earth');
%         
%         %%% Find Earth position at arrival time
%         [rSaturn_A, vSaturn_A] = JulianDate_2_HeliocentricState(JD_A, 'Saturn');
%         
%         %%% Use Lambert solver to find solution
%         tof_sec = (JD_A - JD_D)*86400;
%         [vSc_EarthDeparture, vSc_SaturnArrival,~] = lambertSolver(rEarth_D,rSaturn_A,tof_sec,0,0,Sun.u);
%         
%         %%% vInf at Earth departure
%         vInf_Earth_p = vSc_EarthDeparture - vEarth_D;
%         VInfs_EGA2_p(xIndex, yIndex) = norm(vInf_Earth_p);
%         
% %         %%% flyby radius at Earth
% %         vInf_Earth_m = vInfs_EGA2_m{xIndex, yIndex};
% %         VInf = (norm(vInf_Earth_m) + norm(vInf_Earth_p))/2;
% %         turningAngle_rad = acos(dot(vInf_Earth_m,vInf_Earth_p)./(norm(vInf_Earth_m)*norm(vInf_Earth_p)));
% %         [ rp ] = calcFlybyRadius( turningAngle_rad, VInf, Venus.u );
% %         flybyRadii_EGA2(xIndex, yIndex) = rp;
% %         
%         %%% vInf at Saturn arrival
%         vInf_Saturn = vSc_SaturnArrival - vSaturn_A;
%         VInfs_SOI_m(xIndex, yIndex) = norm(vInf_Saturn);
%         
%         %%% TOFs
%         TOFs_days(xIndex, yIndex) = JD_A - JD_D;
%         
%         %%% Starting state for this leg of the trajectory [r0; v0]
%         X0_sc_leg3{xIndex, yIndex} = [rEarth_D; vSc_EarthDeparture];
% 
%     end
% end
% 
% 
% % ------------------------------------
% %%% Pork Chop Plot
% % ------------------------------------
% figure('position',[440 234 762 564]); hold all
% vInf_EGA2_contours = [8.8, 9, 9.5, 10, 11, 12, 20];
% [cs1, h1] = contour(JD_Departures_EGA2-JD_Departures_EGA2(1),JD_Arrivals_SOI-JD_Arrivals_SOI(1),VInfs_EGA2_p,vInf_EGA2_contours,'ShowText','on','color',colors.std.blue,'linewidth',1.2);
% 
% vInf_SOI_contours = [5.8, 5.85, 5.9, 6, 7];
% [cs2, h2] = contour(JD_Departures_EGA2-JD_Departures_EGA2(1),JD_Arrivals_SOI-JD_Arrivals_SOI(1),VInfs_SOI_m,vInf_SOI_contours,'ShowText','on','color',colors.std.ltred,'linewidth',1.2);
% % 
% TOF_contours = [1160, 1180, 1200, 1220, 1240];
% [cs3,h3] = contour(JD_Departures_EGA2-JD_Departures_EGA2(1),JD_Arrivals_SOI-JD_Arrivals_SOI(1),TOFs_days,TOF_contours,'ShowText','on','color',colors.std.black,'linewidth',1.2);
% 
% p1 = plot(30, 30,'o','markeredgecolor',colors.std.black,'markerfacecolor',colors.std.cyan,'markersize',10);
% 
% [yr, mo, d, hr, min, sec] = julian2calendar(JD_Departures_EGA2_vec(1));
% calendar0_EGA2 = sprintf('Days Past %d/%d/%d',mo,d,yr);
% 
% [yr, mo, d, hr, min, sec] = julian2calendar(JD_Arrivals_SOI_vec(1));
% calendar0_SOI = sprintf('Days Past %d/%d/%d',mo,d,yr);
% 
% PlotBoi2(calendar0_EGA2,calendar0_SOI,16)
% grid off
% clabel(cs1 ,h1,'color',colors.std.blue)
% clabel(cs2 ,h2,'color',colors.std.ltred)
% clabel(cs3 ,h3,'color',colors.std.black)
% 
% legend([h1, h2, h3],'V_\infty @ EGA2, km/s','V_\infty @ SOI, km/s','TOF, days')
% 
% 
% end % run_study_EGA2_2_SOI == 1

































