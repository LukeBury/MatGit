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

%% =======================================================================
%%% Setup
% ========================================================================
% ------------------------------------
%%% Bodies
% ------------------------------------
Sun       = bodies.sun;
Venus     = bodies.venus;
Earth     = bodies.earth;
Saturn    = bodies.saturn;
Mimas     = bodies.mimas;
Titan     = bodies.titan;
Enceladus = bodies.enceladus;
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
%%% Julian Date guess ranges
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

%% =======================================================================
%%% Set requirements
% ========================================================================
%%% C3 requirement on Earth Launch
req_c3_Launch_ub = 25; % km^2/s^2

%%% Tolerance on abs(|vInf_m| - |vInf_p|) at VGA
req_VInf_VGA_diff = 0.1; % km/s

%%% Flyby radius for Venus
req_altitude_VGA = 250; % km
req_flybyRadius_VGA_lb = Venus.R + req_altitude_VGA; % km

%%% Upper bound of VInf_SOI
req_vInf_SOI_ub = 9; % km/s

%%% Upper bound of total travel time
req_TOF_yrs_ub = 15;

%%% Tolerance on abs(|vInf_m| - |vInf_p|) at EGA2
req_VInf_EGA2_diff = 0.1; % km/s

%%% Earth flyby altitude lower bound
req_altitude_EGA = 300; % km
req_flybyRadius_EGA_lb = Earth.R + req_altitude_EGA; % km

%% =======================================================================
%%% Loop through dates
% ========================================================================

% ------------------------------------
%%% Preparing to loop through Launch dates
% ------------------------------------
solutionIndex = 0;
solutionsStruc{10000000} = [];

sum1 = 0;
sum2 = 0;
minValue = 100000000;
% ------------------------------------
%%% Loop through Launch dates
% ------------------------------------
for k = 1:length(JDs_Launch)
    JD_Launch = JDs_Launch(k);
    
    % ------------------------------------
    %%% Loop through VGA dates
    % ------------------------------------
    for j = 1:length(JDs_VGA)
        JD_VGA = JDs_VGA(j);
        
        %%% Make sure VGA is after Launch
        if (JD_VGA - JD_Launch) < 10
            continue
        end
        
        %%% Get states of planets
        [rEarth_Launch, vEarth_Launch] = JulianDate_2_HeliocentricState(JD_Launch, 'Earth');
        [rVenus_VGA, vVenus_VGA]       = JulianDate_2_HeliocentricState(JD_VGA, 'Venus');

        %%% Lamberts between planets
        [vSc_Launch, vSc_VGA_m,~] = lambertSolver(rEarth_Launch,rVenus_VGA,(JD_VGA - JD_Launch)*86400,0,0,Sun.u);
        
        %%% Determine Launch c3
        vInf_Launch = vSc_Launch - vEarth_Launch; % km/s
        VInf_Launch = norm(vInf_Launch);
        c3_Launch = norm(vInf_Launch)^2; % km^2/s^2
        
        %%% Check c3_Launch requirement
        if c3_Launch > req_c3_Launch_ub
            continue
        end
        
        % ------------------------------------
        %%% Loop through EGA dates
        % ------------------------------------
        for i = 1:length(JDs_EGA1)
            JD_EGA1 = JDs_EGA1(i);
            JD_EGA2 = JDs_EGA2(i);
            
            %%% Make sure EGA1 is after VGA
            if (JD_EGA1 - JD_VGA) < 40
                continue
            end
            
            %%% Get states of planets
            [rEarth_EGA1, vEarth_EGA1] = JulianDate_2_HeliocentricState(JD_EGA1, 'Earth');

            %%% Lamberts between planets
            [vSc_VGA_p, vSc_EGA1_m,~] = lambertSolver(rVenus_VGA,rEarth_EGA1,(JD_EGA1 - JD_VGA)*86400,0,0,Sun.u);
        
            %%% Check VInf_VGA_diff requirement
            vInf_VGA_m = vSc_VGA_m - vVenus_VGA; % km/s
            vInf_VGA_p = vSc_VGA_p - vVenus_VGA; % km/s
            
            VInf_VGA_m = norm(vInf_VGA_m); % km/s
            VInf_VGA_p = norm(vInf_VGA_p); % km/s
            
            if abs(VInf_VGA_p - VInf_VGA_m) > req_VInf_VGA_diff
                continue
            end
            
            %%% Check Venus flyby altitude requirement            
            turningAngle_VGA_rad = acos(dot(vInf_VGA_m,vInf_VGA_p)./(VInf_VGA_m*VInf_VGA_p));
            [ flybyRadius_VGA ] = calcFlybyRadius( turningAngle_VGA_rad, (VInf_VGA_m+VInf_VGA_p)/2, Venus.u );
            
            if flybyRadius_VGA < req_flybyRadius_VGA_lb
                continue
            end
            
            % ------------------------------------
            %%% Loop through SOI dates
            % ------------------------------------
            for q = 1:length(JDs_SOI)
                JD_SOI = JDs_SOI(q);
                
                %%% Make sure SOI is after EGA2
                if (JD_SOI - JD_EGA2) < 365.25
                    continue
                end

                %%% Check TOF requirement
                if (JD_SOI - JD_Launch)/days_per_year > req_TOF_yrs_ub
                    continue
                end
                
                %%% Get states of planets
                [rEarth_EGA2, vEarth_EGA2] = JulianDate_2_HeliocentricState(JD_EGA2, 'Earth');
                [rSaturn_SOI, vSaturn_SOI] = JulianDate_2_HeliocentricState(JD_SOI, 'Saturn');
                
                %%% Lamberts between planets
                [vSc_EGA2_p, vSc_SOI_m,~] = lambertSolver(rEarth_EGA2,rSaturn_SOI,(JD_SOI - JD_EGA2)*86400,0,0,Sun.u);
                
                %%% Check VInf_SOI requirement
                vInf_SOI_m = vSc_SOI_m - vSaturn_SOI; % km/s
                VInf_SOI_m = norm(vInf_SOI_m);        % km/s
                
                if VInf_SOI_m > req_vInf_SOI_ub
                    continue
                end
                
                % ------------------------------------
                %%% Study resonant orbit
                % ------------------------------------
                %%% Lambert estimate from EGA1 to EGA2
                [vSc_EGA1_p, vSc_EGA2_m, ~] = lambertSolver(rEarth_EGA1, rEarth_EGA2, resonance_period_secs, 0, 0, Sun.u); % km/s
                
                %%% Necessary pre-determined velocity values for resonance study
                vInf_EGA1_m = vSc_EGA1_m - vEarth_EGA1; % km/s
                VInf_EGA1_m = norm(vInf_EGA1_m);        % km/s
                
                vInf_EGA2_p = vSc_EGA2_p - vEarth_EGA2; % km/s
                VInf_EGA2_p = norm(vInf_EGA2_p);        % km/s
                
                VSc_EGA1_p = norm(vSc_EGA1_p);   % km/s
                VEarth_EGA1 = norm(vEarth_EGA1); % km/s
                theta_EGA1_rad = acos((-VSc_EGA1_p^2 + VInf_EGA1_m^2 + VEarth_EGA1^2)/(2*VInf_EGA1_m*VEarth_EGA1));
                
                VSc_EGA2_p = norm(vSc_EGA2_p);   % km/s
                VEarth_EGA2 = norm(vEarth_EGA2); % km/s
                
                %%% Options for resonance study
                n_phis = 360;
                phis_rad = linspace(0,2*pi,n_phis);
                
                %%% Preallocating for this batch of phis
                rps_EGA1     = zeros(n_phis,1);
                rps_EGA2     = zeros(n_phis,1);
                smallerRps   = zeros(n_phis,1);
                vInfs_EGA1_p = zeros(n_phis,3);
                vInfs_EGA2_m = zeros(n_phis,3);
                
                %%% Looping through phi values
                for kk = 1:n_phis
                    %%% Set phi
                    phi_rad = phis_rad(kk);
                    
                    % ---------------------
                    %%% EGA 1
                    % ---------------------
                    %%% Compute v-inf plus in VNC frame
                    vInf_EGA1_p_VNC = VInf_EGA1_m * [cos(pi-theta_EGA1_rad); sin(pi-theta_EGA1_rad)*cos(phi_rad); -sin(pi-theta_EGA1_rad)*sin(phi_rad)];
                    
                    %%% Rotate from VNC to inertial
                    vHat_EGA1 = vEarth_EGA1 / VEarth_EGA1;
                    hHat_EGA1 = cross(rEarth_EGA1,vEarth_EGA1)/norm(cross(rEarth_EGA1,vEarth_EGA1));
                    nHat_EGA1 = hHat_EGA1;
                    cHat_EGA1 = cross(vHat_EGA1,nHat_EGA1);
                    T_VNC2Ecliptic_EGA1 = [vHat_EGA1, nHat_EGA1, cHat_EGA1];
                    vInf_EGA1_p = T_VNC2Ecliptic_EGA1 * vInf_EGA1_p_VNC;
                    VInf_EGA1_p = norm(vInf_EGA1_p);
                    vInfs_EGA1_p(kk,:) = vInf_EGA1_p';

                    %%% Calculate turning angle
                    turningAngle_EGA1_rad = acos(dot(vInf_EGA1_m,vInf_EGA1_p)/(VInf_EGA1_m*VInf_EGA1_p));

                    %%% Calculate EGA1 flyby radius
                    rps_EGA1(kk) = (Earth.u/(VInf_EGA1_m*VInf_EGA1_p))*((1/cos((pi-turningAngle_EGA1_rad)/2)) - 1);

                    % ---------------------
                    %%% EGA 2
                    % ---------------------
                    %%% Theta
                    vInf_EGA2_m = vInf_EGA1_p + vEarth_EGA1 - vEarth_EGA2;
                    VInf_EGA2_m = norm(vInf_EGA2_m);
                    vInfs_EGA2_m(kk,:) = vInf_EGA2_m';
                    
                    %%% Calculate turning angle
                    turningAngle_EGA2_rad = acos(dot(vInf_EGA2_m,vInf_EGA2_p)/(VInf_EGA2_m*VInf_EGA2_p));

                    %%% Calculate EGA2 flyby radius
                    rps_EGA2(kk) = (Earth.u/(VInf_EGA2_m*VInf_EGA2_p))*((1/cos((pi-turningAngle_EGA2_rad)/2)) - 1);
                    
                    %%% Storing smaller flyby radius value between EGA1 and EGA2
                    smallerRps(kk) = smallerValue(rps_EGA1(kk),rps_EGA2(kk));
                end
                

                
                %%% Grab phi for maximizing the smallest flyby radius
                myPhi_index = find(smallerRps == max(smallerRps));
                myPhi_rad = phis_rad(myPhi_index);
                myPhi_deg = myPhi_rad*180/pi;
                
                %%% Check Earth flyby altitude requirement
                if smallerRps(myPhi_index) < req_flybyRadius_EGA_lb
                    continue
                end
                
                %%% Use Phi index to get vInfs and rps for EGA1 and EGA2
                flybyRadius_EGA1 = rps_EGA1(myPhi_index);
                flybyRadius_EGA2 = rps_EGA2(myPhi_index);
                vInf_EGA1_p      = vInfs_EGA1_p(myPhi_index,:);
                vInf_EGA2_m      = vInfs_EGA2_m(myPhi_index,:);
                VInf_EGA1_p      = norm(vInf_EGA1_p);
                VInf_EGA2_m      = norm(vInf_EGA2_m);
                
                %%% Check VInf_EGA2_diff requirement
                sum1 = sum1 + 1;
                
                if abs(VInf_EGA2_p - VInf_EGA2_m) > req_VInf_EGA2_diff
                    if abs(VInf_EGA2_p - VInf_EGA2_m) < minValue
                        minValue = abs(VInf_EGA2_p - VInf_EGA2_m);
                    end
                    continue
                end
                sum2 = sum2 + 1;
                
                %%% Calculating Launch RLA and DLA
                RLA_deg = atan2(vInf_Launch(2),vInf_Launch(1))*180/pi; % deg
                DLA_deg = asin(vInf_Launch(3)/VInf_Launch)*180/pi;     % deg
                
                %%% Calculating the dv for the trans-venus-injection maneuver
                V2_TVI = sqrt(VInf_Launch^2 + 2*Earth.u/(Earth.R+300));
                DV_TVI = V2_TVI - sqrt(Earth.u / (Earth.R+300));
                
                %%% Calculating B-Plane parameters before storing
                [BT_VGA, BR_VGA, B_VGA, theta_VGA_rad] = BPlaneTargeter_vInfs(vInf_VGA_m, vInf_VGA_p, Venus.u, 1);
                [BT_EGA1, BR_EGA1, B_EGA1, theta_EGA1_rad] = BPlaneTargeter_vInfs(vInf_EGA1_m, vInf_EGA1_p, Earth.u, 1);
                [BT_EGA2, BR_EGA2, B_EGA2, theta_EGA2_rad] = BPlaneTargeter_vInfs(vInf_EGA2_m, vInf_EGA2_p, Earth.u, 1);
                
                %%% Store solution
                solutionIndex = solutionIndex + 1;
                
                if solutionIndex == 9
                    %%% Plotting
                    figure; hold all
                    p4 = fill([82.39, 82.39, 90.828, 90.828], [0, 8000,8000,0], 'g');
                    p1 = plot(phis_rad.*(180/pi),rps_EGA1,'linewidth',2,'color',colors.std.blue);
                    p2 = plot(phis_rad.*(180/pi),rps_EGA2,'linewidth',2,'color',colors.std.ltred);
                    p3 = plot([0,360],[Earth.R+300, Earth.R+300],'--k','linewidth',1);
                    p5 = plot([myPhi_deg,myPhi_deg],[0, 8000],'m','linewidth',1.5);
                    PlotBoi2('$\Phi$, $^\circ$','Flyby Radius, $km$',18,'LaTex')
                    xlim([0 360])
%                     legend([p1 p2 p3],'EGA1','EGA2','300 km Altitude')
                    legend([p1 p2 p3 p4 p5],'EGA1','EGA2','300 km Altitude','Range of Acceptable \Phi', 'Selected \Phi')

                end
                
                solutionsStruc{solutionIndex}.JD_Launch = JD_Launch;
                solutionsStruc{solutionIndex}.JD_VGA    = JD_VGA;
                solutionsStruc{solutionIndex}.JD_EGA1   = JD_EGA1;
                solutionsStruc{solutionIndex}.JD_EGA2   = JD_EGA2;
                solutionsStruc{solutionIndex}.JD_SOI    = JD_SOI;
                
                solutionsStruc{solutionIndex}.c3_Launch        = c3_Launch;
                solutionsStruc{solutionIndex}.RLA_deg          = RLA_deg;
                solutionsStruc{solutionIndex}.DLA_deg          = DLA_deg;
                solutionsStruc{solutionIndex}.DV_TVI           = DV_TVI;
                solutionsStruc{solutionIndex}.VInf_VGA_diff    = abs(VInf_VGA_p - VInf_VGA_m);
                solutionsStruc{solutionIndex}.VInf_EGA2_diff   = abs(VInf_EGA2_p - VInf_EGA2_m);
                solutionsStruc{solutionIndex}.DV               = solutionsStruc{solutionIndex}.VInf_VGA_diff + solutionsStruc{solutionIndex}.VInf_EGA2_diff;
                solutionsStruc{solutionIndex}.flybyRadius_VGA  = flybyRadius_VGA;
                solutionsStruc{solutionIndex}.flybyRadius_EGA1 = flybyRadius_EGA1;
                solutionsStruc{solutionIndex}.flybyRadius_EGA2 = flybyRadius_EGA2;
                solutionsStruc{solutionIndex}.VInf_VGA_m       = VInf_VGA_m;
                solutionsStruc{solutionIndex}.VInf_VGA_p       = VInf_VGA_p;
                solutionsStruc{solutionIndex}.VInf_EGA1_m      = VInf_EGA1_m;
                solutionsStruc{solutionIndex}.VInf_EGA1_p      = VInf_EGA1_p;
                solutionsStruc{solutionIndex}.VInf_EGA2_m      = VInf_EGA2_m;
                solutionsStruc{solutionIndex}.VInf_EGA2_p      = VInf_EGA2_p;
                solutionsStruc{solutionIndex}.VInf_SOI         = VInf_SOI_m;
                solutionsStruc{solutionIndex}.phi_EGA_rad      = myPhi_rad;
                
                solutionsStruc{solutionIndex}.leg1_X0 = [rEarth_Launch; vSc_Launch];
                solutionsStruc{solutionIndex}.leg2_X0 = [rVenus_VGA; vSc_VGA_p];
                solutionsStruc{solutionIndex}.leg3_X0 = [rEarth_EGA1; vSc_EGA1_p];
                solutionsStruc{solutionIndex}.leg4_X0 = [rEarth_EGA2; vSc_EGA2_p];
                
                solutionsStruc{solutionIndex}.BT_VGA  = BT_VGA;
                solutionsStruc{solutionIndex}.BR_VGA  = BR_VGA;
                solutionsStruc{solutionIndex}.BT_EGA1 = BT_EGA1;
                solutionsStruc{solutionIndex}.BR_EGA1 = BR_EGA1;
                solutionsStruc{solutionIndex}.BT_EGA2 = BT_EGA2;
                solutionsStruc{solutionIndex}.BR_EGA2 = BR_EGA2;
                
            end
        end
    end
end


% ---------------------
%%% Reorganize data
% ---------------------
%%% Deleting empty cells
if length(solutionsStruc(~cellfun('isempty',solutionsStruc))) == length(solutionsStruc)
    warning('Need to preallocate more space to solutionsStruc')
end
solutionsStruc = solutionsStruc(~cellfun('isempty',solutionsStruc));

%%% Reassigning results to get them out of cells
clear solutions
solutions = [solutionsStruc{:}];
solutions = orderfields(solutions);

%% =======================================================================
%%% Post Processing Data
% ========================================================================

% ---------------------
%%% Pick my solution
% ---------------------
solutions_DV = [solutions.DV];

mySolutionIndex = find(solutions_DV == min(solutions_DV));

mySolution = solutions(mySolutionIndex);

% ---------------------
%%% SOI estimated burn to achieve Mimas resonant orbit
% ---------------------
% %%% Burn at periapsis
% Tp_MimasResonantOrbit_3 = Mimas.Tp * 700;
% a_MimasResonantOrbit_3 = (((Tp_MimasResonantOrbit_3/(2*pi))^2)*Saturn.u)^(1/3);
% rp_MimasResonantOrbit_3 = Mimas.a;
% ra_MimasResonantOrbit_3 = 2*a_MimasResonantOrbit_3 - rp_MimasResonantOrbit_3;
% e_MimasResonantOrbit_3 = (ra_MimasResonantOrbit_3 - rp_MimasResonantOrbit_3)/(ra_MimasResonantOrbit_3 + rp_MimasResonantOrbit_3);
% [ vSc_MimasResonantOrbit_rp ] = visviva_v( rp_MimasResonantOrbit_3, a_MimasResonantOrbit_3, Saturn.u);
% 
% VMimas_SCI = sqrt(Mimas.a / Saturn.u); % km/s
% vSc_MimasCrossing_SCI = sqrt(mySolution.VInf_SOI^2 + 2*Saturn.u/Mimas.a);
% 
% mySolution.DV_SOI_MimasOrbit = vSc_MimasCrossing_SCI - vSc_MimasResonantOrbit_rp;
% 
% saturnSphereOfInfluence = Saturn.a * (Saturn.mass / Sun.mass)^(2/5);
% if ra_MimasResonantOrbit_3 < saturnSphereOfInfluence
%     fprintf('Within Saturn''s Sphere of Influence\n')
% else
%     warning('Not in Saturn''s Sphere of Influence\n')
% end
% 
% ra_MimasResonantOrbit_3
% mySolution.DV_SOI_MimasOrbit



% %%% Burn at periapsis
% % Tp_SOIOrbit_3 = Mimas.Tp * 700;
% rp_SaturnOrbit = 20000 + Saturn.R; % From Cassini
% ra_SaturnOrbit = Titan.a * 16;
% a_SaturnOrbit = (rp_SaturnOrbit + ra_SaturnOrbit) / 2;
% e_SaturnOrbit = (ra_SaturnOrbit - rp_SaturnOrbit)/(ra_SaturnOrbit + rp_SaturnOrbit);
% Tp_SaturnOrbit =2*pi*sqrt((a_SaturnOrbit^3)/Saturn.u);
% [ vSc_SaturnOrbit_rp ] = visviva_v( rp_SaturnOrbit, a_SaturnOrbit, Saturn.u);
% vSc_SOIOrbitCrossing_rp = sqrt(mySolution.VInf_SOI^2 + 2*Saturn.u/rp_SaturnOrbit);
% 
% mySolution.DV_SOI_SaturnOrbit = vSc_SOIOrbitCrossing_rp - vSc_SaturnOrbit_rp;
% mySolution.DV_SOI_SaturnOrbit
% 
% 
% mySolution.SaturnOrbit_ra = ra_SaturnOrbit;
% mySolution.SaturnOrbit_rp = rp_SaturnOrbit;
% mySolution.SaturnOrbit_a = a_SaturnOrbit;
% mySolution.SaturnOrbit_e = e_SaturnOrbit;
% mySolution.SaturnOrbit_Tp_sec = Tp_SaturnOrbit;
% mySolution.SaturnOrbit_Tp_day = Tp_SaturnOrbit/86400;
% mySolution.SaturnOrbit_Tp_yr = Tp_SaturnOrbit/86400/days_per_year;
% 
% mySolution = orderfields(mySolution)



%%% Burn at periapsis
Tp_SaturnOrbit = days_per_year*0.5*86400;
a_SaturnOrbit = (((Tp_SaturnOrbit/(2*pi))^2)*Saturn.u)^(1/3);
rp_SaturnOrbit = 20000 + Saturn.R; % From Cassini
ra_SaturnOrbit = 2*a_SaturnOrbit - rp_SaturnOrbit;
e_SaturnOrbit = (ra_SaturnOrbit - rp_SaturnOrbit)/(ra_SaturnOrbit + rp_SaturnOrbit);
% Tp_SaturnOrbit =2*pi*sqrt((a_SaturnOrbit^3)/Saturn.u);
[ vSc_SaturnOrbit_rp ] = visviva_v( rp_SaturnOrbit, a_SaturnOrbit, Saturn.u);
vSc_SOIOrbitCrossing_rp = sqrt(mySolution.VInf_SOI^2 + 2*Saturn.u/rp_SaturnOrbit);

mySolution.DV_SOI_SaturnOrbit = vSc_SOIOrbitCrossing_rp - vSc_SaturnOrbit_rp;
mySolution.DV_SOI_SaturnOrbit


mySolution.SaturnOrbit_ra = ra_SaturnOrbit;
mySolution.SaturnOrbit_rp = rp_SaturnOrbit;
mySolution.SaturnOrbit_a = a_SaturnOrbit;
mySolution.SaturnOrbit_e = e_SaturnOrbit;
mySolution.SaturnOrbit_Tp_sec = Tp_SaturnOrbit;
mySolution.SaturnOrbit_Tp_day = Tp_SaturnOrbit/86400;
mySolution.SaturnOrbit_Tp_yr = Tp_SaturnOrbit/86400/days_per_year;

mySolution = orderfields(mySolution)
% ---------------------
%%% Plot my solution
% ---------------------

%%% Choosing ode tolerance
tol = 1e-13;
%%% Setting integrator options
options = odeset('RelTol',tol,'AbsTol',tol);
n_dt = 1000;
time0_leg1 = linspace(mySolution.JD_Launch, mySolution.JD_VGA,n_dt) .*86400;
time0_leg2 = linspace(mySolution.JD_VGA,    mySolution.JD_EGA1,n_dt).*86400;
time0_leg3 = linspace(mySolution.JD_EGA1,   mySolution.JD_EGA2,n_dt).*86400;
time0_leg4 = linspace(mySolution.JD_EGA2,   mySolution.JD_SOI,n_dt) .*86400;

[~, X_sc_leg1] = ode113(@Int_2BI, time0_leg1, mySolution.leg1_X0, options, Sun.u);
[~, X_sc_leg2] = ode113(@Int_2BI, time0_leg2, mySolution.leg2_X0, options, Sun.u);
[~, X_sc_leg3] = ode113(@Int_2BI, time0_leg3, mySolution.leg3_X0, options, Sun.u);
[~, X_sc_leg4] = ode113(@Int_2BI, time0_leg4, mySolution.leg4_X0, options, Sun.u);

[~, X_Earth] = ode113(@Int_2BI, linspace(0,Earth.Tp,1000), [rEarth_Launch; vEarth_Launch], options, Sun.u);
[~, X_Venus] = ode113(@Int_2BI, linspace(0,Venus.Tp,1000), [rVenus_VGA; vVenus_VGA], options, Sun.u);
[~, X_Saturn] = ode113(@Int_2BI, linspace(0,Saturn.Tp,1000), [rSaturn_SOI; vSaturn_SOI], options, Sun.u);



nominalTrajColors = colorScale([colors.std.ltred; colors.std.ltblue],4);
lw = 4;
    
figure; hold all
p1 = plot3(X_sc_leg1(:,1),X_sc_leg1(:,2),X_sc_leg1(:,3),'color',nominalTrajColors(1,:),'linewidth',lw);
p2 = plot3(X_sc_leg2(:,1),X_sc_leg2(:,2),X_sc_leg2(:,3),'color',nominalTrajColors(2,:),'linewidth',lw);
p3 = plot3(X_sc_leg3(:,1),X_sc_leg3(:,2),X_sc_leg3(:,3),'color',nominalTrajColors(3,:),'linewidth',lw);
p4 = plot3(X_sc_leg4(:,1),X_sc_leg4(:,2),X_sc_leg4(:,3),'color',nominalTrajColors(4,:),'linewidth',lw);

[rSaturn_SOI, vSaturn_SOI] = JulianDate_2_HeliocentricState(mySolution.JD_EGA2, 'Saturn');

plot3(X_Earth(:,1),X_Earth(:,2),X_Earth(:,3),'k')
plot3(X_Venus(:,1),X_Venus(:,2),X_Venus(:,3),'k')
plot3(X_Saturn(:,1),X_Saturn(:,2),X_Saturn(:,3),'k')
plot3(rSaturn_SOI(1),rSaturn_SOI(2),rSaturn_SOI(3),'o','markerfacecolor',colors.std.black,'markersize',10)

legend([p1 p2 p3 p4],'Launch to VGA','VGA to EGA1','EGA1 to EGA2','EGA2 to SOI')
axis equal
PlotBoi3('X, km','Y, km', 'Z, km', 18, 'LaTex')
xlim([-11, 9].*1e8)
ylim([-2, 14].*1e8)

toc

