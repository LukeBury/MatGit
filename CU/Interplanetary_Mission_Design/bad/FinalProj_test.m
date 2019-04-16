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
JDs_Launch = linspace(JD_2025,JD_2025 + days_per_year, 75);
JDs_VGA    = linspace(2460751.5, 2460956.5, 75);
JDs_EGA1   = linspace(2460957.5, 2461456.5, 75);
JDs_EGA2   = JDs_EGA1 + resonance_period_days;
JDs_SOI    = linspace(2463060.737444906, JD_2041, 75);


%% =======================================================================
%%% Set requirements
% ========================================================================
% %%% C3 requirement on Earth Launch
% req_c3_Launch_ub = 25; % km^2/s^2
% 
% %%% Tolerance on abs(|vInf_m| - |vInf_p|) at VGA
% req_VInf_VGA_diff = 0.2; % km/s
% 
% %%% Flyby radius for Venus
% req_altitude_VGA = 150; % km
% req_flybyRadius_VGA_lb = Venus.R + req_altitude_VGA; % km
% 
% %%% Upper bound of VInf_SOI
% req_vInf_SOI_ub = 9; % km/s
% 
% %%% Upper bound of total travel time
% req_TOF_yrs_ub = 15;
% 
% %%% Tolerance on abs(|vInf_m| - |vInf_p|) at EGA2
% req_VInf_EGA2_diff = 0.2; % km/s
% 
% %%% Earth flyby altitude lower bound
% req_altitude_EGA = 150; % km
% req_flybyRadius_EGA_lb = Earth.R + req_altitude_EGA; % km

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
                end
                
                %%% Max rps of both flybys
                max_rp_EGA1_index = find(rps_EGA1 == max(rps_EGA1));
                max_rp_EGA2_index = find(rps_EGA2 == max(rps_EGA2));

                [ smallest_rpMax_EGA ] = smallerValue( rps_EGA1(max_rp_EGA1_index), rps_EGA2(max_rp_EGA2_index) ); % km
                
                %%% Check Earth flyby altitude requirement
                if smallest_rpMax_EGA < req_flybyRadius_EGA_lb
                    continue
                end
                
                %%% Check VInf_EGA2_diff requirement
                sum1 = sum1 + 1;
                
                if abs(VInf_EGA2_p - VInf_EGA2_m) > req_VInf_EGA2_diff
                    if abs(VInf_EGA2_p - VInf_EGA2_m) < minValue
                        minValue = abs(VInf_EGA2_p - VInf_EGA2_m);
                    end
                    continue
                end
                sum2 = sum2 + 1;
                
                %%% Store solution
                solutionIndex = solutionIndex + 1;
                
                solutionsStruc{solutionIndex}.JD_Launch = JD_Launch;
                solutionsStruc{solutionIndex}.JD_VGA    = JD_VGA;
                solutionsStruc{solutionIndex}.JD_EGA1   = JD_EGA1;
                solutionsStruc{solutionIndex}.JD_EGA2   = JD_EGA2;
                solutionsStruc{solutionIndex}.JD_SOI    = JD_SOI;
                
                solutionsStruc{solutionIndex}.c3_Launch           = c3_Launch;
                solutionsStruc{solutionIndex}.VInf_VGA_diff       = abs(VInf_VGA_p - VInf_VGA_m);
                solutionsStruc{solutionIndex}.VInf_EGA2_diff      = abs(VInf_EGA2_p - VInf_EGA2_m);
                solutionsStruc{solutionIndex}.DV                  = solutionsStruc{solutionIndex}.VInf_VGA_diff + solutionsStruc{solutionIndex}.VInf_EGA2_diff;
                solutionsStruc{solutionIndex}.flybyRadius_VGA     = flybyRadius_VGA;
                solutionsStruc{solutionIndex}.flybyRadius_EGA_min =  smallest_rpMax_EGA;
                solutionsStruc{solutionIndex}.VInf_SOI            = VInf_SOI_m;
                
                solutionsStruc{solutionIndex}.leg1_X0 = [rEarth_Launch; vSc_Launch];
                solutionsStruc{solutionIndex}.leg2_X0 = [rVenus_VGA; vSc_VGA_p];
                solutionsStruc{solutionIndex}.leg3_X0 = [rEarth_EGA1; vSc_EGA1_p];
                solutionsStruc{solutionIndex}.leg4_X0 = [rEarth_EGA2; vSc_EGA2_p];
                

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

% ---------------------
%%% Pick my solution
% ---------------------
solutions_DV = [solutions.DV];

mySolutionIndex = find(solutions_DV == min(solutions_DV));

mySolution = solutions(mySolutionIndex)

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

nominalTrajColors = colorScale([colors.std.ltred; colors.std.ltblue],4);
lw = 4;
    
figure; hold all
p1 = plot3(X_sc_leg1(:,1),X_sc_leg1(:,2),X_sc_leg1(:,3),'color',nominalTrajColors(1,:),'linewidth',lw);
p2 = plot3(X_sc_leg2(:,1),X_sc_leg2(:,2),X_sc_leg2(:,3),'color',nominalTrajColors(2,:),'linewidth',lw);
p3 = plot3(X_sc_leg3(:,1),X_sc_leg3(:,2),X_sc_leg3(:,3),'color',nominalTrajColors(3,:),'linewidth',lw);
p4 = plot3(X_sc_leg4(:,1),X_sc_leg4(:,2),X_sc_leg4(:,3),'color',nominalTrajColors(4,:),'linewidth',lw);

legend([p1 p2 p3 p4],'Launch to VGA','VGA to EGA1','EGA1 to EGA2','EGA2 to SOI')
axis equal
PlotBoi3('X, km','Y, km', 'Z, km', 18, 'LaTex')




