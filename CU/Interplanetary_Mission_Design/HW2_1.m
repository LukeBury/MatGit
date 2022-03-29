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
%%% Problem 1
% ========================================================================
% ------------------------------------
%%% Setup
% ------------------------------------
%%% Launch date to julian date
year_0   = 2018;
month_0  = 5;
day_0    = 1;
hour_0   = 0;
minute_0 = 0;
second_0 = 0;
[ JD_Launch ] = calendar2julian(year_0, month_0, day_0, hour_0, minute_0, second_0);

% %%% Finding nominal Hohmann transfer time from Earth to Mars
% a_HohTrans = (bodies.earth.a + bodies.mars.a)/2;     % km
% T_HohTrans = pi * sqrt((a_HohTrans^3)/bodies.sun.u); % sec
% 
% %%% Setting transfer window vectors
% delta_transTime = 60*86400;
% T_trans_min = T_HohTrans - delta_transTime;
% T_trans_max = T_HohTrans + delta_transTime;

n_trajsPerType = 1500;
transferTimes_type1 = linspace(120*86400, 243*86400, n_trajsPerType);
transferTimes_type2 = linspace(243*86400, 375*86400, n_trajsPerType);

%%% Finding initial earth position
[L_rad, a, e, i_rad, Omega_rad, Pi_rad, w_rad, M_rad, ta_rad] = ...
    getPlanetElements_Meeus(JD_Launch, 'Earth', 'radians');
[r0_Earth, v0_Earth] = OE2ECI(a, e, i_rad, Omega_rad, w_rad, ta_rad, bodies.sun.u);

%%% Preallocating resultsfor type-1 and type-2
C3_departures_type1 = zeros(length(transferTimes_type1),1);
C3_departures_type2 = zeros(length(transferTimes_type2),1);
vInf_arrivals_type1 = zeros(length(transferTimes_type1),1);
vInf_arrivals_type2 = zeros(length(transferTimes_type2),1);

% ------------------------------------
%%% Type-1 Loop
% ------------------------------------
%%% Looping through type-1 possible trajectories
for kk = 1:length(transferTimes_type1)
    % ----------------------
    %%% Find TOF and final position of mars
    % ----------------------
    %%% time of flight
    TOF_sec = transferTimes_type1(kk);
    TOF_days = TOF_sec/86400;
    
    %%% Find Julida date of arrival
    JD_Arrival = JD_Launch + TOF_days; 
    
    %%% Find Mars position at arrival time
    [L_rad, a, e, i_rad, Omega_rad, Pi_rad, w_rad, M_rad, ta_rad] = ...
        getPlanetElements_Meeus(JD_Arrival, 'Mars', 'radians');
    [rf_Mars, vf_Mars] = OE2ECI(a, e, i_rad, Omega_rad, w_rad, ta_rad, bodies.sun.u);
    
    % ----------------------
    %%% Computing Lambert solution
    % ----------------------
    %%% Compute Lambert solution
    nOrbits = 0;
    [V1, V2,exitflag] = lambertSolver(r0_Earth, rf_Mars, TOF_sec, nOrbits, 0, bodies.sun.u);
    
    % ----------------------
    %%% Gathering Results
    % ----------------------
    %%% Compute C3
    C3_depart = (norm(V1 - v0_Earth))^2; % km^2/s^2
    
    %%% Compute arrival vInf
    vInf_arrival = norm(V2 - vf_Mars); % km/s
    
    %%% Store results
    C3_departures_type1(kk) = C3_depart;    % km^2/s^2
    vInf_arrivals_type1(kk) = vInf_arrival; % km/s
end

% ------------------------------------
%%% Type-2 Loop
% ------------------------------------
%%% Looping through type-1 possible trajectories
for kk = 1:length(transferTimes_type2)
    % ----------------------
    %%% Find TOF and final position of mars
    % ----------------------
    %%% time of flight
    TOF_sec = transferTimes_type2(kk);
    TOF_days = TOF_sec/86400;
    
    %%% Find Julida date of arrival
    JD_Arrival = JD_Launch + TOF_days; 
    
    %%% Find Mars position at arrival time
    [L_rad, a, e, i_rad, Omega_rad, Pi_rad, w_rad, M_rad, ta_rad] = ...
        getPlanetElements_Meeus(JD_Arrival, 'Mars', 'radians');
    [rf_Mars, vf_Mars] = OE2ECI(a, e, i_rad, Omega_rad, w_rad, ta_rad, bodies.sun.u);
    
    % ----------------------
    %%% Computing Lambert solution
    % ----------------------
    %%% Compute Lambert solution
    nOrbits = 0;
    [V1, V2,exitflag] = lambertSolver(r0_Earth, rf_Mars, TOF_sec, nOrbits, 0, bodies.sun.u);
   
    % ----------------------
    %%% Gathering Results
    % ----------------------
    %%% Compute C3
    C3_depart = (norm(V1 - v0_Earth))^2; % km^2/s^2
    
    %%% Compute arrival vInf
    vInf_arrival = norm(V2 - vf_Mars); % km/s
    
    %%% Store results
    C3_departures_type2(kk) = C3_depart;    % km^2/s^2
    vInf_arrivals_type2(kk) = vInf_arrival; % km/s
end

% ------------------------------------
%%% Finding minimum values and associated transfer times
% ------------------------------------
%%% Type 1
C3_min_type1   = min(C3_departures_type1)
vInf_min_type1 = min(vInf_arrivals_type1)

transTime_min_C3_type1   = transferTimes_type1(find(C3_departures_type1 == C3_min_type1));
transTime_min_vInf_type1 = transferTimes_type1(find(vInf_arrivals_type1 == vInf_min_type1));

%%% Type 2
C3_min_type2   = min(C3_departures_type2)
vInf_min_type2 = min(vInf_arrivals_type2)

transTime_min_C3_type2   = transferTimes_type2(find(C3_departures_type2 == C3_min_type2));
transTime_min_vInf_type2 = transferTimes_type2(find(vInf_arrivals_type2 == vInf_min_type2));

% ------------------------------------
%%% Finding associated Mars arrival dates
% ------------------------------------
%%% Type 1
JD_Arrival_min_C3_type1   = JD_Launch + transTime_min_C3_type1/86400;
JD_Arrival_min_vInf_type1 = JD_Launch + transTime_min_vInf_type1/86400;

[Year, Mon, Day, hr, min, s] = julian2calendar(JD_Arrival_min_C3_type1);
date_C3min_type1 = [Year, Mon, Day, hr, min, s]'

[Year, Mon, Day, hr, min, s] = julian2calendar(JD_Arrival_min_vInf_type1);
date_vInfmin_type1 = [Year, Mon, Day, hr, min, s]'

%%% Type 2
JD_Arrival_min_C3_type2   = JD_Launch + transTime_min_C3_type2/86400;
JD_Arrival_min_vInf_type2 = JD_Launch + transTime_min_vInf_type2/86400;

[Year, Mon, Day, hr, min, s] = julian2calendar(JD_Arrival_min_C3_type2);
date_C3min_type2 = [Year, Mon, Day, hr, min, s]'

[Year, Mon, Day, hr, min, s] = julian2calendar(JD_Arrival_min_vInf_type2);
date_vInfmin_type2 = [Year, Mon, Day, hr, min, s]'

% ------------------------------------
%%% Plotting
% ------------------------------------
%%% Type 1
figure('position',[440 444 739 354])
subplot(1,2,1); hold all
plot(transferTimes_type1/86400, C3_departures_type1,'linewidth',2,'color',colors.std.purp)
plot(transTime_min_C3_type1/86400, C3_min_type1,'o',...
    'markeredgecolor',colors.std.black,'markerfacecolor',colors.std.ltred,'markersize',10)
PlotBoi2('Transfer Time, $days$','$C_3$, $km^2/s^2$',18,'LaTex')
ylim([0 100])

subplot(1,2,2) ; hold all
plot(transferTimes_type1/86400, vInf_arrivals_type1,'linewidth',2,'color',colors.std.purp)
plot(transTime_min_vInf_type1/86400, vInf_min_type1,'o',...
    'markeredgecolor',colors.std.black,'markerfacecolor',colors.std.ltred,'markersize',10)
PlotBoi2('Transfer Time, $days$','$v_\infty$, $km/s$',18,'LaTex')

text(83, 31.5,'Type 1','FontSize',13,'FontName','Times New Roman')

%%% Type 2
figure('position',[440 444 739 354])
subplot(1,2,1); hold all
plot(transferTimes_type2/86400, C3_departures_type2,'linewidth',2,'color',colors.std.purp)
plot(transTime_min_C3_type2/86400, C3_min_type2,'o',...
    'markeredgecolor',colors.std.black,'markerfacecolor',colors.std.ltred,'markersize',10)
PlotBoi2('Transfer Time, $days$','$C_3$, $km^2/s^2$',18,'LaTex')
ylim([0 100])

subplot(1,2,2) ; hold all
plot(transferTimes_type2/86400, vInf_arrivals_type2,'linewidth',2,'color',colors.std.purp)
plot(transTime_min_vInf_type2/86400, vInf_min_type2,'o',...
    'markeredgecolor',colors.std.black,'markerfacecolor',colors.std.ltred,'markersize',10)
PlotBoi2('Transfer Time, $days$','$v_\infty$, $km/s$',18,'LaTex')
ylim([0 25])

text(203, 26.5,'Type 2','FontSize',13,'FontName','Times New Roman');








% %%% 
% nOrbits = 0;
% [V1, V2,exitflag] = lambertSolver(r1vec, r2vec, tf, nOrbits, 0, bodies.sun.u)



% [L_deg, a, e, i_deg, Omega_deg, Pi_deg, w_deg, M_deg, ta_deg] = ...
%     getPlanetElements_Meeus(JDE, bodyName, 'degrees');
% 
% yr = 2000;
% mo = 1;
% day = 1;
% hr = 12;
% min = 0;
% sec = 0;
% 
% [ JD ] = calendar2julian(yr, mo, day, hr, min, sec)
% 
% 
% % [yr,mo,d,hr,min,s,dayweek,dategreg] = julian2greg(JD)
% 
% % [Year, Mon, Day, hr, min, s] = julian2greg2(JD)
% [Year, Mon, Day, hr, min, s] = julian2greg2(2450383.09722222)
