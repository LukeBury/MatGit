clear
clc
close all
mbinPath = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/mbin';
addpath(genpath(mbinPath))

% ========================================================================
%%% Info
% ========================================================================
% Using a patched-conics approach, computes the total dV necessary for a
% 2-burn transfer from a 300 km circular Earth orbit to a stationary
% landing on Europa and Enceladus. Burn 1 is at Earth departure, Burn 2 is
% at Europa/Enceladus landing. There are four "patched" frames: The Earth
% frame, the Sun frame, the Jupiter/Saturn frame, and the Europa/Enceladus
% frame

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
%%% Define bodies
% ------------------------------------
Sun       = bodies.sun;
Earth     = bodies.earth;
Jupiter   = bodies.jupiter;
Saturn    = bodies.saturn;
Europa    = bodies.europa;
Enceladus = bodies.enceladus;

% ========================================================================
%%% Calculations
% ========================================================================
% ------------------------------------
%%% Circular orbit velocities
% ------------------------------------
%%% Spacecraft about Earth in 300 km orbit
EarthAltitude = 300; % km
rSc_EarthCirc = Earth.R + EarthAltitude;
[ VSc_EarthCirc ]     = visviva_v( rSc_EarthCirc, rSc_EarthCirc, Earth.u); % km/s

%%% Earth about Sun
[ VEarth_Sun ]        = visviva_v( Earth.a, Earth.a, Sun.u);               % km/s

%%% Jupiter about Sun
[ VJupiter_Sun ]      = visviva_v( Jupiter.a, Jupiter.a, Sun.u);           % km/s

%%% Saturn about Sun
[ VSaturn_Sun ]       = visviva_v( Saturn.a, Saturn.a, Sun.u);             % km/s

%%% Europa about Jupiter
[ VEuropa_Jupiter ]   = visviva_v( Europa.a, Europa.a, Jupiter.u);         % km/s

%%% Enceladus about Saturn
[ VEnceladus_Saturn ] = visviva_v( Enceladus.a, Enceladus.a, Saturn.u);    % km/s

% ------------------------------------
%%% Transfer orbit parameters
% ------------------------------------
%%% Semi major axis of transfer orbits
a_trans_Jup = (Earth.a + Jupiter.a) / 2; % km
a_trans_Sat = (Earth.a + Saturn.a) / 2;  % km

%%% Periapsis velocity of transfer orbits
[ VSc_transP_Jup ] = visviva_v( Earth.a, a_trans_Jup, Sun.u); % km/s
[ VSc_transP_Sat ] = visviva_v( Earth.a, a_trans_Sat, Sun.u); % km/s

%%% Apoapsis velocity of transfer orbits
[ VSc_transA_Jup ] = visviva_v( Jupiter.a, a_trans_Jup, Sun.u); % km/s
[ VSc_transA_Sat ] = visviva_v( Saturn.a, a_trans_Sat, Sun.u);  % km/s

% ------------------------------------
%%% V-Infinities of transfer orbits
% ------------------------------------
%%% V-infinities at Earth departure
VSc_infEarth_Jup = VSc_transP_Jup - VEarth_Sun; % km/s
VSc_infEarth_Sat = VSc_transP_Sat - VEarth_Sun; % km/s

%%% V-infinities at planet arrival
VSc_infJupiter = VJupiter_Sun - VSc_transA_Jup; % km/s
VSc_infSaturn  = VSaturn_Sun - VSc_transA_Sat; % km/s

% ------------------------------------
%%% Post-dV1 velocities necessary to reach V-Infinity
% ------------------------------------
VSc_dV1_p_Jup = sqrt(VSc_infEarth_Jup^2 + 2*Earth.u/rSc_EarthCirc); % km/s
VSc_dV1_p_Sat = sqrt(VSc_infEarth_Sat^2 + 2*Earth.u/rSc_EarthCirc); % km/s

% ------------------------------------
%%% Planet-centric velocities at moon arrivals
% ------------------------------------
VSc_EuropaCrossing    = sqrt(VSc_infJupiter^2 + 2*Jupiter.u/Europa.a);  % km/s
VSc_EnceladusCrossing = sqrt(VSc_infSaturn^2 + 2*Saturn.u/Enceladus.a); % km/s

% ------------------------------------
%%% V-Infinities at moon arrivals
% ------------------------------------
VSc_infEuropa    = VSc_EuropaCrossing - VEuropa_Jupiter;      % km/s
VSc_infEnceladus = VSc_EnceladusCrossing - VEnceladus_Saturn; % km/s

% ------------------------------------
%%% Moon-centric velocity at landing site
% ------------------------------------
VSc_EuropaImpact    = sqrt(VSc_infEuropa^2 + 2*Europa.u/Europa.R);          % km/s
VSc_EnceladusImpact = sqrt(VSc_infEnceladus^2 + 2*Enceladus.u/Enceladus.R); % km/s

% ========================================================================
%%% Calculating Delta-Vs
% ========================================================================
% ------------------------------------
%%% dV1 (leaving Earth)
% ------------------------------------
VSc_dV1_m_Jup = VSc_EarthCirc; % km/s
VSc_dV1_m_Sat = VSc_EarthCirc; % km/s
dV1_Eur = VSc_dV1_p_Jup - VSc_dV1_m_Jup; % km/s
dV1_Enc = VSc_dV1_p_Sat - VSc_dV1_m_Sat; % km/s

% ------------------------------------
%%% dV2 (landing)
% ------------------------------------
%%% Europa
dV2_Europa_m = VSc_EuropaImpact; % km/s
dV2_Europa_p = 0;                % km/s

dV2_Eur = dV2_Europa_m - dV2_Europa_p; % km/s

%%% Enceladus
dV2_Enceladus_m = VSc_EnceladusImpact; % km/s
dV2_Enceladus_p = 0;                   % km/s

dV2_Enc = dV2_Enceladus_m - dV2_Enceladus_p; % km/s

% ------------------------------------
%%% Total simplistic mission Delta-V
% ------------------------------------
dV_EuropaLandingFromHohmann    = dV1_Eur + dV2_Eur % km/s
dV_EnceladusLandingFromHohmann = dV1_Enc + dV2_Enc % km/s





































