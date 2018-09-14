clear
clc
% ------------------------------------------------------------------------
%%% Givens
% ------------------------------------------------------------------------
u_Sun = 1.32712428E11; % km^3 / s^2

a_Venus   = 108208601; % km
a_Earth   = 149598023; % km
a_Mars    = 227939186; % km
a_Jupiter = 778298361; % km

%%% Seconds-to-years conversion
sec2yrs = 1/(365.25*24*60*60); % years/seconds

%% =======================================================================
fprintf('===================== PROBLEM 1 =====================\n\n')

% ------------------------------------------------------------------------
%%% Orbital Velocities
% ------------------------------------------------------------------------
v_Venus   = sqrt(u_Sun/a_Venus); % km/s
v_Earth   = sqrt(u_Sun/a_Earth); % km/s
v_Mars    = sqrt(u_Sun/a_Mars); % km/s
v_Jupiter = sqrt(u_Sun/a_Jupiter); % km/s

% ------------------------------------------------------------------------
%%% Venus Transfer
% ------------------------------------------------------------------------
%%% Calculating semi-major axis, v1, and v2 of transfer orbit
at_Venus = (a_Venus + a_Earth) / 2; % km
vt_1 = sqrt(2*u_Sun/a_Earth - u_Sun/at_Venus); % km/s
vt_2 = sqrt(2*u_Sun/a_Venus - u_Sun/at_Venus); % km/s

%%% Calculating dVs
dv1_Venus = v_Earth - vt_1; % km/s
dv2_Venus = vt_2 - v_Venus; % km/s
dvTotal_Venus = dv1_Venus + dv2_Venus; % km/s

%%% Calculating time period
T_Venus = pi * sqrt(at_Venus^3 / u_Sun); % seconds

fprintf('a_Venus         = %9d km\n', a_Venus)
fprintf('dv1_Venus       = %2.4f km/s\n', dv1_Venus)
fprintf('dv2_Venus       = %2.4f km/s\n', dv2_Venus)
fprintf('dvTotal_Venus   = %2.4f km/s\n', dvTotal_Venus)
fprintf('T_Venus         = %1.5f years\n\n', T_Venus * sec2yrs)
% ------------------------------------------------------------------------
%%% Mars Transfer
% ------------------------------------------------------------------------
%%% Calculating semi-major axis, v1, and v2 of transfer orbit
at_Mars = (a_Mars + a_Earth) / 2; % km
vt_1 = sqrt(2*u_Sun/a_Earth - u_Sun/at_Mars); % km/s
vt_2 = sqrt(2*u_Sun/a_Mars - u_Sun/at_Mars); % km/s

%%% Calculating dVs
dv1_Mars = vt_1 - v_Earth; % km/s
dv2_Mars = v_Mars - vt_2; % km/s
dvTotal_Mars = dv1_Mars + dv2_Mars; % km/s

%%% Calculating time period
T_Mars = pi * sqrt(at_Mars^3 / u_Sun); % seconds

fprintf('a_Mars          = %9d km\n', a_Mars)
fprintf('dv1_Mars        = %2.4f km/s\n', dv1_Mars)
fprintf('dv2_Mars        = %2.4f km/s\n', dv2_Mars)
fprintf('dvTotal_Mars    = %2.4f km/s\n', dvTotal_Mars)
fprintf('T_Mars          = %1.5f years\n\n', T_Mars * sec2yrs)
% ------------------------------------------------------------------------
%%% Jupiter Transfer
% ------------------------------------------------------------------------
%%% Calculating semi-major axis, v1, and v2 of transfer orbit
at_Jupiter = (a_Jupiter + a_Earth) / 2; % km
vt_1 = sqrt(2*u_Sun/a_Earth - u_Sun/at_Jupiter); % km/s
vt_2 = sqrt(2*u_Sun/a_Jupiter - u_Sun/at_Jupiter); % km/s

%%% Calculating dVs
dv1_Jupiter = vt_1 - v_Earth; % km/s
dv2_Jupiter = v_Jupiter - vt_2; % km/s
dvTotal_Jupiter = dv1_Jupiter + dv2_Jupiter; % km/s

%%% Calculating time period
T_Jupiter = pi * sqrt(at_Jupiter^3 / u_Sun); % seconds

fprintf('a_Jupiter       = %9d km\n', a_Jupiter)
fprintf('dv1_Jupiter     = %2.4f km/s\n', dv1_Jupiter)
fprintf('dv2_Jupiter     = %2.4f km/s\n', dv2_Jupiter)
fprintf('dvTotal_Jupiter = %2.4f km/s\n', dvTotal_Jupiter)
fprintf('T_Jupiter       = %1.5f years\n', T_Jupiter * sec2yrs)

%% =======================================================================
fprintf('\n===================== PROBLEM 2 =====================\n\n')

%%% Calculating mean motions
n_Mars  = sqrt(u_Sun / a_Mars^3); % rad/s
n_Earth = sqrt(u_Sun / a_Earth^3); % rad/s

%%% Calculating phase angle and synodic period
phaseAngle = pi - T_Mars * n_Mars; % rads
synodicPeriod = 2 * pi / (n_Earth - n_Mars); % sec

fprintf('Phasing Angle  = %2.4f rads\n', phaseAngle)
fprintf('Synodic Period = %1.5f years\n', synodicPeriod * sec2yrs)

