clear
clc
% ------------------------------------------------------------------------
%%% Givens
% ------------------------------------------------------------------------
u_Earth = 398600; % km^3 / s^2
rE = 6378.1363; % km
a1 = rE + 6000; % km
phase = 30 * pi / 180; % rads

% ------------------------------------------------------------------------
%%% Calculations
% ------------------------------------------------------------------------
%%% Mean motion of target satellite
w_tgt = sqrt(u_Earth / a1 ^ 3); % rad/s

%%% Time between t0 and t1
dt = (2 * pi + phase) / w_tgt; % sec

%%% Phasing orbit semi major axis, rp, and ra
a_phase = (u_Earth * (dt / (2 * pi))^2 )^(1/3); % km
rp = a1; % km

%%% Calculating velocities
v_circ = sqrt(u_Earth / a1); % km/s
v_tran = sqrt(2 * u_Earth / rp - u_Earth/a_phase); % km/s

%%% Calculating dVs
dv1 = v_tran - v_circ; % km/s
% dV2 = dV1 because the satellite is returning to its original orbit
dv2 = dv1; % km/s
dv_Total = dv1 + dv2; % km/s

%%% Printing results
fprintf('dv_Total      = %2.4f km/s\n', dv_Total)
fprintf('Transfer time = %5.4f sec\n', dt)
fprintf('              = %1.4f hrs\n', dt/(60*60))

