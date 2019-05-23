%%% Inputs
% 1) Normalized time vector [1 x nt]
% 2) Normalized radius of secondary body
% 3-4) Gravitational parameter of primary/secondary body (km^3/s^2)
% 5-6) Latitude/Longitude of s/c starting position (deg)
% 7-8) Azimuth and elevation of s/c initial velocity (deg)
% 9) Normalized magnitude of s/c initial velocity
%%% Outputs
% 1) Normalized time vector [nt x 1]
% 2) Normalized s/c states [nt x 6]
function [Times_n, States_BCR_n] = normalCR3BP_getTraj(time, rad2_n, u, lat0, lon0, az, el, v_mag_n)
% ------------------------------------------------------------------------
%%% System Parameters
% ------------------------------------------------------------------------
%%% Rotating frame cooridinates
rB1_BCR_n = [-u, 0, 0];
rB2_BCR_n = [1-u, 0, 0];

% ------------------------------------------------------------------------
%%% Defining Particle State
% ------------------------------------------------------------------------
%%% Initial Particle Position
[rH0_ECEF_n] = latlon2surfECEF(lat0, lon0, rad2_n);
rH0_BCR_n = rH0_ECEF_n + rB2_BCR_n;

%%% Initial Partical Velocity
[vH0_BCR_n] = angleVelocity(v_mag_n, az, el, lat0, lon0);

% ------------------------------------------------------------------------
%%% Propagating State
% ------------------------------------------------------------------------
%%% Choosing ode45 tolerance
tol = 1e-13;

%%% Setting integrator options
options = odeset('Events',@event_Impact_CR3Bn,'RelTol',tol,'AbsTol',tol);

%%% Setting Initial State Vector (ECEF)
X0_n = [rH0_BCR_n, vH0_BCR_n]'; % km, km/s

%%% Propagating the State
prms.u = u;
prms.R2_n = rad2_n;
% [Times_n,States_BCR_n] = ode45(@normalCR3BP_Int,time,X0_n,options,u,rB1_BCR_n,rB2_BCR_n,rad2_n);
[Times_n,States_BCR_n] = ode45(@Int_CR3Bn,time,X0_n,options,prms);
end

