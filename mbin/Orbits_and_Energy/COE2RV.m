function [r, v] = COE2RV(a, e, i, raan, w, ta, mu)
%%% Description
%       Transforms classical orbital elements to an inertial state
%       
% --------------------------------------------------------------
%%% Inputs
%       a    - [1x1] semi-major axis (km) 
%       e    - [1x1] eccentricity 
%       i    - [1x1] inclination (deg)
%       raan - [1x1] right ascension of the ascending node (deg)
%       w    - [1x1] argument of periapsis (deg)
%       ta   - [1x1] true anomaly (deg)
%       mu   - [1x1] gravitational parameter of primary body (km^3/s^2)
% --------------------------------------------------------------
%%% Outputs
%       r - [3x1] inertial position vector (km)
%       v - [3x1] inertial velocity vector (km/s)
% --------------------------------------------------------------
%%% Author
%       Luke Bury, luke.bury@colorado.edu
% ===============================================================
%%% Semiparameter
p = a * (1 - e^2); % (km)

% %%% Angular momentum
% h = sqrt(p*mu); % (km^2/s)

%%% Position in PQW frame 
r1_PQW = p * cosd(ta) / (1 + e * cosd(ta)); % (km)
r2_PQW = p * sind(ta) / (1 + e * cosd(ta)); % (km)
r3_PQW = 0; % (km)
r_PQW = [r1_PQW; r2_PQW; r3_PQW]; % (km)

%%% Velocity in PQW frame
v1_PQW = -sqrt(mu/p) * sind(ta); % (km/s)
v2_PQW = sqrt(mu/p) * (e + cosd(ta)); % (km/s)
v3_PQW = 0; % (km/s)
v_PQW = [v1_PQW; v2_PQW; v3_PQW]; % (km/s)

%%% Creating rotation matrices to rotate state from PQW to inertial
% R1(-i)
R1_i = [1 0 0; 0 cosd(-i) sind(-i); 0 -sind(-i) cosd(-i)];
% R3(-raan)
R3_raan = [cosd(-raan) sind(-raan) 0; -sind(-raan) cosd(-raan) 0; 0 0 1];
% R3(-w)
R3_w = [cosd(-w) sind(-w) 0; -sind(-w) cosd(-w) 0; 0 0 1];

%%% Rotating state from PQW to inertial
r = R3_raan * R1_i * R3_w * r_PQW;
v = R3_raan * R1_i * R3_w * v_PQW;
end
