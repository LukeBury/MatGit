function [BT, BR, B, theta_rad] = BPlaneTargeter_vInfs(vInf_minus, vInf_plus, mu, option_vInfMagnitude)
%%% Description
%       Calculates B-Plane targeting parameters from v-infinity-minus and
%       v-infinity-plus.
%       
%       *** NOTE: The equations require a |vInfinity|^2 term ... for this,
%       the magnitude of vInf_minus and vInf_plus are multiplied instead
% --------------------------------------------------------------
%%% Inputs
%       vInf_minus - Initial position vector (km/s) [3x1]
%       vInf_plus  - Target position vector (km/s) [3x1]
%       mu         - Gravitational parameter of central body (km^3/s^2)
%       option_vInfMagnitude - The equations require a |vInf|^2 term.
%           This can be dealt with in two ways: 
%           (1) |vInf|^2 = |vInf_minus|*|vInf_plus|
%           or
%           (2) |vInf|^2 = |vInf_minus|^2
% --------------------------------------------------------------
%%% Outputs
%       BT       - B-dot-T vector (km)
%       BR       - B-dot-R vector (km)
%       B         - B vector (km)
%       theta_rad - Theta vector (angle of B vector from B_t), (rad)
% --------------------------------------------------------------
%%% Author
%       Luke Bury, luke.bury@colorado.edu
% ===============================================================
%%% Magnitudes of vInfinities
VInf_minus = norm(vInf_minus); % km/s
VInf_plus  = norm(vInf_plus); % km/s

%%% Handling vInf^2 option
if option_vInfMagnitude == 1
    vInfSquared = VInf_minus*VInf_plus; % km/s
elseif option_vInfMagnitude == 2
    vInfSquared = VInf_minus*VInf_minus; % km/s
else
    warning('Incorrect Input')
    return
end

%%% Setting up vectors that define the B-Plane
k_hat = [0; 0; 1];
S_hat = vInf_minus ./ VInf_minus;
T_hat = cross(S_hat,k_hat)/norm(cross(S_hat,k_hat));
R_hat = cross(S_hat,T_hat)/norm(cross(S_hat,T_hat));
h_hat = cross(vInf_minus,vInf_plus)/norm(cross(vInf_minus,vInf_plus));
B_hat = cross(S_hat,h_hat)/norm(cross(S_hat,h_hat));

%%% Turning angle of flyby
turningAngle_rad = acos(dot(vInf_minus,vInf_plus)/(vInfSquared));

%%% Flyby radius
rFlyby = (mu / (vInfSquared))*((1)/(cos((pi-turningAngle_rad)/2)) - 1);

%%% 'b' vector
b = (mu / (vInfSquared)) * sqrt((1 + vInfSquared*rFlyby/mu)^2 - 1);

%%% 'B' vector
B = B_hat.*b;

%%% B-dot-T and B-dot-R
BT = dot(B,T_hat);
BR = dot(B,R_hat);

%%% Angle from B-dot-T to B
theta_rad = acos(dot(T_hat,B_hat));
end



