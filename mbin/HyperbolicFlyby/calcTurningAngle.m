function [ d_rad ] = calcTurningAngle( rp, vInf, u )
%%% Calculates turning angle of hyperbolic flyby
%%% Inputs:
%       1) rp - radius (altitude) of flyby (km)
%       2) vInf - magnitude of v-infinity of flyby (km/s)
%       3) u - gravitational parameter of body (km^3/s^2)
%%% Outputs: 
%       1) d - turning angle of flyby (rad)
%=========================================================================
%%% Calculating turning angle
d_rad = pi - 2*acos(1/((rp*vInf^2)/u + 1)); % (rad)

end

