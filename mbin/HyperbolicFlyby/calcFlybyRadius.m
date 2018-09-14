function [ rp ] = calcFlybyRadius( d, vInf, u )
%%% Calculates flyby radius (altitude) of hyperbolic flyby
%%% Inputs:
%       1) d - turning angle of flyby (rad)
%       2) vInf - v-infinity of flyby (km/s)
%       3) u - gravitational parameter of body (km^3/s^2)
%%% Outputs: 
%       1) rp - flyby radius of hyperbolic flyby
%=========================================================================
%%% Calculating flyby radius
rp = (u/(vInf^2))*(1/(cos((pi-d)/2)) - 1); % km

end

