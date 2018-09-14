function [ vecOut ] = azEl_2_Vec( refVec, az, el )
%%% Calculates vector from azimuth, elevation, magnitude
%%% Inputs:
%       1) refVec - reference vector from which azimuth and elevation are
%       measured
%       2) azimuth - planar angle from +y (deg)
%       3) elevation - elevation angle from x-y plane (deg)
%%% Outputs: 
%       1) vecOut - vector in desired direction with desired magnitude
%=========================================================================
v1 = R1(refVec, el*pi/180);

vecOut = R3(v1, -az*pi/180);
end

