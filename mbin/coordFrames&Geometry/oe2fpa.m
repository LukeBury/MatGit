function [ fpa ] = oe2fpa( e, v)
%%% Inputs:
% 1) e - [nx1] eccentricty
% 2) v - [nx1] true anomaly
%%% Outputs:
% 1) fpa - [nx1] (rad) flight path angle
% =========================================================================
n = size(e,1);
fpa = zeros(n,1);
for kk = 1:n
    fpa(kk) = atan2(e(kk)*sin(v(kk)),1+e(kk)*cos(v(kk))); % rad
end
end

