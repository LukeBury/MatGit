function [ a ] = visviva_a( r, v, u)
%%% Calculates semimajor axis based on a position and velocity from the vis viva equation
%%% Inputs:
% 1) r - [nx1] distance to primary % km
% 2) v - [nx1] velocity at instance % km/s
% 2) u - gravitational parameter % km^3/s^2
%%% Outputs:
% 1) a - [nx1] semimajor axis of orbit % km/s
% =========================================================================
if size(r,1) ~= size(v,1)
    warning('Dimension error')
    return
end
n = size(r,1);
a = zeros(n,1);
for kk = 1:n
    a(kk) = 1/(2/r(kk) - v(kk)^2/u);
end
end

