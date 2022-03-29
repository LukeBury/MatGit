function [ v ] = visviva_v( r, a, u)
%%% Calculates velocity based on a position from the vis viva equation
%%% Inputs:
% 1) r - [nx1] distance to primary % km
% 2) a - [nx1] semimajor axis % km
% 2) u - gravitational parameter % km^3/s^2
%%% Outputs:
% 1) v - [nx1] velocity at r in orbit % km/s
% =========================================================================
if size(r,1) ~= size(a,1)
    warning('Dimension error')
    return
end
n = size(r,1);
v = zeros(n,1);
for kk = 1:n
    v(kk) = sqrt(u*(2/r(kk) - 1/a(kk)));
end
end

