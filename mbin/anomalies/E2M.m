% E = Eccentric Anomaly (rads)
% e = eccentricity
% M = Mean Anomaly (rads)
function [M] = E2M(E,e)
M = E - e*sin(E);
end