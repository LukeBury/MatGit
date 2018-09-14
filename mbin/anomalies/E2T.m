% E = Eccentric anomaly (rad)
% e = eccentricity
function [v] = E2T(E,e)
v = 2 * atan(sqrt((1+e)/(1-e)) * tan(E/2));
while v < 0
    v = v + 2*pi;
end
while v > 2*pi
    v = v - 2*pi;
end
end