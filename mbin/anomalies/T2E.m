% v = True anomaly (rad)
% e = eccentricity
function [E] = T2E(v,e)
E = 2 * atan(tan(v/2)/sqrt((1+e)/(1-e)));
while E < 0
    E = E + 2*pi;
end
while E > 2*pi
    E = E - 2*pi;
end
end