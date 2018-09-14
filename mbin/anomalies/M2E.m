% M = Mean anomaly (rad)
% e = eccentricity
function [E] = M2E(M,e)
tol = 1e-10;
E0 = M;
E1 = E0 + (M - E0 + e * sin(E0))/(1 - e * cos(E0));
while abs(E1-E0) > tol
    E0 = E1;
    E1 = E0 + (M - E0 + e * sin(E0))/(1 - e * cos(E0));
end
E = E1;
while E < 0
    E = E + 2*pi;
end
while E > 2*pi
    E = E - 2*pi;
end
end
