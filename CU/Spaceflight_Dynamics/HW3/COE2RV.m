% Input angles in DEGREES
function [r, v] = COE2RV(a, e, i, raan, w, ta, u)
d2r = pi/180;
% Semiparameter (km)
p = a * (1 - e^2);

% Angular momentum (km^2/s)
h = sqrt(p*u);
% r_pqw (km)
r1_pqw = p * cos(ta) / (1 + e * cos(ta));
r2_pqw = p * sin(ta) / (1 + e * cos(ta));
r3_pqw = 0;
r_pqw = [r1_pqw; r2_pqw; r3_pqw]
% v_pqw (km/s)
v1_pqw = -sqrt(u/p) * sin(ta);
v2_pqw = sqrt(u/p) * (e + cos(ta));
v3_pqw = 0;
v_pqw = [v1_pqw; v2_pqw; v3_pqw];
% R1(-i)
R1_i = [1 0 0; 0 cos(-i) sin(-i); 0 -sin(-i) cos(-i)];
% R3(-raan)
R3_raan = [cos(-raan) sin(-raan) 0; -sin(-raan) cos(-raan) 0; 0 0 1];
% R3(-w)
R3_w = [cos(-w) sin(-w) 0; -sin(-w) cos(-w) 0; 0 0 1];
% Rotating for final r & v
r = R3_raan * R1_i * R3_w * r_pqw;
v = R3_raan * R1_i * R3_w * v_pqw;
end
