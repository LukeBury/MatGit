% r = [x; y; z] (km)
% v = [xdot; ydot; zdot] (km/s)
% u = gravitational parameter of body being orbited (km^3 / s^2)
function [a,e,i,raan,w,ta] = RV2COE(r,v,u)
R = norm(r);
V = norm(v);
% For converting between degrees and radians
r2d = 180/pi;
d2r = pi/180;
% Semi Major Axis (km)
a = 1 / (2/R - (V^2)/u);
% Mean Motion (rad/s)
n = sqrt(u/(a^3));
% Eccentricity vector
e_vec = ((V^2 - u/R) * r - dot(r,v) * v) / u;
% Eccentricity
e = norm(e_vec);
% Radius of apoapsis and periapsis (km)
ra = a * (1 + e);
rp = a * (1 - e);
% Angular momentum vector (m^2/s)
h = cross(r,v);
% Defining axes
I = [1;0;0];
J = [0;1;0];
K = [0;0;1];
% Inclination (deg)
i = acos(dot(K,h)/norm(h));
i = i * r2d;
% raan (deg)
node = cross(K,h);
raan = acos(dot(I,node)/norm(node));
if node(2) < 0
    raan = 2*pi - raan;
end
raan = raan * r2d;
% Argument of periapse (deg)
w = acos(dot(node,e_vec)/(norm(node)*e));
if e_vec(3) < 0
    w = 2*pi - w;
end
w = w * r2d;
% True anomaly
ta = acos(dot(e_vec,r)/(e * R));
if dot(r,v) < 0
    ta = 2*pi - ta;
end
ta = ta * r2d;
end