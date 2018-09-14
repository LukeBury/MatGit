%%% Inputs
% 1) ECI S/C State [r (km); v (km/s)] [6x1]
% 2) Gravitational Parameter of Planet (km^3/s^2)
% 3) vInfinity Vector (km/s) [3x1]
%%% Outputs
% 1) BdT [3x1]
% 2) BdR [3x1]
% 3) ECI-to-BPlane DCM [3x3]
% 4) Linear Time of Flight
function [ BdT, BdR, DCM, LTOF ] = BPlaneCalculator(state, u, vInf)
if size(state) == [1,6]
    state = state';
end
if size(vInf) == [1,3]
    vInf = vInf';
end
r = state(1:3); % km
R = norm(r); % km
v = state(4:6); % km/s
V = norm(v); % km/s

% Semimajor axis
a = 1 / (2/R - (V^2)/u); % km
% Eccentricity Vector
eVec = ((V^2 - u/R) * r - dot(r,v) * v) / u;
% Eccentricity
e = norm(eVec);
% Turning Angle
% d = acos(1/e); % rads
% Semiminor axis
b = a * (e^2 - 1); % km

PHat = eVec/e;
WHat = cross(r,v)/norm(cross(r,v));
QHat = cross(WHat,PHat);

SHat = v/V;

kHat = [0, 0, 1]';
THat = cross(SHat,kHat)/norm(cross(SHat,kHat));

RHat = cross(SHat,THat);

BHat = cross(SHat,WHat);
BVec = b*BHat;

BdT = dot(BVec,THat);
BdR = dot(BVec,RHat);

ta = dot(r/R,PHat);
f = acosh(1 + (norm(vInf)^2)*a*(1-e^2)/(u*(1+e*cos(ta)))); % rads
% f = acosh(1 + (norm(vInf)^2)*a*(1-e^2)/(u*(1+e*ta))); % rads
LTOF = u*(sinh(f)-f)/(norm(vInf)^3); % sec
DCM = [SHat';THat';RHat'];

end

