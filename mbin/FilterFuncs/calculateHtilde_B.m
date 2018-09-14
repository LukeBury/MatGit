%%% Inputs
% 1) Current Reference State [1x6]
% 2) Current Station State [1x6]
% 3) Current Station Number (1, 2, or 3)
% 4) Range Bias Associated w/ Current Station (km)
%%% Outputs
% 1) Htilde Matrix (based on Range and Range Rate equations)
function [ Htilde ] = calculateHtilde_B(refState, station, stationNum, pb)
Htilde = zeros(2,16);
rS = refState(1:3);
vS = refState(4:6);
rStat = station(1:3);
vStat = station(4:6);
p = norm(rS-rStat) + pb;

dPdr = (rS - rStat) / norm(rS - rStat);
dPdv = zeros(1,3);
dRPdr = (vS - vStat)/p - (dot((rS-rStat),(vS-vStat))*(rS-rStat))/(p^3);
dRPdv = (rS - rStat) / norm(rS - rStat);

Htilde(:,1:6) = [dPdr, dPdv;
         dRPdr, dRPdv];

x = rS(1); y = rS(2); z = rS(3);
dx = vS(1); dy = vS(2); dz = vS(3);
xS = rStat(1); yS = rStat(2); zS = rStat(3);
dxS = vStat(1); dyS = vStat(2); dzS = vStat(3);

if stationNum == 1
    Htilde(1,11:16) = [1, 0, 0, 0, 0, 0];
    Htilde(2,11:16) = [0, 0, 0,-((conj(x) - conj(xS))*(dx - dxS) + (conj(y) - conj(yS))*(dy - dyS) + (conj(z) - conj(zS))*(dz - dzS))/(pb + (abs(x - xS)^2 + abs(y - yS)^2 + abs(z - zS)^2)^(1/2))^2, 0, 0];
elseif stationNum == 2
    Htilde(1,11:16) = [0, 1, 0, 0, 0, 0];
    Htilde(2,11:16) = [0, 0, 0, 0,-((conj(x) - conj(xS))*(dx - dxS) + (conj(y) - conj(yS))*(dy - dyS) + (conj(z) - conj(zS))*(dz - dzS))/(pb + (abs(x - xS)^2 + abs(y - yS)^2 + abs(z - zS)^2)^(1/2))^2, 0];
elseif stationNum == 3
    Htilde(1,11:16) = [0, 0, 1, 0, 0, 0];
    Htilde(2,11:16) = [0, 0, 0, 0, 0, -((conj(x) - conj(xS))*(dx - dxS) + (conj(y) - conj(yS))*(dy - dyS) + (conj(z) - conj(zS))*(dz - dzS))/(pb + (abs(x - xS)^2 + abs(y - yS)^2 + abs(z - zS)^2)^(1/2))^2];
end
  
end
