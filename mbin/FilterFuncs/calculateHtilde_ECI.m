%%% Inputs
% 1) Current Reference State [1x6]
% 2) Current Station State [1x6]
%%% Outputs
% 1) Htilde Matrix (based on Range and Range Rate equations)
function [ Htilde ] = calculateHtilde(refState, station)
rS = refState(1:3);
vS = refState(4:6);
rStat = station(1:3);
vStat = station(4:6);
p = norm(rS-rStat);

dPdr = (rS - rStat) / norm(rS - rStat);
dPdv = zeros(1,3);
dRPdr = (vS - vStat)/p - (dot((rS-rStat),(vS-vStat))*(rS-rStat))/(p^3);
dRPdv = (rS - rStat) / norm(rS - rStat);

Htilde = [dPdr, dPdv;
         dRPdr, dRPdv];
end
