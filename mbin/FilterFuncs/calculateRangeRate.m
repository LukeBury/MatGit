function [ rangeRate ] = calculateRangeRate(refState, station)
rS = refState(1:3);
vS = refState(4:6);
rStat = station(1:3);
vStat = station(4:6);
p = norm(rS-rStat);

rangeRate = dot((rS-rStat),(vS-vStat))/p;
end

