%%% Inputs
% 1) Current Reference State [1x6]
% 2) Current Station State [1x6]
% 3) Measurements [2x1]
%%% Outputs
% 1) Current Pre-Fit Residual
function [ yi ] = calculatePreFitResiduals(refState, station, measurements)
rS = refState(1:3);
vS = refState(4:6);
rStat = station(1:3);
vStat = station(4:6);
p = norm(rS-rStat);
[ dp ] = calculateRangeRate(refState, station);

yi = zeros(size(measurements));

yi(1) = measurements(1) - p;
yi(2) = measurements(2) - dp;
end

