%%% Inputs
% 1) Satellite wrt Europa (Europa-centered inertial)
% 2) Europa wrt Jupiter (JCI)
% 3) Jupiter wrt Sun (HCI)
% 4) Earth wrt Sun (HCI)
% 5) Station wrt Earth (ECI)
% 6) Measurements [2x1]
% 7) Range Rate [1x1]
%%% Outputs
% 1) Current Pre-Fit Residual
function [ yi ] = calculatePreFitResiduals_EurCI(sE, EJ, JSun, EaSun, stEa, measurements)
%%% Parsing States into relevant vectors
rsE = sE(1:3); % km
rEJ = EJ(1:3); % km
rJSun = JSun(1:3); % km
rEaSun = EaSun(1:3); % kmz
rstEa = stEa(1:3); % km

% ECI Range and Range-Rate from Europa-Cenetered-Inertial States
p = norm(rsE + rEJ + rJSun - rEaSun - rstEa); % km
dp = calculateRangeRate_EurCI(sE, EJ, JSun, EaSun, stEa);
yi = zeros(size(measurements));

yi(1) = measurements(1) - p;
yi(2) = measurements(2) - dp;
end

