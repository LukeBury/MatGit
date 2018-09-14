%%% Inputs (full states, [1x6])
% 1) Satellite wrt Europa (Europa-centered inertial)
% 2) Europa wrt Jupiter (JCI)
% 3) Jupiter wrt Sun (HCI)
% 4) Earth wrt Sun (HCI)
% 5) Station wrt Earth (ECI)
%%% Outputs
% 1) Range Rate [1x1] % km/s
function [ rangeRate ] = calculateRangeRate_EurCI(sE,EJ,JSun,EaSun,stEa)
%%% Parsing States into relevant vectors
rsE = sE(1:3); % km
vsE = sE(4:6); % km/s
rEJ = EJ(1:3); % km
vEJ = EJ(4:6); % km/s
rJSun = JSun(1:3); % km
vJSun = JSun(4:6); % km/s
rEaSun = EaSun(1:3); % km
vEaSun = EaSun(4:6); % km/s
rstEa = stEa(1:3); % km
vstEa = stEa(4:6); % km/s

% ECI Range and Range-Rate from Europa-Cenetered-Inertial States
range = norm(rsE + rEJ + rJSun - rEaSun - rstEa); 
rangeRate = dot((rsE + rEJ + rJSun - rEaSun - rstEa),(vsE + vEJ + vJSun - vEaSun - vstEa))/range;
end

