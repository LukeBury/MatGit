clear
clc
addpath('../../bin')
% ------------------------------------------------------------------------
%%% Givens
% ------------------------------------------------------------------------
rE = 6378.1363; % km
alt = 500; % km
uE = 398600.4415; % km^3 / s^2
deg2rad = pi/180;
rad2deg = 180/pi;

% ------------------------------------------------------------------------
%%% Determining fuel required for 10 deg inclination change
% ------------------------------------------------------------------------
%%% Creating circular orbit
a = rE + alt; % km
e = 0;
i = 0; 
raan = 0;
w = 0; 
ta = 0; 

%%% Determining velocity of non-inclined circular orbit
[r1, v1] = OE2ECI(a, e, i, raan, w, ta, uE)

%%% Determining velocity of inclined circular orbit
i2 = i + 10*deg2rad;
[r2, v2] = OE2ECI(a, e, i2, raan, w, ta, uE)

%%% Calculating required dV (aka, fuel amount of satellite)
dV = norm(v1-v2) % km/s

% ------------------------------------------------------------------------
%%% Comparing dV to apoapsis
% ------------------------------------------------------------------------
ras = zeros(length(0:.001:dV),1);
dVs = zeros(length(0:.001:dV),1);
index = 1;
v1 = norm(v1);
% for dv = 0:.001:dV%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% bad
for dv = 0:.0001:dV
    %%% Determining new orbit apoapsis
    rp = a; % km
    % New velocity
    v = v1 + dv; % km/s
    % New semi-major axis
    at = 1/(2/rp - (v^2)/uE); % km
    % New eccentricity
    et = 1 - rp/at;
%     ra = at*(1+e); % km%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% bad
    ra = at*(1+et); % km
    
    %%% Storing info
    ras(index,1) = ra; % km
    dVs(index,1) = dv; % km/s
    index = index + 1;
    
    %%% If using half fuel, store distance
%     if index == round(length(0:.001:dV)/2)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if index == round(length(0:.0001:dV)/2)  
        halfFuelTime = pi*sqrt((at^3)/uE); % sec
        halfFuelDist = ra; % km
    end

end

%%% Calculating furthest point and time to get there
far = max(ras(:,1)); % km
time = pi*sqrt((at^3)/uE); % sec

close all
plot(dVs,ras,'linewidth',2)
PlotBoi2('\DeltaV, km/s', 'Apoapsis, km',18)
grid on

fprintf('The furthest the satellite can reach is: %f km\n',far)
fprintf('The time it takes to get there is: %f seconds\n',time)
fprintf('Or: %f hours\n\n', time/3600)

fprintf('Using half its fuel, the satellite would make it to: %f km\n',...
    halfFuelDist)
fprintf('in %f hours', halfFuelTime/3600)






















