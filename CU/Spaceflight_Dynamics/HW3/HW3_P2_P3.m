clear
clc
addpath('../../bin')
u = 398600.4418; % km^3 / s^2
[a, e, i, raan, w, M, n, year, day] = tle2RV('ISS_TLE.txt');
E = M2E(M,e); % rads
ta = E2T(E,e); % rads

%%% #2, ECI r & v
[r, v] = COE2RV(a, e, i, raan, w, ta, u) % km & km/s

%%% #3
% tsp = 1 hr + epoch of TLE
time_since_periapse = 60 * 60 + (E-e*sin(E))/n; % seconds
M2 = n * time_since_periapse; % rads
E2 = M2E(M2,e); % rads
ta2 = E2T(E2,e); % rads

[r2, v2] = COE2RV(a, e, i, raan, w, ta2, u) % km & km/s

% Marielle's
%r_pqw = [-1.839749708612435   6.514869048605542   0]e3
%r2 = [-6.119752050727281  -0.407224837445325  -2.865484336059941]e3