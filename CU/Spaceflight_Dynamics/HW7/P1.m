clear
clc
addpath('../../bin')
% ------------------------------------------------------------------------
%%% Givens
% ------------------------------------------------------------------------
r_ECIi = [5492; 3984.001; 2.955]; % km
v_ECIi = [-3.931; 5.498; 3.665]; % km
uE = 398600.4415; % km^3/s^2

% ------------------------------------------------------------------------
%%% Preparing for two-body propagation
% ------------------------------------------------------------------------
%%% Calculating Initial orbital elements
[a,e,i,raan,w,ta] = ECI2OE(r_ECIi,v_ECIi,uE);

% Initial eccentric anomaly
Ei = T2E(ta,e); % rads
% Initial mean anomaly
Mi = E2M(Ei,e); % rads
% Mean motion
n = sqrt(uE/(a^3)); % rad/s
% ------------------------------------------------------------------------
%%% Propagating two-body orbit
% ------------------------------------------------------------------------
ti = 0; % sec
tf = 1000000; % sec

for dt = ti:tf
    %%% Calculating new true anomaly
    Mt = Mi + n*dt; % rads
    Et = M2E(Mt,e); % rads
    tat = E2T(Et,e); % rads
    
    %%% Calculating new ECI state
    [r_ECIt, v_ECIt] = OE2ECI(a, e, i, raan, w, tat, uE);

    if dt == 100
        fprintf('At t = 100 sec:\n')
        r1_100 = r_ECIt;
        table(r1_100)
    end
end
fprintf('At t = 1,000,000 sec:\n')
r1_1000000 = r_ECIt;
table(r1_1000000)

clearvars -except r1_100 r1_1000000