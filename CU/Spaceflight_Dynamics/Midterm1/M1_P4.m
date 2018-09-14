clear
clc
addpath('../../bin')
% ------------------------------------------------------------------------
%%% Givens
% ------------------------------------------------------------------------
rad2deg = 180/pi;
uE = 398600.4415; % km^3 / s^2
r_ECI = [1275.6; -6378.1; 0]; % km
v_ECI = [-0.158; 0.790; 9.09]; % km

% ------------------------------------------------------------------------
%%% 4a
% ------------------------------------------------------------------------
[a,e,i,raan,w,ta] = ECI2OE(r_ECI,v_ECI,uE);

%%% Printing results
fprintf('-------------------- 4a --------------------\n')
fprintf('a    = %2.4f km\n', a)
fprintf('e    = %2.4f\n', e)
fprintf('i    = %2.4f  deg\n', i*rad2deg)
fprintf('raan = %2.4f deg\n', raan*rad2deg)
fprintf('w    = %2.4f  deg\n', w*rad2deg)
fprintf('ta   = %2.4f deg\n\n', ta*rad2deg)

% ------------------------------------------------------------------------
%%% 4b
% ------------------------------------------------------------------------
time = 1000; % sec
finalTime = time; % sec
% Mean motion
n = sqrt(uE/(a^3)); % rad/s

%%% Determining initial mean anomaly
Ei = T2E(ta,e); % rads
Mi = E2M(Ei,e); % rads

pos = zeros(time,3);
vel = zeros(time,3);
for dt = 1:time
    %%% Calculate new true anomaly
    Mt = Mi + n*dt; % rads
    Et = M2E(Mt,e); % rads
    tat = E2T(Et,e); % rads
    
    %%% Calculate new ECI state
    [r_ECI_t, v_ECI_t] = OE2ECI(a, e, i, raan, w, tat, uE);
    
    %%% Storing new ECI state
    pos(dt,:) = r_ECI_t; % km
    vel(dt,:) = v_ECI_t; % km/s
    
    %%% Calculate latitude
    lat_t = pi/2 - acos(dot([0 0 1],r_ECI_t)/norm(r_ECI_t));
    lat_t*180/pi;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% Check if 60 deg North is crossed
    if lat_t > (60*pi/180)
        finalTime = dt;
        break
    end

%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Post-Exam Testing
%     [lat, lon] = ECEF2latlon(r_ECI_t);
%     if lat > (60*pi/180)
%         finalTime = dt;
%         break
%     end
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%%% Chopping off extra zeros from pos(t) and vel(t)
pos = pos(1:dt,:);
vel = vel(1:dt,:);

%%% Printing results
fprintf('-------------------- 4b --------------------\n')
fprintf('Time object takes to cross 60 North: %d seconds\n',finalTime)
fprintf('Or: %3.3f minutes\n\n',finalTime/60)

% ------------------------------------------------------------------------
%%% 4c
% ------------------------------------------------------------------------
%%% Printing results
fprintf('-------------------- 4c --------------------\n')
fprintf('Ending ECI position: [%f, %f, %f] km\n',pos(end,:))
fprintf('Ending ECI velocity: [%f, %f, %f] km/s\n\n',vel(end,:))

% ------------------------------------------------------------------------
%%% 4d
% ------------------------------------------------------------------------
time = 86400/2; % sec
finalTime = time; % sec
rE = 6378.1363; % km

pos2 = zeros(time,3);
range = zeros(time,1);
for dt = 1:time
    %%% Calculate new true anomaly
    Mt = Mi + n*dt; % rads
    Et = M2E(Mt,e); % rads
    tat = E2T(Et,e); % rads
    
    %%% Calculate new ECI state
    [r_ECI_t, v_ECI_t] = OE2ECI(a, e, i, raan, w, tat, uE);
    
    %%% Storing new ECI state
    pos2(dt,:) = r_ECI_t; % km
    range(dt,1) = norm(pos2(dt,:)); % km
    
    %%% Check if object impacted
    if norm(pos2(dt,:)) < rE
        finalTime2 = dt;
        break
    end
end

close all
plot(1:finalTime, range(:,1)-rE)
ylim([0 150])
PlotBoi('Time, sec','Altitude, km')

%%% Calculating lowest altitude
lowest = min(range(:,1)) - rE; % km

%%% Printing Results
fprintf('-------------------- 4d --------------------\n')
fprintf('Lowest altitude = %f km\n', lowest)



