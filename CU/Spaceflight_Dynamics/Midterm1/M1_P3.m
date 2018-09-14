clear
clc
close all
addpath('../../bin')
rad2deg = 180/pi;
deg2rad = pi/180;
% ------------------------------------------------------------------------
%%% Givens
% ------------------------------------------------------------------------
alt = 380; % km
rE = 6378.1363; % km
a = alt + rE; % km
uE = 398600.4415; % km^3 / s^2

%%% Initial satellite states
xA0 = 0; yA0 = 0; zA0 = 0; % m
xB0 = 0; yB0 = 0; zB0 = 0; % m

xAdot0 = 0; yAdot0 = 0; zAdot0 = 0.2; % m/s
xBdot0 = 0; yBdot0 = 0.4; zBdot0 = 0; % m/s

% ------------------------------------------------------------------------
%%% 3a
% ------------------------------------------------------------------------
%%% Determine wt
dt = 8*60; % sec
w = sqrt(uE / (a^3)); % rad/s
wt = w*dt; % rads

%%% Determining position of satellite A
xA_t1 = (xAdot0/w)*sin(wt) - (3*xA0 + 2*yAdot0/w)*cos(wt) + (4*xA0 + 2*yAdot0/w); % m
yA_t1 = (6*xA0 + 4*yAdot0/w)*sin(wt) + 2*xAdot0/w*cos(wt) - (6*w*xA0 + 3*yAdot0)*dt +...
    (yA0 - 2*xAdot0/w); % m
zA_t1 = zA0*cos(wt) + zAdot0/w*sin(wt); % m

%%% Determining position of satellite B
xB_t1 = (xBdot0/w)*sin(wt) - (3*xB0 + 2*yBdot0/w)*cos(wt) + (4*xB0 + 2*yBdot0/w); % m
yB_t1 = (6*xB0 + 4*yBdot0/w)*sin(wt) + 2*xBdot0/w*cos(wt) - (6*w*xB0 + 3*yBdot0)*dt +...
    (yB0 - 2*xBdot0/w); % m
zB_t1 = zB0*cos(wt) + zBdot0/w*sin(wt); % m

fprintf('------------------- 3a -------------------\n')
fprintf('After 8 minutes, \n')
fprintf('Sat A: RSW = [%f, %f, %f] meters\n', xA_t1, yA_t1, zA_t1)
fprintf('Sat b: RSW = [%f, %f, %f] meters\n\n', xB_t1, yB_t1, zB_t1)

% ------------------------------------------------------------------------
%%% 3b (trick: make B the origin?)
% ------------------------------------------------------------------------
%%% Correct 8 minute position so that satellite B is now the origin
xA_t1 = xA_t1 - xB_t1; % m
yA_t1 = yA_t1 - yB_t1; % m
zA_t1 = zA_t1 - zB_t1; % m

xB_t1 = 0; % m
yB_t1 = 0; % m
zB_t1 = 0; % m

%%% Correct 8 minute velocity (to get vA relative to vB)
xAdot1 = xAdot0*cos(wt) + (3*w*xA0 + 2*yAdot0)*sin(wt); % m/s
yAdot1 = (6*w*xA0 + 4*yAdot0)*cos(wt) - 2*xAdot0*sin(wt) - (6*w*xA0 + 3*yAdot0); % m/s
zAdot1 = -zA0*w*sin(wt) + zAdot0*cos(wt); % m/s

xBdot1 = xBdot0*cos(wt) + (3*w*xB0 + 2*yBdot0)*sin(wt); % m/s
yBdot1 = (6*w*xB0 + 4*yBdot0)*cos(wt) - 2*xBdot0*sin(wt) - (6*w*xB0 + 3*yBdot0); % m/s
zBdot1 = -zB0*w*sin(wt) + zBdot0*cos(wt); % m/s

xAdot1 = xAdot1 - xBdot1; % m/s
yAdot1 = yAdot1 - yBdot1; % m/s
zAdot1 = zAdot1 - zBdot1; % m/s

%%% Deterimine 0'-min velocity of A needed to rendezvous with B at 12'-min
dt = 12*60; % sec
wt = w*dt; % rads

yAdot_i = ((6*xA_t1*(wt - sin(wt)) - yA_t1)*w*sin(wt) - 2*w*xA_t1*(4-3*cos(wt))*(1-cos(wt)))/...
    ((4*sin(wt) - 3*wt)*sin(wt) + 4*(1-cos(wt))^2);
xAdot_i = -(w*xA_t1*(4-3*cos(wt)) + 2*(1 - cos(wt))*yAdot_i)/sin(wt);
zAdot_i = -zA_t1*w*cot(wt);

%%% Determine dV of A
dV_A = [(xAdot_i - xAdot1); (yAdot_i - yAdot1);(zAdot_i - zAdot1)];

%%% Printing result
fprintf('------------------- 3b -------------------\n')
fprintf('Necessary dV of Satellite A to achieve rendezvous trajectory:\n')
fprintf('dRSW = [%f, %f, %f] m/s\n\n', dV_A(1), dV_A(2), dV_A(3))

% ------------------------------------------------------------------------
%%% 3c
% ------------------------------------------------------------------------
%%% Determine new velocity of satellite A at 12' minutes
dt = 12*60; % sec
wt = w*dt; % rads
xAdot1 = xAdot_i*cos(wt) + (3*w*xA_t1 + 2*yAdot_i)*sin(wt); % m/s
yAdot1 = (6*w*xA_t1 + 4*yAdot_i)*cos(wt) - 2*xAdot_i*sin(wt) - (6*w*xA_t1 + 3*yAdot_i); % m/s
zAdot1 = -zA_t1*w*sin(wt) + zAdot_i*cos(wt); % m/s

%%% Printing result
fprintf('------------------- 3c -------------------\n')
fprintf('Necessary dV of Satellite A rendezvous at t=20 min:\n')
fprintf('dRSW = [%f, %f, %f] m/s\n', -xAdot1, -yAdot1, -zAdot1)


% ------------------------------------------------------------------------
%%% 3d
% ------------------------------------------------------------------------
time = 3600*5; % 2 hours (in seconds)
%%% Initializing matrices to store positional and range data
pos = zeros(time,3);
posNorm = zeros(time,1);
for dt = 1:time
wt = w*dt;
%%% Calculating position of Satellite A after entering rendezvous
%%% trajectory
xA_t = (xAdot_i/w)*sin(wt) - (3*xA_t1 + 2*yAdot_i/w)*cos(wt) + (4*xA_t1 + 2*yAdot_i/w); % m
yA_t = (6*xA_t1 + 4*yAdot_i/w)*sin(wt) + 2*xAdot_i/w*cos(wt) - (6*w*xA_t1 + 3*yAdot_i)*dt +...
    (yA_t1 - 2*xAdot_i/w); % m
zA_t = zA_t1*cos(wt) + zAdot_i/w*sin(wt); % m
pos(dt,:) = [xA_t, yA_t, zA_t];
posNorm(dt,1) = norm(pos(dt,:));
end
%%% Plotting distance between satellites over time
plot((1:time)./60,posNorm,'linewidth',2)
PlotBoi('Time, minutes','Distance Between SatA & SatB, meters')






