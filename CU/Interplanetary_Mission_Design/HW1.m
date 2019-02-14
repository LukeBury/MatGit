clear
clc
close all
mbinPath = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/mbin';
addpath(genpath(mbinPath))

% ========================================================================
%%% Importing Data
% ========================================================================
%%% General data on solar system bodies
bodies = getBodyData(mbinPath);

%%% Color options/schemes
colors = get_colors();

% ========================================================================
%%% Run Switches
% ========================================================================
run_p2 = 1;

% ========================================================================
%%% Setting up Data
% ========================================================================
%%% Loading bodies
Earth = bodies.earth;
Mars = bodies.mars;
Sun = bodies.sun;

%%% Resetting hw-specific parameters
Sun.u   = 1.327e11;        % km^3/s^2
Earth.u = 3.986e5;         % km^3/s^2
Mars.u  = 4.305e4;         % km^3/s^2

AU_km   = 1.4959787e8;     % km

Earth.a = 1 * AU_km;       % km
Mars.a  = 1.52368 * AU_km; % km

Earth.R = 6378.1363;       % km
Mars.R  = 3397.2;          % km

%%% Semi-major axis of circ orbits at Earth and Mars
r_Earth = Earth.R + 400;   % km
r_Mars  = Mars.R + 400;    % km
% ========================================================================
%%% Problem 1
% ========================================================================
fprintf('========================\nProblem 1\n========================\n')
%%% Heliocentric velocities of planets
v_Earth_HC = sqrt(Sun.u / Earth.a); % km/s
v_Mars_HC  = sqrt(Sun.u / Mars.a);  % km/s

%%% Spacecraft velocity at Earth and Mars orbits
v_circ_Earth = sqrt(Earth.u / r_Earth); % km/s
v_circ_Mars  = sqrt(Mars.u / r_Mars);   % km/s

%%% Transfer orbit velocities
a_transfer = (Earth.a + Mars.a) / 2;
v_transfer_per = visviva_v( Earth.a, a_transfer, Sun.u) % km/s
v_transfer_apo = visviva_v( Mars.a, a_transfer, Sun.u)  % km/s

%%% Spacecraft v-infinity values
vInf_Earth = v_transfer_per - v_Earth_HC; % km/s
vInf_Mars  = v_Mars_HC - v_transfer_apo;  % km/s

%%% Energy of hyperbolic orbits with respect to Earth/Mars
energy_vInf_Earth = (vInf_Earth^2) / 2; % km^2 / s^2
energy_vInf_Mars  = (vInf_Mars^2) / 2;  % km^2 / s^2

%%% Close-approach hyperbolic trajectory velocities
v_postDV1 = sqrt(2*energy_vInf_Earth + 2*Earth.u/r_Earth); % km/s
v_preDV2  = sqrt(2*energy_vInf_Mars + 2*Mars.u/r_Mars);    % km/s

%%% Delta-V magnitudes
dV_Earth = v_postDV1 - v_circ_Earth % km/s
dV_Mars  = v_preDV2 - v_circ_Mars   % km/s
dV_Total = dV_Earth + dV_Mars       % km/s

%%% Transfer time of flight
transfer_time_sec = pi * sqrt((a_transfer^3) / Sun.u) % sec
transfer_time_days = transfer_time_sec / (3600*24)    % days


% ========================================================================
%%% Problem 2
% ========================================================================
fprintf('========================\nProblem 2\n========================\n')
if run_p2 == 1
% -------------------------------------------------
% Setting initial/final states
% -------------------------------------------------
%%% Given states of planets
r0_Earth = [-578441.002878924, -149596751.684464, 0]; % km
v0_Earth = [29.7830732658560, -0.115161262358529, 0]; % km/s
X0_Earth = [r0_Earth, v0_Earth]';

rf_Mars = [-578441.618274359, 227938449.869731, 0];    % km
vf_Mars = [-24.1281802482527, -0.0612303173808154, 0]; % km/s
Xf_Mars = [rf_Mars, vf_Mars]';

%%% Initial position of Mars
r0_Mars = R3(rf_Mars, -Mars.meanMot*transfer_time_sec);

%%% Initial spacecraft state
r0_Sc = [0, -Earth.a, 0];       % km
v0_Sc = [v_transfer_per, 0, 0]; % km/s
X0_Sc = [r0_Sc, v0_Sc]';

% -------------------------------------------------
% Preparing for Integration
% -------------------------------------------------
%%% Creating time vector
t0 = 0;             % sec
tf = transfer_time_sec; % sec
time = linspace(t0,tf,1000);

%%% Choosing ode tolerance
tol = 2.22045e-14;

%%% Setting integrator options
options = odeset('RelTol',tol,'AbsTol',tol);

% -------------------------------------------------
% Integrating both cases
% -------------------------------------------------
%%% Integrating spacecraft state with 2-body-inertial integrator
[time, X_Sc_2BI] = ode113(@Int_2BI_HW1,time, X0_Sc, options, Sun.u);

%%% Integrating spacecraft state with perturbed 2BI integrator
[time, X_Sc_2BI_pert] = ode113(@Int_2BI_pert_HW1,time, X0_Sc, options, Sun.u, r0_Earth, r0_Mars, Earth, Mars);


% -------------------------------------------------
% Creating time history of planet positions
% -------------------------------------------------
r_Earth_HC = zeros(length(time),3);
r_Mars_HC  = zeros(length(time),3);

for kk = 1:length(time)
    r_Earth_HC(kk,:) = R3(r0_Earth, Earth.meanMot*time(kk));
    r_Mars_HC(kk,:)  = R3(r0_Mars, Mars.meanMot*time(kk));
end

% -------------------------------------------------
% Differencing nominal and perturbed trajectories
% -------------------------------------------------
r_diff = X_Sc_2BI_pert(:,1:3) - X_Sc_2BI(:,1:3);
v_diff = X_Sc_2BI_pert(:,4:6) - X_Sc_2BI(:,4:6);

% -------------------------------------------------
% Plotting
% -------------------------------------------------
figure; hold all
axis equal
plot3(r_Earth_HC(1,1),r_Earth_HC(1,2),r_Earth_HC(1,3),'ob','linewidth',1,'markersize',10);
plot3(r_Mars_HC(1,1),r_Mars_HC(1,2),r_Mars_HC(1,3),'or','linewidth',1,'markersize',10);
plot3(r_Earth_HC(end,1),r_Earth_HC(end,2),r_Earth_HC(end,3),'xb','linewidth',1,'markersize',10);
plot3(r_Mars_HC(end,1),r_Mars_HC(end,2),r_Mars_HC(end,3),'xr','linewidth',1,'markersize',10);
p_Earth = plot3(r_Earth_HC(:,1),r_Earth_HC(:,2),r_Earth_HC(:,3),'b','linewidth',2);
p_Mars = plot3(r_Mars_HC(:,1),r_Mars_HC(:,2),r_Mars_HC(:,3),'r','linewidth',2);
p_2BI  = plot3(X_Sc_2BI(:,1),X_Sc_2BI(:,2),X_Sc_2BI(:,3),'--k','linewidth',2);
p_pert = plot3(X_Sc_2BI_pert(:,1),X_Sc_2BI_pert(:,2),X_Sc_2BI_pert(:,3),'--m','linewidth',2);
plotBody2(bodies.sun.R,[0,0,0],colors.std.ylw,colors.std.black,1,1)
PlotBoi3('X, $km$','Y, $km$','Z, $km$',18,'LaTex')
legend([p_Earth, p_Mars, p_2BI, p_pert],'Earth','Mars','2-Body Dynamics','Perturbed Dynamics')

figure('position',[440 470 809 328])
subplot(2,2,1)
plot(time.*86400,r_diff(:,1),'linewidth',2)
PlotBoi2('','$\Delta$ X-Pos, $km$',18,'LaTex')
subplot(2,2,3)
plot(time.*86400,r_diff(:,2),'linewidth',2)
PlotBoi2('Time, $Days$','$\Delta$ Y-Pos, $km$',18,'LaTex')
subplot(2,2,2)
plot(time.*86400,v_diff(:,1),'linewidth',2)
PlotBoi2('','$\Delta$ X-Vel, $km/s$',18,'LaTex')
subplot(2,2,4)
plot(time.*86400,v_diff(:,2),'linewidth',2)
PlotBoi2('Time, $Days$','$\Delta$ Y-Vel, $km/s$',18,'LaTex')
end










% ========================================================================
%%% Functions
% ========================================================================
function [dX] = Int_2BI_HW1(t,X,u)
%%% For numerical integration of ECI state under influence of standard
%%% 2-body gravity
%%% Inputs:
%          t - time vector
%          X - initial state [6x1]
%          u - gravitational parameter of primary body

%%% Preallocate state output
dX = zeros(6,1);

%%% Distances from primary body to spacecraft
r = [X(1); X(2); X(3)]; % km
r_mag = norm(r);

%%% 2B accelerations
a_2B = r * (-u / (r_mag^3));

%%% Output the derivative of the state
dX(1:3) = [X(4); X(5); X(6)];    % km/s
dX(4:6) = [a_2B(1); a_2B(2); a_2B(3)]; % km/s^2

end

function [dX] = Int_2BI_pert_HW1(t,X,u,r0_Earth,r0_Mars,Earth,Mars)
%%% For numerical integration of ECI state under influence of standard
%%% 2-body gravity and 3rd body perturbations
%%% Inputs:
%          t - time vector
%          X - initial state [6x1]
%          u - gravitational parameter of primary body
%          u - g
%          u - g
%          u - g
%          u - g

%%% Preallocate state output
dX = zeros(6,1);

%%% Compute current HC position of Earth
rEarth_Sun = R3(r0_Earth, Earth.meanMot*t); % km
rSun_Earth = -rEarth_Sun;

%%% Compute current HC positions of Mars
rMars_Sun  = R3(r0_Mars, Mars.meanMot*t);   % km
rSun_Mars = -rMars_Sun;

%%% Distances from primary body to spacecraft
rSc_Sun = [X(1), X(2), X(3)];
rSc_Sun_mag = norm(rSc_Sun);

%%% Computing position of spacecraft with respect ot planets
rSc_Earth = rSun_Earth + rSc_Sun; % km
rSc_Mars  = rSun_Mars  + rSc_Sun; % km

%%% 2B accelerations
a_2B = rSc_Sun * (-u / (rSc_Sun_mag^3));

%%% 3B accelerations - Earth
a_3B_Earth = -Earth.u * (rSc_Earth/(norm(rSc_Earth)^3) - rSun_Earth/(norm(rSun_Earth)^3));

%%% 3B accelerations - Mars
a_3B_Mars  = -Mars.u * (rSc_Mars/(norm(rSc_Mars)^3) - rSun_Mars/(norm(rSun_Mars)^3));

%%% Output the derivative of the state
dX(1:3) = [X(4); X(5); X(6)];    % km/s
dX(4:6) = [a_2B(1) + a_3B_Earth(1) + a_3B_Mars(1);...
           a_2B(2) + a_3B_Earth(2) + a_3B_Mars(2);...
           a_2B(3) + a_3B_Earth(3) + a_3B_Mars(3)]; % km/s^2

end


