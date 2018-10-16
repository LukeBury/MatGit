clear
clc
close all
mbinPath = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/mbin';
moonFuncsPath = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/CU/Research/Moon_Landing/Moon_Landing_funcs';
addpath(genpath(mbinPath))
addpath(genpath(moonFuncsPath))

% ========================================================================
%%% Importing Data
% ========================================================================
%%% General data on solar system bodies
bodies = getBodyData(mbinPath);

%%% Color options/schemes
colors = get_colors();

% ========================================================================
%%% Testing
% ========================================================================


primary = bodies.jupiter;
secondary = bodies.europa;

L123 = EquilibriumPoints(secondary.MR,1:3); % [3x3] of L1, L2, L3 normalized BCR coordinates

%%% Normalizing constants
rNorm = secondary.a;         % n <-> km
tNorm = 1/secondary.meanMot; % n <-> sec
vNorm = rNorm / tNorm;       % n <-> km/sec

X0_n = [1.02046170, 0.00483566, -0.00102579, -0.01164610, -0.00122405, -0.00248909]; % Impact
% X0_n = [1.02046170, 0.00483566, -0.00102579, 0.01164610, -0.00122405, -0.00248909]; % L1 escape
% X0_n = [1.005, 0.00483566, -0.00102579, 0.01164610*10, -0.00122405, -0.00248909]; % L2 escape

%%% Selecting time vector
t_i = 0; % sec
t_f = 4*pi;
dt = t_f/10000;
time0_n = t_i:dt:t_f;

%%% Choosing ode45 tolerance
tol = 1e-13;

%%% Setting integrator options
options_Impact = odeset('Events',@event_ImpactEscape_CR3Bn,'RelTol',tol,'AbsTol',tol);

%%% Setting extra parameters
extras.u = secondary.MR;
extras.R2_n = secondary.R_n;
extras.L1x = L123(1,1);
extras.L2x = L123(2,1);

[time_n, X_BCR_n, time_eventImpact, X_eventImpact, index_eventImpact] = ode113(@Int_CR3Bn,...
            time0_n, X0_n, options_Impact, extras);
        
% r_SCR_n = X_BCR_n(:,1:3) - [1-secondary.MR, 0, 0];
%%% Creating SCR position
% rImpact_SCR_n = X_eventImpact(1,1:3) - [1-secondary.MR, 0, 0];


figure; hold all
plotBodyTexture3(secondary.R_n, [1-secondary.MR, 0, 0], secondary.img)
plot3(X_BCR_n(:,1),X_BCR_n(:,2),X_BCR_n(:,3),'linewidth',1.5)
PlotBoi3('x','y','z',16)
axis equal
% 
% a = dlmread('/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/M.iGS_EurL2_200mps_1200km_156v0s_land.txt',',',1,0);
% 
% plot3(a(:,2),a(:,3),a(:,4),'r','linewidth',1.5)
% 
% diff = abs(a(2:end,2:4) - X_BCR_n(:,1:3))
% 
% 
% 
% 
% 
% %%% Creating SCR position vector
% rImpactHat_SCR_n = rImpact_SCR_n./norm(rImpact_SCR_n);
% 
% %%% Velocity unit vector at impact
% vHatImpact_n = X_eventImpact(end,4:6)./norm(X_eventImpact(end,4:6));
% 
% %%% Angle between velocity and surface
% A = R3(rImpactHat_SCR_n,pi/2);
% B = vHatImpact_n;
% impactAngle = acos(dot(A,B)/(norm(A)*norm(B)));
% if impactAngle > pi/2
%     impactAngle = pi - impactAngle;
% end
% 
% %%% Going with degrees
% impactAngle = impactAngle*180/pi
% 
% 
% 
% %%% Creating SCR position vector
% rImpactHat_SCR_n = rImpact_SCR_n./norm(rImpact_SCR_n);
% 
% %%% Velocity unit vector at impact
% vHatImpact_n = X_eventImpact(end,4:6)./norm(X_eventImpact(end,4:6));
% 
% %%% Angle between velocity and surface
% A = R3(rImpactHat_SCR_n,pi/2);
% B = vHatImpact_n;
% impactAngle = acos(dot(A,B)/(norm(A)*norm(B)));
% if impactAngle > pi/2
%     impactAngle = pi - impactAngle;
% end
% 
% %%% Going with degrees
% impactAngle = impactAngle*180/pi
% 
% 
% 
% 
% %%% Creating SCR position vector
% rImpact_SCR_n = a(end,2:4) - [1-secondary.MR, 0, 0];
% rImpactHat_SCR_n = rImpact_SCR_n./norm(rImpact_SCR_n);
% 
% %%% Velocity unit vector at impact
% vHatImpact_n = a(end,5:7)./norm(a(end,5:7));
% 
% %%% Angle between velocity and surface
% A = R3(rImpactHat_SCR_n,pi/2);
% B = vHatImpact_n;
% impactAngle = acos(dot(A,B)/(norm(A)*norm(B)));
% if impactAngle > pi/2
%     impactAngle = pi - impactAngle;
% end
% 
% %%% Going with degrees
% impactAngle = impactAngle*180/pi
% 
