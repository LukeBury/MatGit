clear
clc
% close all
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

n = 10;
data{n} = [];
parfor ii = 1:n
    data{ii}.d = 5 + ii;
end

dataMatrix = zeros(n,1);
for kk = 1:length(data)
    dataMatrix(kk) = data{kk}.d;
end
dataMatrix

% primary = bodies.jupiter;
% secondary = bodies.europa;
% 
% L123 = EquilibriumPoints(secondary.MR,1:3); % [3x3] of L1, L2, L3 normalized BCR coordinates
% 
% %%% Normalizing constants
% rNorm = secondary.a;         % n <-> km
% tNorm = 1/secondary.meanMot; % n <-> sec
% vNorm = rNorm / tNorm;       % n <-> km/sec
% 
% % % % X0_n = [1.02046170, 0.00483566, -0.00102579, -0.01164610, -0.00122405, -0.00248909]; % Impact
% % % % X0_n = [1.02046170, 0.00483566, -0.00102579, 0.01164610, -0.00122405, -0.00248909]; % L1 escape
% % % % X0_n = [1.005, 0.00483566, -0.00102579, 0.01164610*10, -0.00122405, -0.00248909]; % L2 escape
% % X0_n = [1.020461701526617, 0.000835866637402, 0.001140712619206, -0.002009572491451, 0.000126431584622, -0.001462926629679]';
% X0_n = [1.020461701526617, 0.000835866637402, 0.001140712619206, -0.001665385408969, -0.000000000000000, -0.001849597877215]';
% 
% %%% Selecting time vector
% t_i = 0; % sec
% t_f = 4*pi;
% dt = t_f/10000;
% time0_n = t_i:dt:t_f;
% 
% %%% Choosing ode45 tolerance
% tol = 1e-13;
% 
% %%% Setting integrator options
% options_Impact = odeset('Events',@event_ImpactEscape_CR3Bn,'RelTol',tol,'AbsTol',tol);
% 
% %%% Setting extra parameters
% extras.u = secondary.MR;
% extras.R2_n = secondary.R_n;
% extras.L1x = L123(1,1);
% extras.L2x = L123(2,1);
% 
% [time_n, X_BCR_n, time_eventImpact, X_eventImpact, index_eventImpact] = ode113(@Int_CR3Bn,...
%             time0_n, X0_n, options_Impact, extras);
%         
% % r_SCR_n = X_BCR_n(:,1:3) - [1-secondary.MR, 0, 0];
% %%% Creating SCR position
% % rImpact_SCR_n = X_eventImpact(1,1:3) - [1-secondary.MR, 0, 0];
% 
% 
% figure; hold all
% plotBodyTexture3(secondary.R_n, [1-secondary.MR, 0, 0], secondary.img)
% plot3(X_BCR_n(:,1),X_BCR_n(:,2),X_BCR_n(:,3),'linewidth',1.5)
% PlotBoi3('x','y','z',16)
% axis equal