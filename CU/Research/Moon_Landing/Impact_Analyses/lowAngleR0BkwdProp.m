clear
clc
close all
mbinPath = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/mbin';
moonFuncsPath = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/CU/Research/Moon_Landing/Moon_Landing_funcs';
addpath(genpath(mbinPath))
addpath(genpath(moonFuncsPath))
tic
 
% ========================================================================
%%% Run Switches
% ========================================================================
plot_fwdTrajs  = 0;
plot_L1escapee = 1;

%%% If bkwd trajs should stop on events
use_eventOnBkwdTraj = 1;
if use_eventOnBkwdTraj == 1
    warning('Some trajectories are getting stopped immediately at neck')
end
% ========================================================================
% Choosing data files
% ========================================================================
%%% Low-angle landing file
% lowImpactAngleFile = 'F.iGS_eurL2_50mps_50km_149v0s_land.txt';
% lowImpactAngleFile = 'F.iGS_eurL2_100mps_50km_149v0s_land.txt';
% lowImpactAngleFile = 'F.iGS_eurL2_200mps_50km_149v0s_land.txt';
% lowImpactAngleFile = 'F.iGS_eurL2_300mps_50km_149v0s_land.txt';

% lowImpactAngleFile = 'F.iGS_encL2_39mps_4km_149v0s_land.txt';
lowImpactAngleFile = 'F.iGS_encL2_39mps_4km_149v0s_land.txt';


%%% Outputs path
MatlabOutputsPath = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/';
        
% ========================================================================
%%% Free variables
% ========================================================================
%%% One trajectory for every (trajCutInt) will be propagated and plotted
%%% .. if trajCutInt = 1, then all trajs will be propagated/plotted
% trajCutInt = 25; % 200
trajCutInt = 50; % 200
%%% "Final" time
t_i = 8*pi;
 
% ========================================================================
%%% Importing Data
% ========================================================================
%%% General data on solar system bodies
bodies = getBodyData(mbinPath);
 
%%% Color options/schemes
colors = get_colors();
 
% ========================================================================
%%% Running
% ========================================================================
% -------------------------------------------------
% Setting conditions
% -------------------------------------------------
%%% bodies
primary = bodies.saturn;
secondary = bodies.enceladus;
 
%%% Normalizing constants
rNorm = secondary.a;         % n <-> km
tNorm = 1/secondary.meanMot; % n <-> sec
vNorm = rNorm / tNorm;       % n <-> km/sec
 
%%% Lagrange Points
L123 = EquilibriumPoints(secondary.MR,1:3); % [3x3] of L1, L2, L3 normalized BCR coordinates
 
% -------------------------------------------------
% Load Low-Impact-Angle Data
% -------------------------------------------------
%%% Load low impact data
lowImpactMat = dlmread([MatlabOutputsPath,lowImpactAngleFile],',',1,0);
nLowImpactTraj = size(lowImpactMat,1);
if lowImpactMat(end,:) == 0
    nLowImpactTraj = 0;
end

% -------------------------------------------------
% Propagate low-impact-angle trajectories
% -------------------------------------------------
%%% Preallocating
lowTrajs_bkwd{nLowImpactTraj} = [];
lowTrajs_fwd{nLowImpactTraj} = [];
 
%%% Setting time vector
t_f = 0; % (ending at 0 ... propagating backwards)
time0_bkwd_n = linspace(t_i, t_f, 10000);
 
%%% Choosing ode tolerance
tol = 1e-13;
 
%%% Setting integrator options
options = odeset('RelTol',tol,'AbsTol',tol);
options_Impact = odeset('Events',@event_Impact_CR3Bn,'RelTol',tol,'AbsTol',tol);
options_ImpactOrLeftXAxis = odeset('Events',@event_ImpactOrLeftXAxis_CR3Bn,'RelTol',tol,'AbsTol',tol);
options_LeftXAxisFromBottom = odeset('Events',@event_LeftXAxisFromBottom_CR3Bn,'RelTol',tol,'AbsTol',tol);
options_ImpactEscape = odeset('Events',@event_ImpactEscape_CR3Bn,'RelTol',tol,'AbsTol',tol);

%%% Setting necessary parameters for integration
prms.u = secondary.MR;
prms.R2_n = secondary.R_n;
prms.L1x = L123(1,1);
prms.L2x = L123(2,1);

for kk = 1:size(lowImpactMat,1)
% for kk = 1:25
    if rem(kk,trajCutInt) == 0
        %%% Propagating trajectory backward
        if use_eventOnBkwdTraj == 0
            [time_bkwd_n, X_BCR_bkwd_n] = ode113(@Int_CR3Bn, time0_bkwd_n, lowImpactMat(kk,2:7)', options, prms);
            
            lowTrajs_bkwd{kk} = [X_BCR_bkwd_n, time_bkwd_n];
        elseif use_eventOnBkwdTraj == 1
            [time_bkwd_n, X_BCR_bkwd_n, time_eventImpact, X_eventImpact, index_eventImpact] = ode113(@Int_CR3Bn, time0_bkwd_n, lowImpactMat(kk,2:7)', options_ImpactOrLeftXAxis, prms);
            
%             if isempty(time_eventImpact) == 0
%                 if X_eventImpact(end,1) < 0
%                     989
%                     lowTrajs_bkwd{kk} = [X_BCR_bkwd_n(1:index_eventImpact,:), time_bkwd_n(1:index_eventImpact)];
%                 else
%                     lowTrajs_bkwd{kk} = [X_BCR_bkwd_n, time_bkwd_n];
%                 end
%             end
            
            if abs(X_BCR_bkwd_n(end,1)) < 1e-10
                X0ofL1escapee = X_BCR_bkwd_n(1,1:6)';
                [time_bkwd_extra_n, X_BCR_bkwd_extra_n, time_extra_eventImpact, X_extra_eventImpact, index_extra_eventImpact] = ode113(@Int_CR3Bn, time0_bkwd_n, X_BCR_bkwd_n(end,1:6)', options_LeftXAxisFromBottom, prms);
                989
                lowTrajs_bkwd{kk} = [[X_BCR_bkwd_n; X_BCR_bkwd_extra_n], [time_bkwd_n; time_bkwd_extra_n]];
            else
                lowTrajs_bkwd{kk} = [X_BCR_bkwd_n, time_bkwd_n];
            end


            
        end

        

        if plot_fwdTrajs == 1
            %%% Propagating trajectory forward
            [time_fwd_n, X_BCR_fwd_n] = ode113(@Int_CR3Bn, fliplr(time0_bkwd_n), lowImpactMat(kk,2:7)', options_ImpactEscape, prms);

            lowTrajs_fwd{kk} = [X_BCR_fwd_n, time_fwd_n];
        end
        
        
        
    end
end

if plot_L1escapee == 1
    X0_L1escapee = [1.020461701526617;-0.002950041600261;0.000897404591331;-0.005397473677259;-0.012122924365303;0.002820666556784];
    [time_bkwd_n, X_BCR_bkwd_L1escapee_n] = ode113(@Int_CR3Bn, time0_bkwd_n, X0_L1escapee, options, prms);
    figure(3); hold all
    
    plot3(X_BCR_bkwd_L1escapee_n(:,1), X_BCR_bkwd_L1escapee_n(:,2), X_BCR_bkwd_L1escapee_n(:,3),'r','linewidth',1.0)
    plotBodyTexture3(secondary.R_n, [1-secondary.MR, 0, 0], secondary.img)
    % plotBodyTexture3(primary.R/rNorm, [-secondary.MR, 0, 0], primary.img)
    plot([L123(2,1), L123(2,1)], [-0.01, 0.01],'color',colors.std.black,'linewidth',1.5)
    plot3(L123(1:2,1), L123(1:2,2), L123(1:2,3),'^','markeredgecolor',colors.std.black,'markerfacecolor',colors.std.grn,'markersize',10)

    PlotBoi3('$x_n$','$y_n$','$z_n$',16,'LaTex')
    view(0,90)
    axis equal
    xlim([L123(1,1)-2*secondary.R_n, L123(2,1)+2*secondary.R_n])
    ylim([-10*secondary.R_n, 10*secondary.R_n])

end
% -------------------------------------------------
% Plotting
% -------------------------------------------------
figure(1); hold all
% figure(2); hold all
for kk = 1:length(lowTrajs_bkwd) 
    if rem(kk,trajCutInt) == 0
        figure(1)
        plot3(lowTrajs_bkwd{kk}(:,1), lowTrajs_bkwd{kk}(:,2), lowTrajs_bkwd{kk}(:,3),'r','linewidth',1.0)
        if plot_fwdTrajs == 1
            plot3(lowTrajs_fwd{kk}(:,1), lowTrajs_fwd{kk}(:,2), lowTrajs_fwd{kk}(:,3),'b','linewidth',1.0)
        end

%         figure(2)
%         p_bkwd = plot3(lowTrajs_bkwd{kk}(:,1), lowTrajs_bkwd{kk}(:,2), lowTrajs_bkwd{kk}(:,3),'r','linewidth',1.0);
%         if plot_fwdTrajs == 1
%             p_fwd = plot3(lowTrajs_fwd{kk}(:,1), lowTrajs_fwd{kk}(:,2), lowTrajs_fwd{kk}(:,3),'b','linewidth',1.0);
%         end
    end
end

figure(1)
plotBodyTexture3(secondary.R_n, [1-secondary.MR, 0, 0], secondary.img)
plotBodyTexture3(primary.R/rNorm, [-secondary.MR, 0, 0], primary.img)
PlotBoi3('$x_n$','$y_n$','$z_n$',16,'LaTex')
view(0,90)
axis equal

% figure(2)
% plotBodyTexture3(secondary.R_n, [1-secondary.MR, 0, 0], secondary.img)
% plot([L123(2,1), L123(2,1)], [-0.01, 0.01],'color',colors.std.black,'linewidth',1.5)
% plot3(L123(1:2,1), L123(1:2,2), L123(1:2,3),'^','markeredgecolor',colors.std.black,'markerfacecolor',colors.std.grn,'markersize',10)
% PlotBoi3('$x_n$','$y_n$','$z_n$',16,'LaTex')
% view(0,90)
% axis equal
% xlim([L123(1,1)-2*secondary.R_n, L123(2,1)+2*secondary.R_n])
% ylim([-10*secondary.R_n, 10*secondary.R_n])
% % if plot_fwdTrajs == 1
% %     legend([p_bkwd, p_fwd],'Backward from X0','Forward from X0')
% % else
% %     legend(p_bkwd,'Backward from X0')
% % end







function [value, isterminal, direction] = event_ImpactOrLeftXAxis_CR3Bn(t,X,prms)
%%% Designed for CR3BP
%%% Event function watching for when "value" = 0
%%% Inputs:
%          t - normalized time vector
%          X - initial state [6x1]
%          prms - (u, R2_n)
%%% Distances to primary (1) and secondary (2) bodies
r2 = sqrt((X(1)+prms.u-1)^2 + X(2)^2 + X(3)^2);

yVal = 10;
if X(1) < 0
%     if abs(X(2)) < 0.01
%         if sign(X(5)) == 1
%             testVal = -X(2);
%         elseif sign(X(5)) == -1
%             testVal = X(2);
%         end
%     end

yVal = X(2);

end
value = [r2 - prms.R2_n, yVal]; % When the surface is impacted
isterminal = [1, 1]; % stops the integration
direction = [-1, -1]; % negative direction only
end




function [value, isterminal, direction] = event_LeftXAxisFromBottom_CR3Bn(t,X,prms)
%%% Designed for CR3BP
%%% Event function watching for when "value" = 0
%%% Inputs:
%          t - normalized time vector
%          X - initial state [6x1]
%          prms - (u, R2_n)

negativeYVal = -X(2);

value = [negativeYVal]; % When the surface is impacted
isterminal = [1]; % stops the integration
direction = [-1]; % negative direction only
end









