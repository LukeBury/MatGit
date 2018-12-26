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
plot_fwdTrajs = 1;

%%% If bkwd trajs should stop on events
use_eventOnBkwdTraj = 0;
if use_eventOnBkwdTraj == 1
    warning('Some trajectories are getting stopped immediately at neck')
end
% ========================================================================
% Choosing data files
% ========================================================================
%%% Low-angle landing file
lowImpactAngleFile = 'F.iGS_eurL2_100mps_50km_149v0s_land.txt';
 
%%% Outputs path
MatlabOutputsPath = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/MatlabOutputs/';
        
% ========================================================================
%%% Free variables
% ========================================================================
%%% One trajectory for every (trajCutInt) will be propagated and plotted
%%% .. if trajCutInt = 1, then all trajs will be propagated/plotted
trajCutInt = 250;

%%% "Final" time
t_i = 4*pi;
 
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
primary = bodies.jupiter;
secondary = bodies.europa;
 
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
time0_bkwd_n = linspace(t_i, t_f, 1000);
 
%%% Choosing ode tolerance
tol = 1e-13;
 
%%% Setting integrator options
options = odeset('RelTol',tol,'AbsTol',tol);
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
            [time_bkwd_n, X_BCR_bkwd_n] = ode113(@Int_CR3Bn, time0_bkwd_n, lowImpactMat(kk,1:6)', options, prms);
        elseif use_eventOnBkwdTraj == 1
            [time_bkwd_n, X_BCR_bkwd_n] = ode113(@Int_CR3Bn, time0_bkwd_n, lowImpactMat(kk,1:6)', options_ImpactEscape, prms);
        end

        lowTrajs_bkwd{kk} = [X_BCR_bkwd_n, time_bkwd_n];

        if plot_fwdTrajs == 1
            %%% Propagating trajectory forward
            [time_fwd_n, X_BCR_fwd_n] = ode113(@Int_CR3Bn, fliplr(time0_bkwd_n), lowImpactMat(kk,1:6)', options_ImpactEscape, prms);

            lowTrajs_fwd{kk} = [X_BCR_fwd_n, time_fwd_n];
        end
    end
end


% -------------------------------------------------
% Plotting
% -------------------------------------------------
figure
for kk = 1:length(lowTrajs_bkwd) 
% for kk = 1:25
    if rem(kk,trajCutInt) == 0
        subplot(1,2,1); hold all
        plot3(lowTrajs_bkwd{kk}(:,1), lowTrajs_bkwd{kk}(:,2), lowTrajs_bkwd{kk}(:,3),'r','linewidth',1.5)
        if plot_fwdTrajs == 1
            plot3(lowTrajs_fwd{kk}(:,1), lowTrajs_fwd{kk}(:,2), lowTrajs_fwd{kk}(:,3),'b','linewidth',1.5)
        end

        subplot(1,2,2); hold all
        p_bkwd = plot3(lowTrajs_bkwd{kk}(:,1), lowTrajs_bkwd{kk}(:,2), lowTrajs_bkwd{kk}(:,3),'r','linewidth',1.5);
        if plot_fwdTrajs == 1
            p_fwd = plot3(lowTrajs_fwd{kk}(:,1), lowTrajs_fwd{kk}(:,2), lowTrajs_fwd{kk}(:,3),'b','linewidth',1.5);
        end
    end
end

subplot(1,2,1); title('Full System')
plotBodyTexture3(secondary.R_n, [1-secondary.MR, 0, 0], secondary.img)
plotBodyTexture3(primary.R/rNorm, [-secondary.MR, 0, 0], primary.img)
PlotBoi3('$x_n$','$y_n$','$z_n$',16,'LaTex')
view(0,90)
axis equal
subplot(1,2,2); title('Zoomed to Secondary')
plotBodyTexture3(secondary.R_n, [1-secondary.MR, 0, 0], secondary.img)
plot([L123(2,1), L123(2,1)], [-0.01, 0.01],'color',colors.std.ltgrn,'linewidth',1.5)
plot3(L123(1:2,1), L123(1:2,2), L123(1:2,3),'^','markeredgecolor',colors.std.mag,'markerfacecolor',colors.std.ltgrn,'markersize',10)
PlotBoi3('$x_n$','$y_n$','$z_n$',16,'LaTex')
view(0,90)
axis equal
xlim([L123(1,1)-2*secondary.R_n, L123(2,1)+2*secondary.R_n])
ylim([-10*secondary.R_n, 10*secondary.R_n])
if plot_fwdTrajs == 1
    legend([p_bkwd, p_fwd],'Backward from X0','Forward from X0')
else
    legend(p_bkwd,'Backward from X0')
end




















