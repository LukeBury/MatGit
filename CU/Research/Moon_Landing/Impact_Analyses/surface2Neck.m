clear
clc
close all
mbinPath = '~/CU_Google_Drive/Documents/MatGit/mbin';
moonFuncsPath = '~/CU_Google_Drive/Documents/MatGit/CU/Research/Moon_Landing/Moon_Landing_funcs';
addpath(genpath(mbinPath))
addpath(genpath(moonFuncsPath))
tic

%%% This script backwards-propagates trajectories from the north pole of
%%% europa at a specified energy level (in terms of L2-flyover speed)
% ========================================================================
%%% Run Switches
% ========================================================================
storePlot_FullTraj = 1;
useEvents          = 1; % stop int on impact or interior escape

print_lowImpactL2Trajs = 1;

onlyLowAngleTrajs  = 0;
impactAngleUpperLimit_deg = 40;
% ========================================================================
%%% Free variables
% ========================================================================
%%% Number of vHat directions to test
n_vHats = 2000;

%%% How fast the SC would be traveling over L2
% dvLp_mps = 200; % Meters per second
% dvLp_mps = 400; % Meters per second
dvLp_mps = 50; % Meters per second

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

%%% Initial position from lat/lon on surface (0-lon points toward primary)
% lat0 = 90;
% lon0 = 0;

% lat0 = -5;
% lon0 = -245;

%%% Europa Plume, https://europa.nasa.gov/resources/31/source-region-for-possible-europa-plumes/
lat0 = 5;
lon0 = -20;

[rSCR0_n] = latlon2SCR(lat0, lon0, secondary.R_n);


rf_n = rSCR0_n + [1-secondary.MR, 0, 0];

%%% Create matrix of initial velocity directions (forward in time)
[vHatHem] = vHatHemisphere(n_vHats,'-y');

vfs = [];
if onlyLowAngleTrajs == 1
    for kk = 1:size(vHatHem,1)
        [impactAngle_kk] = calcImpactAngle(rSCR0_n,vHatHem(kk,:),'degrees');
        
        if impactAngle_kk < impactAngleUpperLimit_deg
            vfs = [vfs; vHatHem(kk,:)];
        end
    end
elseif onlyLowAngleTrajs == 0
    vfs = vHatHem;
end

% -------------------------------------------------
% Determing |v0| based on desired JC (in terms of L2-flyover)
% -------------------------------------------------
%%% Putting this in terms of JC
dJC_vel_kps = dvLp_mps/1000;
dJC_Lp = (dJC_vel_kps/vNorm)^2;
    
%%% Jacobi constant of Lagrange point
[JC_Lp] = JacobiConstantCalculator(secondary.MR,L123(2,:),[0,0,0]);

%%% s/c starting JC (JC_scDesired) is lower than JC_Lp because the dJC_Lp is accounted for
%%% in velocity
JC_scInitial = JC_Lp-dJC_Lp;

%%% With JC0 defined, starting velocity is a function of position. So first
%%% we must calculate the JC of the stationary starting position
JC_initialPos = JacobiConstantCalculator(secondary.MR,rf_n,[0,0,0]);

%%% Starting velocity is found from difference between s/c JC (JC_scDesired) and the
%%% JC of the stationary starting position (JC_initialPos)
dJC_forInitialVelocity = JC_initialPos - JC_scInitial;

%%% Find necessary velocity magnitude
if dJC_forInitialVelocity < 0
    warning('Spacecraft starting in a forbidden zone')
elseif dJC_forInitialVelocity >= 0
    v0Mag = sqrt(abs(dJC_forInitialVelocity));
end

% -------------------------------------------------
% Preparing for integration
% -------------------------------------------------
%%% Setting time vector
t_f = 0;
n_dt = 1000;
time0_bkwds_n = linspace(t_i,t_f,n_dt);

%%% Choosing ode tolerance
tol = 2.22045e-14;

%%% Setting integrator options
options = odeset('RelTol',tol,'AbsTol',tol);
options_ImpactEscape = odeset('Events',@event_ImpactEscape_CR3Bn,'RelTol',tol,'AbsTol',tol);

% -------------------------------------------------
% Propagating states
% -------------------------------------------------
%%% Preallocating
trajs{size(vfs,1)} = [];
L2_escapees_X0T = [];
%%% Looping through conditions
for kk = 1:size(vfs,1)
    %%% Initial state
    vf_n = vfs(kk,:).*v0Mag;
    
    Xf_n = [rf_n, vf_n]';
    
    
    %%% Setting extra parameters
    extras.u = secondary.MR;
    extras.R2_n = secondary.R_n;
    extras.L1x = L123(1,1);
    extras.L2x = L123(2,1);
    
    time_eventImpact_bkwds = [];
    X_eventImpact_bkwds = [];
    index_eventImpact = [];
    %%% Integrating
    if useEvents == 1
        [time_bkwds_n, X_BCR_bkwds_n, time_eventImpact_bkwds, X_eventImpact_bkwds, index_eventImpact] = ode113(@Int_CR3Bn,...
                time0_bkwds_n, Xf_n, options_ImpactEscape, extras);
    elseif useEvents == 0
        [time_bkwds_n, X_BCR_bkwds_n] = ode113(@Int_CR3Bn, time0_bkwds_n, Xf_n, options, extras);
    end
    
%     if X_BCR_bkwds_n(end,1) > 1.02
%         kk
%         X_eventImpact_bkwds
%     end
%     
    if useEvents == 1
        if isempty(X_eventImpact_bkwds) == 0
            if X_eventImpact_bkwds(end,1) == L123(2,1)
                L2_escapees_X0T = [L2_escapees_X0T; X_eventImpact_bkwds(end,:), time0_bkwds_n-time_eventImpact_bkwds(end)];
            end
        end
    end
    
    %%% Storing
    if storePlot_FullTraj == 0
        trajs{kk} = [X_BCR_bkwds_n(1,:); X_BCR_bkwds_n(end,:)];
    elseif storePlot_FullTraj == 1
        trajs{kk} = X_BCR_bkwds_n;
    end
end

if print_lowImpactL2Trajs == 1
    if isempty(L2_escapees_X0T) == 0
        if useEvents == 1
            for kk = 1:size(L2_escapees_X0T,1)
                fprintf('[%1.16f, %1.16f, %1.16f, %1.16f, %1.16f, %1.16f]\n', L2_escapees_X0T(kk,1:6))
            end
        end
        
        minIdx = find(L2_escapees_X0T(:,7) == min(L2_escapees_X0T(:,7)));
        fprintf('\nMinimum-time trajectory: %1.6f\n', L2_escapees_X0T(minIdx,7))
        fprintf('[%1.16f, %1.16f, %1.16f, %1.16f, %1.16f, %1.16f]\n\n', L2_escapees_X0T(minIdx,1:6))
        maxIdx = find(L2_escapees_X0T(:,7) == max(L2_escapees_X0T(:,7)));
        fprintf('\nMaximum-time trajectory: %1.6f\n', L2_escapees_X0T(maxIdx,7))
        fprintf('[%1.16f, %1.16f, %1.16f, %1.16f, %1.16f, %1.16f]\n\n', L2_escapees_X0T(maxIdx,1:6))
    end
end


% -------------------------------------------------
% Find Y-Z contour points
% -------------------------------------------------
%%% Create grid of starting locations based on y-z neck
ys = linspace(-0.02, 0.02, 500);
zs = linspace(-0.03, 0.03, 500);
[Y_yz,Z_yz] = meshgrid(ys,zs);

%%% Calculating JCs across y-z grid
JCs_yz_Lpoint = zeros(size(Y_yz));
for yk = 1:size(Y_yz,1)
    for zk = 1:size(Y_yz,2)
        %%% Zero-Velocity Curve
        zv = JacobiConstantCalculator(secondary.MR,[1-secondary.MR, Y_yz(yk,zk), Z_yz(yk,zk)] ,[0, 0, 0]);
        JCs_yz_Lpoint(yk,zk) = zv;
    end
end

%%% Get points of y-z contour in 3D space
[ yzContourPoints ] = getContourPoints( Y_yz, Z_yz, JCs_yz_Lpoint, JC_scInitial );


figure; hold all
plotBodyTexture3(secondary.R_n,[1-secondary.MR,0,0],secondary.img)
if storePlot_FullTraj == 0
    for kk = 1:length(trajs)
        plot3(trajs{kk}(2,1),trajs{kk}(2,2),trajs{kk}(2,3),'r.','markersize',10)
    end
elseif storePlot_FullTraj == 1
    for kk = 1:length(trajs)
        plot3(trajs{kk}(:,1),trajs{kk}(:,2),trajs{kk}(:,3),'r-','linewidth',0.5)
    end
end
plot3(L123(1,1),L123(1,2),L123(1,3),'^','markersize',8,'markeredgecolor',colors.std.black,'markerfacecolor',colors.std.ltgrn)
plot3(L123(2,1),L123(2,2),L123(2,3),'^','markersize',8,'markeredgecolor',colors.std.black,'markerfacecolor',colors.std.ltgrn)
% plot3(ones(size(yzContourPoints,2)),yzContourPoints(1,:),yzContourPoints(2,:),'linewidth',1.5,'color',colors.std.black)
PlotBoi3('$x_n$','$y_n$','$z_n$',16,'LaTex')
axis equal
plotCR3BP_Neck(secondary,L123,JC_scInitial,600,200,colors.std.black,2)    


toc

figure; hold all
for kk = 1:length(trajs)
maxRad = max(rowNorm(trajs{kk}(:,1:3)-[1-secondary.MR,0,0]));
plot(dot([0,0,1],vfs(kk,:)),maxRad,'b.','markersize',10)
end



figure; hold all
for kk = 1:length(trajs)
JCkk = JacobiConstantCalculator(secondary.MR,[trajs{kk}(end,1),trajs{kk}(end,2),trajs{kk}(end,3)],[trajs{kk}(end,4),trajs{kk}(end,5),trajs{kk}(end,6)]);
plot(kk,abs(JCkk - JC_scInitial),'b.','markersize',10)
end
% 
% kk = 2;
% JCs = JacobiConstantCalculator(secondary.MR,[trajs{kk}(:,1),trajs{kk}(:,2),trajs{kk}(:,3)],[trajs{kk}(:,4),trajs{kk}(:,5),trajs{kk}(:,6)]);
% figure; 
% subplot(1,2,1); hold all
% plot((JCs-JC_scInitial),'.')
% subplot(1,2,2); hold all
% plot(percentchange(JCs),'.')




























