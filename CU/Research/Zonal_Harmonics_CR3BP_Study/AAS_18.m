clear
clc
close all
mbinPath = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/mbin';
moonFuncsPath = '/Users/lukebury/CU_Google_Drive/Documents/MatGit/CU/Research/Moon_Landing/Moon_Landing_funcs';
addpath(genpath(mbinPath))
addpath(genpath(moonFuncsPath))

% ========================================================================
%%% Run Switches
% ========================================================================
run_TrajectorySim              = 1;
run_POcomparison               = 0;
run_J2J4ForceComp              = 0;
run_J2ForceCompMesh            = 0;
run_EquillibriumMovementTable  = 0;
run_zeroVelocity               = 0;
run_multSimulation             = 0;
        %%% Secondary Functions
        firstPeriapsisAnalysis_mult = 1;
        %%% Plots and Extras
        plotBCR_n_mult  = 1;
        fpaVSspeed_mult = 0;
        latlonPlot_mult = 0;
run_SurfaceJCAnalysis          = 0;

% ========================================================================
%%% Importing Data
% ========================================================================
%%% General data on solar system bodies
bodies = getBodyData(mbinPath);

%%% Color options/schemes
colors = get_colors();

% ========================================================================
%%% Trajectories
% ========================================================================
if run_TrajectorySim == 1
% -------------------------------------------------
% ICS
% -------------------------------------------------
% for elements with two rows - top row is for run without J2, bottom row is
% for run with J2
ICs = {'1:Primary','2:Secondary','3:L-Pt [2x1]-(1,2,3,4(1-J2),5(2-J2),6(3-J2))','4:Position Nudge Vector','5:Velocity Nudge Vector','6:dt','7:tf', '8:description';...
    % IC_i = 2; exaggerated differences
    bodies.jupiter, bodies.europa, [2;2], [0,0,0;0,0,0], [-1, -20, 0;-1, -20, 0], 60, 3600*24*3.9,'exaggerated differences';... % IC_i = 2
    
    % IC_i = 3; from noJ2 L2, no J2 barely moves and J2 impacts
    bodies.jupiter, bodies.europa, [2;2], [0,0,0;0,0,0], [0,0,0;0,0,0],               10, 3600*24*1.57,'from noJ2 L2, no J2 barely moves and J2 impacts';... % IC_i = 3
    
    % IC_i = 4; from J2 L2, J2 barely moves and no J2 departs
    bodies.jupiter, bodies.europa, [5;5], [0,0,0;0,0,0], [0,0,0;0,0,0],               600, 3600*24*5,'from J2 L2, J2 barely moves and no J2 departs';... % IC_i = 4
    
    % IC_i = 5; from no J2 L2, J2 orbits and impacts, no J2 departs
    bodies.jupiter, bodies.europa, [2;2], [-1,1,0;-1,1,0], [1,0,0;1,0,0],             60, 3600*24*5,'from no J2 L2, J2 orbits and impacts, no J2 departs';... % IC_i = 5
    
    % IC_i = 6; impacts w/o J2, orbits w/ J2
    bodies.jupiter, bodies.europa, [2;2], [-1.238964789451243,0.29555244,0.013817949;-1.238964789451243,0.29555244,0.013817949]*(1e4), [-836.97,-10,0;-836.97,-10,0], 10, 3600*2,'impacts w/o J2, orbits w/ J2';... % IC_i = 6
    
    % IC_i = 7; each start from respective L2 and barely move
    bodies.jupiter, bodies.io, [2;5], [0,0,0;0,0,0], [0,0,0;0,0,0],                                                       600, 3600*24*3,'each start from respective L2 and barely move';... % IC_i = 7
    
    % IC_i = 8; J2 impacts, No J2 escapes
    bodies.jupiter, bodies.europa, [2;2], [-7.505614e3,0,0;-7.505614e3,0,0], [0,0,716.424;0,0,716.424].*1.098,             60, 3600*24*6.5,'J2 impacts, No J2 escapes';... % IC_i = 8
    
    % IC_i = 9; J2 invariant complex trajectory
    bodies.jupiter, bodies.europa, [2;2], [-1.0627214e4,0,0;-1.0627214e4,0,0], [0,0,1.013177e3;0,0,1.013177e3].*1.2        60, 3600*24*4,'J2 invariant complex trajectory';... % IC_i = 9
    
    % IC_i = 10; chaotic difference from IC 11
    bodies.jupiter, bodies.europa, [2;2], [-1.0627214e4,0,0;-1.0627214e4,0,0], [0,0,1.013177e3;0,0,1.013177e3].*1.23       60, 3600*24*4,'chaotic difference from IC 11';... % IC_i = 10
    
    % IC_i = 11; chaotic difference from IC 10
    bodies.jupiter, bodies.europa, [2;2], [-1.0627214e4,0,0;-1.0627214e4,0,0], [0,0,1.013177e3;0,0,1.013177e3].*1.24       60, 3600*24*8,'chaotic difference from IC 10';... % IC_i = 11

%     % IC_i = 12
%     bodies.jupiter, bodies.europa, [2;2], [.5-bodies.europa.MR,sqrt(3)/2,0;0,0,0], [.5-bodies.europa.MR,sqrt(3)/2,0;0,0,0]       60, 3600*24*8,'testing';... % IC_i = 12
   };

IC_i = 8; % For journal paper - 2,3,5,8,9

% 2 
% 3,4
% 9, 10 - J2 invariant
% -------------------------------------------------
% Setup
% -------------------------------------------------
%%% Selecting bodies
primary   = ICs{IC_i,1};
secondary = ICs{IC_i,2};

%%% Normalizing constants
rNorm = secondary.a;         % n <-> km
tNorm = 1/secondary.meanMot; % n <-> sec
vNorm = rNorm / tNorm;       % n <-> km/sec

%%% Acquire Lagrange points
L123_noJ2n = EquilibriumPoints(secondary.MR,1:3); % [3x3] of L1, L2, L3 normalized BCR coordinates
L123_J2n = EquilibriumPoints_J2(secondary.MR,primary.J2,0,primary.R/rNorm,secondary.R_n);
L123 = [L123_noJ2n; L123_J2n]; % stacking

%%% Nudges from a set state
positionNudge   = ICs{IC_i,4}(1,:)./rNorm; % enter in km ... normalized
positionNudgeJ2 = ICs{IC_i,4}(2,:)./rNorm; % enter in km ... normalized
velocityNudge   = ICs{IC_i,5}(1,:)./1000./vNorm; % enter in m/s ... normalized
velocityNudgeJ2 = ICs{IC_i,5}(2,:)./1000./vNorm; % enter in m/s ... normalized

%%% Initial States (w/o J2)
x0_n = L123(ICs{IC_i,3}(1,1),:) + positionNudge; % equilibrium point plus nudge
% x0_n = [1.002, 0.004404, 0.0002059];
xdot0_n = velocityNudge;
X0_n = [x0_n xdot0_n]';

%%% Initial States (w/ J2)
x0_J2n = L123(ICs{IC_i,3}(2,1),:) + positionNudgeJ2; % equilibrium point plus nudge
% x0_J2n = [1.002, 0.004404, 0.0002059];
xdot0_J2n = velocityNudgeJ2;
X0_J2n = [x0_J2n xdot0_J2n]';

%%% Setting time vector and normalizing 
ti = 0; % sec
dt = ICs{IC_i,6};
tf = ICs{IC_i,7};
time0_n = [ti:dt:tf] ./ tNorm;

%%% Choosing ode45 tolerance
tol = 1e-12; % journal - this was 10

%%% Setting integrator options
options = odeset('Events',@event_Impact_CR3Bn,'RelTol',tol,'AbsTol',tol);
options_J2 = odeset('Events',@event_Impact_CR3Bn_J2,'RelTol',tol,'AbsTol',tol);

%%% Setting prms
prms.u = secondary.MR;
prms.R1_n = primary.R/rNorm;
prms.R2_n = secondary.R_n;
prms.J21 = primary.J2;
prms.J22 = 0;

%%% Propagating the States without J2
[time_n, X_BCR_n] = ode45(@Int_CR3Bn, time0_n, X0_n, options, prms);
    
%%% Propagating the States with J2
[time_J2_n, X_BCR_J2_n] = ode45(@Int_CR3Bn_J2, time0_n, X0_J2n, options_J2, prms);


%%% Determining shorter simulation, and appending NaNs for plotting/comparison purposes
if length(time_n) < length(time_J2_n) || length(time_n) == length(time_J2_n)
    timeShort = time_n; timeLong = time_J2_n;
    diff = size(X_BCR_J2_n,1) - size(X_BCR_n,1);
    X_BCR_n = [X_BCR_n;NaN(diff,6)];
elseif length(time_n) > length(time_J2_n)
    timeShort = time_J2_n; timeLong = time_n;
    diff = size(X_BCR_n,1) - size(X_BCR_J2_n,1);
    X_BCR_J2_n = [X_BCR_J2_n;NaN(diff,6)];
end

%%% Finding and printing minimum altitudes
minAlt_i = find(rowNorm((X_BCR_n(:,1:3) - [1-secondary.MR,0,0])) == min(rowNorm((X_BCR_n(:,1:3) - [1-secondary.MR,0,0]))));
minAltJ2_i = find(rowNorm((X_BCR_J2_n(:,1:3) - [1-secondary.MR,0,0])) == min(rowNorm((X_BCR_J2_n(:,1:3) - [1-secondary.MR,0,0]))));
% minAlt = (norm((X_BCR_n(minAlt_i,1:3) - [1-secondary.MR,0,0]))-secondary.R_n)*rNorm
% minAltJ2 = (norm((X_BCR_J2_n(minAltJ2_i,1:3) - [1-secondary.MR,0,0]))-secondary.R_n)*rNorm

% -------------------------------------------------
% Altitudes & Lat/Lon
% -------------------------------------------------
%%% No-J2 Loop
altitudes = zeros(size(X_BCR_n,1),1);
latLons   = zeros(size(X_BCR_n,1),2);
for kk = 1:size(X_BCR_n,1)
    %%% Turning BCR to SCR
    r_SCR = X_BCR_n(kk,1:3) - [1-secondary.MR, 0, 0];
    
    %%% Calculating altitude
    altitudes(kk) = (norm(r_SCR) - secondary.R_n).*rNorm;
    
    %%% Calculating latitude and longitude
    [lat, lon] = ECEF2latlon(r_SCR,'degrees');
    if mod(kk,20) == 0
        latLons(kk,:) = [lat, lon];
    else
        latLons(kk,:) = [NaN, NaN];
    end
    
end

% %%% J2 Loop
% for kk = 1:size(X_BCR_J2_n,1)
%     d
% end


% -----------------------------
% Rotating Plot
% -----------------------------
% figure
% figure('position',[-2285 62 1078 296])
% subplot(3,2,[1,3,5]); hold all
% plot3(X_BCR_n(:,1),X_BCR_n(:,2),X_BCR_n(:,3),'b','linewidth',2)
% plot3(X_BCR_J2_n(:,1),X_BCR_J2_n(:,2),X_BCR_J2_n(:,3),'r','linewidth',2)
% plot3(X_BCR_n(1,1),X_BCR_n(1,2),X_BCR_n(1,3),'bo','linewidth',2,'markersize',8)
% plot3(X_BCR_J2_n(1,1),X_BCR_J2_n(1,2),X_BCR_J2_n(1,3),'ro','linewidth',2,'markersize',8)
% plot3(X_BCR_n(end,1),X_BCR_n(end,2),X_BCR_n(end,3),'bx','linewidth',2,'markersize',8)
% plot3(X_BCR_J2_n(end,1),X_BCR_J2_n(end,2),X_BCR_J2_n(end,3),'rx','linewidth',2,'markersize',8)
% % plotBodyTexture3(secondary.R_n, [1-secondary.MR, 0, 0],secondary.img);
% % plotBody3(secondary.R_n, [1-secondary.MR, 0, 0], secondary.color)
% PlotBoi3('$x_n$','$y_n$','$z_n$',20,'LaTex')
% legend('Without J_2','W/ J2')
% axis equal
% view(0,90)
% grid on
% subplot(3,2,2)
% plot(timeLong,(X_BCR_J2_n(:,1)-X_BCR_n(:,1)).*rNorm,'k','linewidth',2)
% PlotBoi2('','$\Delta$x, km',20,'LaTex')
% subplot(3,2,4)
% plot(timeLong,(X_BCR_J2_n(:,2)-X_BCR_n(:,2)).*rNorm,'k','linewidth',2)
% PlotBoi2('','$\Delta$y, km',20,'LaTex')
% subplot(3,2,6)
% plot(timeLong,(X_BCR_J2_n(:,3)-X_BCR_n(:,3)).*rNorm,'k','linewidth',2)
% PlotBoi2('Normalized Time','$\Delta$z, km',20,'LaTex')
figure('position',[-2285 62 1078 296])
subplot(3,2,[1,3,5]); hold all
plot3(X_BCR_n(:,1),X_BCR_n(:,2),X_BCR_n(:,3),'b','linewidth',2)
plot3(X_BCR_J2_n(:,1),X_BCR_J2_n(:,2),X_BCR_J2_n(:,3),'r','linewidth',2)
plot3(X_BCR_n(1,1),X_BCR_n(1,2),X_BCR_n(1,3),'bo','linewidth',2,'markersize',8)
plot3(X_BCR_J2_n(1,1),X_BCR_J2_n(1,2),X_BCR_J2_n(1,3),'ro','linewidth',2,'markersize',8)
plot3(X_BCR_n(length(timeShort),1),X_BCR_n(length(timeShort),2),X_BCR_n(length(timeShort),3),'bx','linewidth',2,'markersize',8)
plot3(X_BCR_J2_n(length(timeShort),1),X_BCR_J2_n(length(timeShort),2),X_BCR_J2_n(length(timeShort),3),'rx','linewidth',2,'markersize',8)
plotBodyTexture3(secondary.R_n, [1-secondary.MR, 0, 0],secondary.img);
% plotBody3(secondary.R_n, [1-secondary.MR, 0, 0], secondary.color)
PlotBoi3('$x_n$','$y_n$','$z_n$',20,'LaTex')
legend('Without J_2','W/ J2')
axis equal
view(0,90)
grid on
subplot(3,2,2)
plot(timeShort,(X_BCR_J2_n(1:length(timeShort),1)-X_BCR_n(1:length(timeShort),1)).*rNorm,'k','linewidth',2)
PlotBoi2('','$\Delta$x, km',20,'LaTex')
subplot(3,2,4)
plot(timeShort,(X_BCR_J2_n(1:length(timeShort),2)-X_BCR_n(1:length(timeShort),2)).*rNorm,'k','linewidth',2)
PlotBoi2('','$\Delta$y, km',20,'LaTex')
subplot(3,2,6)
plot(timeShort,(X_BCR_J2_n(1:length(timeShort),3)-X_BCR_n(1:length(timeShort),3)).*rNorm,'k','linewidth',2)
PlotBoi2('Normalized Time','$\Delta$z, km',20,'LaTex')

% figure; hold all
figure('position',[-2285 62 1078 296]); hold all
plot3(X_BCR_n(:,1),X_BCR_n(:,2),X_BCR_n(:,3),'b','linewidth',2)
plot3(X_BCR_J2_n(:,1),X_BCR_J2_n(:,2),X_BCR_J2_n(:,3),'r','linewidth',2)
plot3(X_BCR_n(1,1),X_BCR_n(1,2),X_BCR_n(1,3),'bo','linewidth',2,'markersize',8)
plot3(X_BCR_J2_n(1,1),X_BCR_J2_n(1,2),X_BCR_J2_n(1,3),'ro','linewidth',2,'markersize',8)
plot3(X_BCR_n(end,1),X_BCR_n(end,2),X_BCR_n(end,3),'bx','linewidth',2,'markersize',8)
plot3(X_BCR_J2_n(end,1),X_BCR_J2_n(end,2),X_BCR_J2_n(end,3),'rx','linewidth',2,'markersize',8)
plotBodyTexture3(secondary.R_n, [1-secondary.MR, 0, 0],secondary.img);
% plotBody3(secondary.R_n, [1-secondary.MR, 0, 0], secondary.color)
PlotBoi3('$x_n$','$y_n$','z_n',20,'LaTex')
legend('Without J_2','With J_2')
axis equal
view(0,90)
grid on

% % -----------------------------
% % Plotting Differences in Position
% % -----------------------------
% figure;
% subplot(3,1,1);
% plot(timeLong.*tNorm./3600,(X_BCR_J2_n(1:length(timeLong),1)-X_BCR_n(1:length(timeLong),1)).*rNorm,'linewidth',2);
% PlotBoi2('','x-Diff, km',12)
% subplot(3,1,2);
% plot(timeLong.*tNorm./3600,(X_BCR_J2_n(1:length(timeLong),2)-X_BCR_n(1:length(timeLong),2)).*rNorm,'linewidth',2);
% PlotBoi2('','y-Diff, km',12)
% subplot(3,1,3);
% plot(timeLong.*tNorm./3600,(X_BCR_J2_n(1:length(timeLong),3)-X_BCR_n(1:length(timeLong),3)).*rNorm,'linewidth',2);
% PlotBoi2('Time, hr','z-Diff, km',12)

end % run_trajectorySim

% ========================================================================
%%% Plotting/Comparing Periodic Orbits
% ========================================================================
if run_POcomparison == 1
% figure('position',[-2331 387 788 187]); hold all
% -----------------------------
% Setup
% -----------------------------
% %%% Europa - L1 PO
% primary = bodies.jupiter;
% secondary = bodies.europa;
% % r0_n = [0.981986265631928; -0.000000001419209; -0.000000009598778];
% % v0_n = [-0.000000008769559; 0.005535983529929; 0.037472301863908];
% % tP = 3.388015997951395;
% r0_n = [0.981460590439603; -0.000000000779596; -0.000000006391289];
% v0_n = [-0.000000004926460; 0.003962110642756; 0.032534117970616];
% tP = 3.298067769777886;
% % % % % r0_n = [0.978809538796240; 0.000000345420468; 0];
% % % % % v0_n = [0.000000264003012; 0.006831168844117; 0];
% % % % % tP = 2.995763088664896;
% % % % % % r0_n = [0.978547828848621; 0.000000678875551; 0];
% % % % % % v0_n = [0.000000540270016; 0.008828949995968; 0];
% % % % % % tP = 3.000336746617772;
% X0_n = [r0_n; v0_n];
% t_J2 = 2.6;
% Lpoint = 1;

% %%% Europa - L2 PO
% primary = bodies.jupiter;
% secondary = bodies.europa;
% r0_n = [1.020027801952722; -0.000000000182340; 0];
% v0_n = [-0.000000000108287; 0.002758777827967; 0];
% X0_n = [r0_n; v0_n];
% tP = 3.077407327172548;
% t_J2 = 2.8;
% Lpoint = 2;


% %%% Enceladus - L1 PO
% primary = bodies.saturn;
% secondary = bodies.enceladus;
% r0_n = [0.996481316200089; 0.000000000645449; 0.000000004116311];
% v0_n = [0.000000003910626; 0.001173058522562; 0.007478505501771];
% tP = 3.432012972191030;
% % r0_n = [0.996321445117066; 0.000000000256312; 0.000000002204430];
% % v0_n = [0.000000001605509; 0.000697606590161; 0.005993809867779];
% % tP = 3.301922568761066;
% X0_n = [r0_n; v0_n];
% t_J2 = 1.4;
% Lpoint = 1;

%%% Enceladus - L2 PO
primary = bodies.saturn;
secondary = bodies.enceladus;
r0_n = [1.003910955372229; 0.000000000371790; 0];
v0_n = [0.000000000220842; 0.000523232459032; 0];
X0_n = [r0_n; v0_n];
tP = 3.042471396742148;
t_J2 = 1.4;
Lpoint = 2;

%%% Normalizing constants
rNorm = secondary.a;         % n <-> km
tNorm = 1/secondary.meanMot; % n <-> sec
vNorm = rNorm / tNorm;       % n <-> km/sec

%%% Acquire Lagrange points
L123_noJ2n = EquilibriumPoints(secondary.MR,1:3); % [3x3] of L1, L2, L3 normalized BCR coordinates

%%% Setting time vector and normalizing 
ti = 0; % sec
dt = tP/10000;
tf = tP;
time0_n = [ti:dt:tf];
tf_J2 = t_J2;
time0_J2_n = [ti:dt:tf_J2];

%%% Choosing ode45 tolerance
tol = 1e-10;

%%% Setting integrator options
options = odeset('RelTol',tol,'AbsTol',tol);

%%% Setting prms
prms.u = secondary.MR;
prms.R1_n = primary.R/rNorm;
prms.R2_n = secondary.R_n;
prms.J21 = primary.J2;
prms.J22 = 0;

%%% Propagating the States without J2
[time_n, X_BCR_n] = ode45(@Int_CR3Bn, time0_n, X0_n, options, prms);
    
%%% Propagating the States with J2
[time_J2_n, X_BCR_J2_n] = ode45(@Int_CR3Bn_J2, time0_J2_n, X0_n, options, prms);

% -----------------------------
% Plotting
% -----------------------------
figure(1); hold all
p1 = plot3(X_BCR_n(:,1),X_BCR_n(:,2),X_BCR_n(:,3),'b','linewidth',2);
p2 = plot3(X_BCR_J2_n(:,1),X_BCR_J2_n(:,2),X_BCR_J2_n(:,3),'r','linewidth',2);
% plot3(X_BCR_n(1,1),X_BCR_n(1,2),X_BCR_n(1,3),'bo','linewidth',2,'markersize',8)
% plot3(X_BCR_J2_n(1,1),X_BCR_J2_n(1,2),X_BCR_J2_n(1,3),'ro','linewidth',2,'markersize',8)
plot3(X_BCR_n(end,1),X_BCR_n(end,2),X_BCR_n(end,3),'bx','linewidth',2,'markersize',8)
plot3(X_BCR_J2_n(end,1),X_BCR_J2_n(end,2),X_BCR_J2_n(end,3),'rx','linewidth',2,'markersize',8)
p3 = plot3(L123_noJ2n(Lpoint,1),L123_noJ2n(Lpoint,2),L123_noJ2n(Lpoint,3),'^','markerfacecolor',colors.std.black,'markeredgecolor',colors.std.black);
plotBodyTexture3(secondary.R_n, [1-secondary.MR, 0, 0],secondary.img);
% plotBody3(secondary.R_n, [1-secondary.MR, 0, 0], secondary.color)
PlotBoi3('$x_n$','$y_n$','$z_n$',20,'LaTex')
legend([p1,p2,p3],'Without J2','With J_2',sprintf('L%1.0f, Without J2',Lpoint))
% legend([p1,p2,p3],'Without J_2','With J_2',sprintf('L-point, Without J_2',Lpoint))
axis equal
view(0,90)
grid on

end

% ========================================================================
%%% Comparison of J2 and J4 forces along arbitrary trajectory
% ========================================================================
if run_J2J4ForceComp == 1
primary = bodies.jupiter;
secondary = bodies.europa;
L123_noJ2n = EquilibriumPoints(secondary.MR,1:3); % [3x3] of L1, L2, L3 normalized BCR coordinates
r0_n = [1.009277654439744; 0; 0];
v0_n = [0;0;0.057238983158620];
X0_n = [r0_n; v0_n];
tf = 11.500597828347050;
Lpoint = 2;

%%% Normalizing constants
rNorm = secondary.a;         % n <-> km
tNorm = 1/secondary.meanMot; % n <-> sec
vNorm = rNorm / tNorm;       % n <-> km/sec

%%% Setting time vector and normalizing 
ti = 0; % sec
dt = tf/10000;
time0_n = [ti:dt:tf];

%%% Choosing ode45 tolerance
tol = 1e-10;

%%% Setting integrator options
options = odeset('RelTol',tol,'AbsTol',tol);

%%% Setting prms
prms.u = secondary.MR;
prms.R1_n = primary.R/rNorm;
prms.R2_n = secondary.R_n;
prms.J21 = primary.J2;
prms.J22 = 0;
prms.J31 = 0;
prms.J32 = 0;
prms.J41 = primary.J4;
prms.J42 = 0;

prmsNA.u = secondary.MR;
prmsNA.R1_n = primary.R/rNorm;
prmsNA.R2_n = secondary.R_n;
prmsNA.J21 = 0;
prmsNA.J22 = 0;
prmsNA.J31 = 0;
prmsNA.J32 = 0;
prmsNA.J41 = 0;
prmsNA.J42 = 0;

prmsJ2.u = secondary.MR;
prmsJ2.R1_n = primary.R/rNorm;
prmsJ2.R2_n = secondary.R_n;
prmsJ2.J21 = primary.J2;
prmsJ2.J22 = 0;
prmsJ2.J31 = 0;
prmsJ2.J32 = 0;
prmsJ2.J41 = 0;
prmsJ2.J42 = 0;

prmsJ4.u = secondary.MR;
prmsJ4.R1_n = primary.R/rNorm;
prmsJ4.R2_n = secondary.R_n;
prmsJ4.J21 = 0;
prmsJ4.J22 = 0;
prmsJ4.J31 = 0;
prmsJ4.J32 = 0;
prmsJ4.J41 = primary.J4;
prmsJ4.J42 = 0;

%%% Propagating the States
[time_n, X_BCR_J2_n] = ode45(@Int_CR3Bn_J2, time0_n, X0_n, options, prms);

[time_n, X_BCR_J2J4_n] = ode45(@Int_CR3Bn_ZH, time0_n, X0_n, options, prms);

aNA = zeros(size(X_BCR_J2_n,1),1);
aJ2 = zeros(size(X_BCR_J2_n,1),1);
aJ4 = zeros(size(X_BCR_J2_n,1),1);
for kk = 1:size(X_BCR_J2_n,1)
    acc = Int_CR3Bn_ZH(0,X_BCR_J2_n(kk,:),prmsNA);
    aNA(kk) = norm(acc(4:6));
    
    acc = Int_CR3Bn_ZH(0,X_BCR_J2_n(kk,:),prmsJ2);
    aJ2(kk) = norm(acc(4:6));
    
    acc = Int_CR3Bn_ZH(0,X_BCR_J2_n(kk,:),prmsJ4);
    aJ4(kk) = norm(acc(4:6));
    
%     clc
%     aNA(kk)
%     aJ2(kk)
%     aJ4(kk)
%     
%     aJ2(kk) - aNA(kk)
%     aJ4(kk) - aNA(kk)
%     (aJ2(kk) - aNA(kk)) / (aJ4(kk) - aNA(kk))
end


figure
subplot(3,2,[1,3,5]); hold all
plot3(X_BCR_J2_n(:,1),X_BCR_J2_n(:,2),X_BCR_J2_n(:,3),'r','linewidth',2);
plot3(X_BCR_J2J4_n(:,1),X_BCR_J2J4_n(:,2),X_BCR_J2J4_n(:,3),'b','linewidth',2);
plotBodyTexture3(secondary.R_n, [1-secondary.MR, 0, 0],secondary.img);
% plotBody3(secondary.R_n, [1-secondary.MR, 0, 0], secondary.color)
legend('J2','J2 & J4')
PlotBoi3('$x_n$','$y_n$','$z_n$',20,'LaTex')
axis equal
view(-18,16)
grid on
subplot(3,2,2)
plot(time_n,(X_BCR_J2J4_n(:,1)-X_BCR_J2_n(:,1)).*rNorm,'k','linewidth',2)
PlotBoi2('','$\Delta x$, km',20,'LaTex')
subplot(3,2,4)
plot(time_n,(X_BCR_J2J4_n(:,2)-X_BCR_J2_n(:,2)).*rNorm,'k','linewidth',2)
PlotBoi2('','$\Delta y$, km',20,'LaTex')
subplot(3,2,6)
plot(time_n,(X_BCR_J2J4_n(:,3)-X_BCR_J2_n(:,3)).*rNorm,'k','linewidth',2)
PlotBoi2('Normalized Time','$\Delta z$, km',20,'LaTex')

% percentDiffJ2J4 = (aJ2-aNA)./(aJ4-aNA);
% figure; hold all
% plot(time_n, percentDiffJ2J4,'.')
% PlotBoi2('Time_n','$F_{J2}$ / $F_{J4}$',20,'LaTex')


figure
subplot(2,2,[1,3]); hold all
plot3(X_BCR_J2_n(:,1),X_BCR_J2_n(:,2),X_BCR_J2_n(:,3),'r','linewidth',2);
plotBodyTexture3(secondary.R_n, [1-secondary.MR, 0, 0],secondary.img);
% plotBody3(secondary.R_n, [1-secondary.MR, 0, 0], secondary.color)
PlotBoi3('$x_n$','$y_n$','$z_n$',20,'LaTex')
axis equal
view(-18,16)
grid on
subplot(2,2,2)
plot(time_n, aJ2-aNA,'k','linewidth',2)
PlotBoi2('','$|a_{J2}|$',20,'LaTex')
subplot(2,2,4)
plot(time_n, aJ4-aNA,'k','linewidth',2)
PlotBoi2('Normalized time','$|a_{J4}|$',20,'LaTex')

end

% ========================================================================
%%% Comparison of J2 no J2 forces across 3B system
% ========================================================================
if run_J2ForceCompMesh == 1
%%% Setting Systems
primary = bodies.jupiter;
secondary = bodies.europa;

%%% Normalizing constants
rNorm = secondary.a;         % n <-> km
tNorm = 1/secondary.meanMot; % n <-> sec
vNorm = rNorm / tNorm;       % n <-> km/sec

% %%% Whole system
% xs = linspace(-.1,1.1,100);
% ys = linspace(-.2,.2,100);

%%% Secondary
xs = linspace(0.9,1.1,500);
% ys = linspace(-2*secondary.R_n,2*secondary.R_n,100);
ys = linspace(-.025,.025,500);

Z = 0;

[X, Y] = meshgrid(xs,ys);
clear xs ys

acc_J2 = zeros(size(X));
acc_NA = zeros(size(X));

prmsNA.u = secondary.MR;
prmsNA.R1_n = primary.R/rNorm;
prmsNA.R2_n = secondary.R_n;
prmsNA.J21 = 0;
prmsNA.J22 = 0;
prmsNA.J31 = 0;
prmsNA.J32 = 0;
prmsNA.J41 = 0;
prmsNA.J42 = 0;

prmsJ2.u = secondary.MR;
prmsJ2.R1_n = primary.R/rNorm;
prmsJ2.R2_n = secondary.R_n;
prmsJ2.J21 = primary.J2;
prmsJ2.J22 = 0;
prmsJ2.J31 = 0;
prmsJ2.J32 = 0;
prmsJ2.J41 = 0;
prmsJ2.J42 = 0;


for xk = 1:size(X,1)
    for yk = 1:size(X,2)

        %%% J2 Zero-Velocity Curve
        acc = Int_CR3Bn_ZH(0,[X(xk,yk),Y(xk,yk),Z,0,0,0],prmsJ2);
        acc_J2(xk,yk) = norm(acc(4:6)).*rNorm./(tNorm^2)*1000; % m/s^2

        %%% Zero-Velocity Curve
        acc = Int_CR3Bn_ZH(0,[X(xk,yk),Y(xk,yk),Z,0,0,0],prmsNA);
        acc_NA(xk,yk) = norm(acc(4:6)).*rNorm./(tNorm^2)*1000; % m/s^2
    end
end

acc_diff = acc_J2 - acc_NA;

figure; hold all
[C_J2,h_J2] = contour(X,Y,acc_diff,'color','k');
plotBody2( secondary.R_n, [1-secondary.MR,0,0], secondary.color ,[0,0,0],1 )
plotBody2( primary.R/rNorm, [-secondary.MR,0,0], primary.color ,[0,0,0],1 )
axis equal; clabel(C_J2,h_J2)
% xlim([0.93, 1.07]); ylim([-0.025, 0.025])
PlotBoi2('x_n - BCR','y_n - BCR',12)


figure; hold all
[C,h] = contourf(X,Y,acc_diff,linspace(min(min(acc_diff)),max(max(acc_diff)),100),'edgecolor','none');
colorbar
% plotBody2( secondary.R_n, [1-secondary.MR,0,0], secondary.color ,[0,0,0],1)
% plotBody2( primary.R/rNorm, [-secondary.MR,0,0], primary.color ,[0,0,0],1 )
axis equal
% xlim([0.93, 1.07]); ylim([-0.02, 0.02])
PlotBoi2('$x_n$','$y_n$',20,'LaTex')
plotBody2(secondary.R_n,[1-secondary.MR,0,0],colors.std.black,colors.std.black,2,0)

figure
surf(acc_diff)
zlim([0 5])

end
% ========================================================================
%%% Analysis of Movement of Linear Equilibrium Points for Different Systems
% ========================================================================
if run_EquillibriumMovementTable == 1
% -----------------------------
% Options
% -----------------------------
%%% Allow for J2 of secondary?
allowJ22 = 1;

% -----------------------------
% Setup
% -----------------------------
%%% Setting 3B systems to loop through
systems = {bodies.earth,bodies.moon;...
          bodies.mars,bodies.phobos;...
          bodies.mars,bodies.deimos;...
          bodies.jupiter,bodies.io;... 
          bodies.jupiter,bodies.europa;...
          bodies.jupiter,bodies.ganymede;...
          bodies.saturn,bodies.enceladus;...
          bodies.saturn,bodies.titan;...
          bodies.neptune,bodies.triton};

%%% Preallocating results
dataTable = {'PRIMARY','SECONDARY','MR','J2,2 CONSIDERED','L1 DIFF_n','L2 DIFF_n','L3 DIFF_n','L1 DIFF, km','L2 DIFF, km','L3 DIFF, km'};

% -----------------------------
% Running Analysis
% -----------------------------
%%% Looping through systems and building table of L123 differences
for kk = 1:size(systems,1)
    %%% Reassigning data for clarity
    primary   = systems{kk,1};
    secondary = systems{kk,2};
    
    %%% Assigning rNorm
    rNorm = secondary.a; % n <-> km
    
    %%% Acquire Lagrange points of both models
    L123_n = EquilibriumPoints(secondary.MR,1:3);
    if allowJ22 == 1 % If allowing J22
        if isfield(secondary,'J2') == 0 % If the body doesn't have a J2 in getBodyData(), set equal to 0
            secondary.J2 = 0;
        end
    else % If NOT allowing J22
        secondary.J2 = 0;
    end
    L123_J2n = EquilibriumPoints_J2(secondary.MR,primary.J2,secondary.J2,primary.R/rNorm,secondary.R_n);
    LsDiff = L123_J2n-L123_n;

    %%% Storing data
    dataTable{kk+1,1} = primary.name;
    dataTable{kk+1,2} = secondary.name;
    dataTable{kk+1,3} = secondary.MR;
    dataTable{kk+1,5} = LsDiff(1);
    dataTable{kk+1,6} = LsDiff(2);
    dataTable{kk+1,7} = LsDiff(3);
    dataTable{kk+1,8} = LsDiff(1)*rNorm;
    dataTable{kk+1,9} = LsDiff(2)*rNorm;
    dataTable{kk+1,10} = LsDiff(3)*rNorm;
    if isfield(secondary,'J2') == 0 || secondary.J2 == 0
        dataTable{kk+1,4} = 'No';
    else
        dataTable{kk+1,4} = 'Yes';
    end
end
dataTable
end

% ========================================================================
%%% Zero Velocity Curves
% ========================================================================
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % check out contourf() given an nxn matrix
if run_zeroVelocity == 1
%%% Setting Systems
primary = bodies.jupiter;
secondary = bodies.europa;

%%% Setting Normalizing factors
rNorm = secondary.a;

% % Whole system
% xs = linspace(-.1,1.1,1000);
% ys = linspace(-.2,.2,1000);

% Secondary
xs = linspace(0.9,1.1,600);
% ys = linspace(-2*secondary.R_n,2*secondary.R_n,100);
ys = linspace(-.025,.025,600);

[X, Y] = meshgrid(xs,ys);
clear xs ys

Z_J2 = zeros(size(X));
Z = zeros(size(X));
for xk = 1:size(X,1)
    for yk = 1:size(X,2)

        %%% J2 Zero-Velocity Curve
        zv_J2 = JacobiConstantCalculator_J2(secondary.MR,[X(xk,yk), Y(xk,yk), 0] ,[0, 0, 0] , primary.R/rNorm, secondary.R_n, primary.J2, 0);
        Z_J2(xk,yk) = zv_J2;

        %%% Zero-Velocity Curve
        zv = JacobiConstantCalculator(secondary.MR,[X(xk,yk), Y(xk,yk), 0] ,[0, 0, 0]);
        Z(xk,yk) = zv;
    end
end
% -----------------------------
% Getting contour levels
% -----------------------------
%%% Acquire Lagrange points
L123_noJ2n = EquilibriumPoints(secondary.MR,1:3); % [3x3] of L1, L2, L3 normalized BCR coordinates
L123_J2n = EquilibriumPoints_J2(secondary.MR,primary.J2,0,primary.R/rNorm,secondary.R_n);
L123 = [L123_noJ2n; L123_J2n]; % stacking

%%% No J2, Surface
[JC_Spx] = JacobiConstantCalculator(secondary.MR,[(1-secondary.MR)+secondary.R_n,0,0],[0,0,0]);
[JC_Smx] = JacobiConstantCalculator(secondary.MR,[(1-secondary.MR)-secondary.R_n,0,0],[0,0,0]);
[JC_Spy] = JacobiConstantCalculator(secondary.MR,[(1-secondary.MR),secondary.R_n,0],[0,0,0]);
[JC_Smy] = JacobiConstantCalculator(secondary.MR,[(1-secondary.MR),-secondary.R_n,0],[0,0,0]);
[JC_Spz] = JacobiConstantCalculator(secondary.MR,[(1-secondary.MR),0,secondary.R_n],[0,0,0]);
[JC_Smz] = JacobiConstantCalculator(secondary.MR,[(1-secondary.MR),0,-secondary.R_n],[0,0,0]);

%%% No J2, Lagrange
[JC_L1] = JacobiConstantCalculator(secondary.MR,L123_noJ2n(1,:),[0,0,0]);
[JC_L2] = JacobiConstantCalculator(secondary.MR,L123_noJ2n(2,:),[0,0,0]);
[JC_L3] = JacobiConstantCalculator(secondary.MR,L123_noJ2n(3,:),[0,0,0]);

%%% J2, Surface
[JC_SpxJ2] = JacobiConstantCalculator_J2(secondary.MR,[(1-secondary.MR)+secondary.R_n,0,0],[0,0,0], primary.R/rNorm, secondary.R_n, primary.J2, 0);
[JC_SmxJ2] = JacobiConstantCalculator_J2(secondary.MR,[(1-secondary.MR)-secondary.R_n,0,0],[0,0,0], primary.R/rNorm, secondary.R_n, primary.J2, 0);
[JC_SpyJ2] = JacobiConstantCalculator_J2(secondary.MR,[(1-secondary.MR),secondary.R_n,0],[0,0,0], primary.R/rNorm, secondary.R_n, primary.J2, 0);
[JC_SmyJ2] = JacobiConstantCalculator_J2(secondary.MR,[(1-secondary.MR),-secondary.R_n,0],[0,0,0], primary.R/rNorm, secondary.R_n, primary.J2, 0);
[JC_SpzJ2] = JacobiConstantCalculator_J2(secondary.MR,[(1-secondary.MR),0,secondary.R_n],[0,0,0], primary.R/rNorm, secondary.R_n, primary.J2, 0);
[JC_SmzJ2] = JacobiConstantCalculator_J2(secondary.MR,[(1-secondary.MR),0,-secondary.R_n],[0,0,0], primary.R/rNorm, secondary.R_n, primary.J2, 0);

%%% J2, Lagrange
[JC_L1J2] = JacobiConstantCalculator_J2(secondary.MR,L123_J2n(1,:),[0,0,0], primary.R/rNorm, secondary.R_n, primary.J2, 0);
[JC_L2J2] = JacobiConstantCalculator_J2(secondary.MR,L123_J2n(2,:),[0,0,0], primary.R/rNorm, secondary.R_n, primary.J2, 0);
[JC_L3J2] = JacobiConstantCalculator_J2(secondary.MR,L123_J2n(3,:),[0,0,0], primary.R/rNorm, secondary.R_n, primary.J2, 0);

% %%% General full system
% levels = [3.001, 3.00361, 3.0038, 3.005,3.006,3.007, 3.009, 3.4, 9, 15];

%%% Secondary plot
levels = [JC_L1, JC_L2, 3.001, 3.005,3.007, 3.009];

figure('position',[93 413 595 382]);
subplot(2,1,1); hold all; title('Without J_2')
[C,h] = contour(X,Y,Z,levels,'color','k');
plotBody2( secondary.R_n, [1-secondary.MR,0,0], secondary.color ,[0,0,0],1)
plotBody2( primary.R/rNorm, [-secondary.MR,0,0], primary.color ,[0,0,0],1 )
axis equal; clabel(C,h)
xlim([0.93, 1.07]); ylim([-0.025, 0.025])
PlotBoi2('','$y_n$',20,'Latex')

subplot(2,1,2); hold all; title('J2')
[C_J2,h_J2] = contour(X,Y,Z_J2,levels,'color','k');
plotBody2( secondary.R_n, [1-secondary.MR,0,0], secondary.color ,[0,0,0],1 )
plotBody2( primary.R/rNorm, [-secondary.MR,0,0], primary.color ,[0,0,0],1 )
axis equal; clabel(C_J2,h_J2)
xlim([0.93, 1.07]); ylim([-0.025, 0.025])
PlotBoi2('$x_n$','$y_n$',20,'Latex')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('position',[93 413 595 382]);
subplot(2,1,1); hold all
[C,h] = contourf(X,Y,Z,levels);
% plotBody2( secondary.R_n, [1-secondary.MR,0,0], secondary.color ,[0,0,0],1)
% plotBody2( primary.R/rNorm, [-secondary.MR,0,0], primary.color ,[0,0,0],1 )
axis equal
xlim([0.93, 1.07]); ylim([-0.02, 0.02])
PlotBoi2('','$y_n$',20,'Latex')
title('(a)')

subplot(2,1,2); hold all
[C_J2,h_J2] = contourf(X,Y,Z_J2,levels);
%plotBody2( secondary.R_n, [1-secondary.MR,0,0], secondary.color ,[0,0,0],1 )
%plotBody2( primary.R/rNorm, [-secondary.MR,0,0], primary.color ,[0,0,0],1 )
axis equal
xlim([0.93, 1.07]); ylim([-0.02, 0.02])
PlotBoi2('$x_n$','$y_n$',20,'Latex')
title('(b)')

cmap = parula(5);
% cmap = jet(5);
colormap(cmap);
colorbar('Position',[.92 .1 .03 .825])
% % hc = gcf;
% % hp = hc.Position;
% % colorbar('Position', [hp(1)+hp(3)+0.01  hp(2)  0.03  hp(2)+hp(3)*0.9])
% colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.03  hp4(2)+hp4(3)*0.9])

end

% ========================================================================
%%% Multiple Run Simulation
% ========================================================================
if run_multSimulation == 1
%%% Setting Systems
primary = bodies.jupiter;
secondary = bodies.europa;

%%% Setting Normalizing factors
rNorm = secondary.a;
tNorm = 1/secondary.meanMot; % n <-> sec
vNorm = rNorm / tNorm;       % n <-> km/sec

%%% Acquire Lagrange points
L123_noJ2n = EquilibriumPoints(secondary.MR,1:3); % [3x3] of L1, L2, L3 normalized BCR coordinates
L123_J2n = EquilibriumPoints_J2(secondary.MR,primary.J2,0,primary.R/rNorm,secondary.R_n);
L123 = [L123_noJ2n; L123_J2n]; % stacking

% %%% Setting time vector and normalizing 
ti_mult = 0; % sec
dt_mult = 60*10; % sec 
tf_mult = 3600*24; % sec 
time0_n_mult = [ti_mult:dt_mult:tf_mult] ./ tNorm;
    
% time0_n_mult = [0:0.004629904709320:6.667062781420788];
    
%%% Set ranges for initial positions
xrange_mult = linspace(0, 0, 1)./rNorm; % enter in km
yrange_mult = linspace(0, 0, 1)./rNorm; % enter in km
zrange_mult = linspace(0, 0, 1)./rNorm; % enter in km

%%% Set ranges for initial velocities
dxrange_mult = linspace(-10, 0, 5)./1000./vNorm; % enter in m/s
dyrange_mult = linspace(0, 0, 1)./1000./vNorm; % enter in m/s
dzrange_mult = linspace(-40, 40, 10)./1000./vNorm; % enter in m/s
    
% %%% Set ranges for initial positions
% xrange_mult = linspace(0, 0, 1); 
% yrange_mult = linspace(0, 0, 1); 
% zrange_mult = linspace(0, 0, 1); 
% 
% %%% Set ranges for initial velocities
% dxrange_mult = linspace(-7.918473792209050e-04, 0, 3);
% dyrange_mult = linspace(0, 0, 1);
% dzrange_mult = linspace(-0.003167389516884, 0.003167389516884, 3);

%%% Create ordered matrices of intial conditions
positionNudgeMatrix = [];
for xx = 1:length(xrange_mult)
    for yy = 1:length(yrange_mult)
        for zz = 1:length(zrange_mult)
            positionNudgeMatrix = [positionNudgeMatrix; xrange_mult(xx), yrange_mult(yy), zrange_mult(zz)];
        end
    end
end

velocityNudgeMatrix = [];
for xx = 1:length(dxrange_mult)
    for yy = 1:length(dyrange_mult)
        for zz = 1:length(dzrange_mult)
            velocityNudgeMatrix = [velocityNudgeMatrix; dxrange_mult(xx), dyrange_mult(yy), dzrange_mult(zz)];
        end
    end
end

%%% Preallocating
r_BCR_mult = [];
r_BCR_multJ2 = [];
endSpeeds_n_mult = [];
endFPAs_mult = [];
latlons_mult = [];
firstPeriapsis_mult = [];
anyLandings = 0;
periapsisFound = 0;

for pp = 1:size(positionNudgeMatrix,1)
    for vv = 1:size(velocityNudgeMatrix,1)

        %%% Adding nudges to the initial state
        x0_n_mult = L123(2,:) + positionNudgeMatrix(pp,:);
        xdot0_n_mult = [0, 0, 0] + velocityNudgeMatrix(vv,:);
        X0_n_mult = [x0_n_mult xdot0_n_mult]';

        %%% Choosing ode45 tolerance
        tol = 1e-10;

        %%% Setting integrator options
        if firstPeriapsisAnalysis_mult == 0 % Stop on impact
            options1 = odeset('Events',@event_Impact_CR3Bn,'RelTol',tol,'AbsTol',tol);
            options2 = odeset('Events',@event_Impact_CR3Bn_J2,'RelTol',tol,'AbsTol',tol);
        elseif firstPeriapsisAnalysis_mult == 1 % Don't Stop on Impact
            options1 = odeset('RelTol',tol,'AbsTol',tol);
            options2 = odeset('RelTol',tol,'AbsTol',tol);
        end
        
        %%% Setting prms
        prms.u = secondary.MR;
        prms.R1_n = primary.R/rNorm;
        prms.R2_n = secondary.R_n;
        prms.J21 = primary.J2;
        prms.J22 = 0;

        %%% Propagating the State with and without J2
        [time_n_mult, X_BCR_n_mult] = ode45(@Int_CR3Bn, time0_n_mult, X0_n_mult, options1, prms);
        [time_n_multJ2, X_BCR_n_multJ2] = ode45(@Int_CR3Bn_J2, time0_n_mult, X0_n_mult, options2, prms);
        
        %%% Storing new state positions
        r_BCR_mult = [r_BCR_mult; X_BCR_n_mult(:,1:3); nan(1,3)];
        r_BCR_multJ2 = [r_BCR_multJ2; X_BCR_n_multJ2(:,1:3); nan(1,3)];

        %%% Creating SCR positions
        r_SCR_n_mult = X_BCR_n_mult(:,1:3) - [1-secondary.MR, 0, 0];

        %%% Analyzing Impacts
        if time_n_mult(end) ~= time0_n_mult(end) % Crash
            anyLandings = 1;

            %%% Calculate end speed and flight path angle
            [ endFPA ] = rv2fpa( X_BCR_n_mult(end,1:3), X_BCR_n_mult(end,4:6)); % rad
            endSpeed_n = norm(X_BCR_n_mult(end,4:6));

            %%% Storing new end speed and flight path angles
            endSpeeds_n_mult = [endSpeeds_n_mult; endSpeed_n];
            endFPAs_mult = [endFPAs_mult; endFPA];

            %%% Acquire Lat/Lon of landing site
            angleUnit = 'degrees';
            [lat, lon] = ECEF2latlon(r_SCR_n_mult(end,1:3),angleUnit);

            %%% Storing Lat/Lon of landing site
            latlons_mult = [latlons_mult; lat, lon];
        end

        %%% Periapsis Analysis
        if firstPeriapsisAnalysis_mult == 1
            secondaryDistances = rowNorm(r_SCR_n_mult);
            for kk = 2:size(secondaryDistances,1)
                if secondaryDistances(kk-1) < secondaryDistances(kk) && norm(secondaryDistances(kk) - secondaryDistances(kk-1)) > tol
                    periapsisFound = 1;
                    firstPeriapsis_mult = [firstPeriapsis_mult; X_BCR_n_mult(kk,1:3)];
                    break
                end
            end
            if periapsisFound == 0 && secondaryDistances(1) > secondaryDistances(end)
                firstPeriapsis_mult = [firstPeriapsis_mult; X_BCR_n_mult(end,1:3)];
            elseif periapsisFound == 0 && secondaryDistances(1) < secondaryDistances(end)
                firstPeriapsis_mult = [firstPeriapsis_mult; X_BCR_n_mult(1,1:3)];
            end
        end
    end
end

% Clearing iterative variables
clear lat lon endFPA endSpeed_n X_BCR_n_mult

%%% Plotting Results of Multiple Runs
if plotBCR_n_mult == 1
    % BCR Trajectories
    figure; hold all
    p1 = plot3(r_BCR_mult(:,1),r_BCR_mult(:,2),r_BCR_mult(:,3),'color',colors.sch.d4_1(2,:),'linewidth',1.5);
    p2 = plot3(r_BCR_multJ2(:,1),r_BCR_multJ2(:,2),r_BCR_multJ2(:,3),'color',colors.sch.d4_1(3,:),'linewidth',1.5);
    plot3(L123(2,1),L123(2,2),L123(2,3),'k^','markersize',5,'linewidth',1.5)
    plot3(L123(1,1),L123(1,2),L123(1,3),'k^','markersize',5,'linewidth',1.5)
%     plotBody3(secondary.R_n, [1-secondary.MR, 0, 0], secondary.color)
    plotBodyTexture3(secondary.R_n, [1-secondary.MR, 0, 0],secondary.img);
    PlotBoi3('X_n','Y_n','Z_n',20)
    axis equal
    xlim([1-secondary.MR-2*secondary.R_n, L123(2,1) + 2*secondary.R_n])
    ylim([-5*secondary.R_n, 5*secondary.R_n])
%     xlim([L123(1,1)-2*secondary.R_n, L123(2,1)+2*secondary.R_n])
%     ylim([-12*secondary.R_n, 12*secondary.R_n])
%     zlim([-5*secondary.R_n, 5*secondary.R_n])
    legend([p1, p2],'Without J_2','With J_2')
    view(180,0)
    if firstPeriapsisAnalysis_mult == 1
        alpha(.6)
    end
end

if anyLandings == 1
    if fpaVSspeed_mult == 1
        % Terminal FPA vs Speed
        figure; hold all
        plot(endFPAs_mult*180/pi,endSpeeds_n_mult, 'bo','markersize',10,'linewidth',1.5)
        PlotBoi2('Flight Path Angle at Crash, °','Normalized Speed at Crash',12)
    end

    if latlonPlot_mult == 1
        % Terminal Lat/Lon
        figure; hold all
        plot(latlons_mult(:,2), latlons_mult(:,1), 'rx','markersize',5,'linewidth',1.5);
        PlotBoi2('Longitude, °', 'Latitude, °',20)
        xlim([-180 180])
        ylim([-90 90])
        grid on
    end
end

if firstPeriapsisAnalysis_mult == 1
    figure
    subplot(1,2,1); hold all
    plot3(firstPeriapsis_mult(:,1),firstPeriapsis_mult(:,2),firstPeriapsis_mult(:,3),'rx','markersize',5,'linewidth',1)
    PlotBoi3('X_n','Y_n','Z_n',20)
    axis equal
    xlim([1-secondary.MR-5*secondary.R_n, 1-secondary.MR+5*secondary.R_n])
    ylim([-5*secondary.R_n, 5*secondary.R_n])
    zlim([-5*secondary.R_n, 5*secondary.R_n])

    subplot(1,2,2); hold all
    plot3(firstPeriapsis_mult(:,1),firstPeriapsis_mult(:,2),firstPeriapsis_mult(:,3),'rx','markersize',5,'linewidth',1)
    plot3(r_BCR_mult(:,1),r_BCR_mult(:,2),r_BCR_mult(:,3),'g','linewidth',1.0)
    plotBody3(secondary.R_n, [1-secondary.MR, 0, 0], secondary.color)
    alpha(.3)
    PlotBoi3('X_n','Y_n','Z_n',20)
    axis equal
    xlim([1-secondary.MR-5*secondary.R_n, 1-secondary.MR+5*secondary.R_n])
    ylim([-5*secondary.R_n, 5*secondary.R_n])
    zlim([-5*secondary.R_n, 5*secondary.R_n])
end
    
end

% ========================================================================
%%% Surface Jacobi Constant Analysis
% ========================================================================
% ***Note for JC motion: S/C CAN reach anywhere with a HIGHER JC ... CANNOT
% reach anywhere with a LOWER JC***
if run_SurfaceJCAnalysis == 1
%%% Selecting bodies
primary = bodies.jupiter;
secondary = bodies.io;

%%% Normalizing constants
rNorm = secondary.a;         % n <-> km
tNorm = 1/secondary.meanMot; % n <-> sec
vNorm = rNorm / tNorm;       % n <-> km/sec

%%% Acquire Lagrange points
L123_noJ2n = EquilibriumPoints(secondary.MR,1:3); % [3x3] of L1, L2, L3 normalized BCR coordinates
L123_J2n = EquilibriumPoints_J2(secondary.MR,primary.J2,0,primary.R/rNorm,secondary.R_n);
L123 = [L123_noJ2n; L123_J2n]; % stacking

%%% No J2, Surface
[JC_Spx] = JacobiConstantCalculator(secondary.MR,[(1-secondary.MR)+secondary.R_n,0,0],[0,0,0]);
[JC_Smx] = JacobiConstantCalculator(secondary.MR,[(1-secondary.MR)-secondary.R_n,0,0],[0,0,0]);
[JC_Spy] = JacobiConstantCalculator(secondary.MR,[(1-secondary.MR),secondary.R_n,0],[0,0,0]);
[JC_Smy] = JacobiConstantCalculator(secondary.MR,[(1-secondary.MR),-secondary.R_n,0],[0,0,0]);
[JC_Spz] = JacobiConstantCalculator(secondary.MR,[(1-secondary.MR),0,secondary.R_n],[0,0,0]);
[JC_Smz] = JacobiConstantCalculator(secondary.MR,[(1-secondary.MR),0,-secondary.R_n],[0,0,0]);

%%% No J2, Lagrange
[JC_L1] = JacobiConstantCalculator(secondary.MR,L123_noJ2n(1,:),[0,0,0]);
[JC_L2] = JacobiConstantCalculator(secondary.MR,L123_noJ2n(2,:),[0,0,0]);
[JC_L3] = JacobiConstantCalculator(secondary.MR,L123_noJ2n(3,:),[0,0,0]);

%%% J2, Surface
[JC_SpxJ2] = JacobiConstantCalculator_J2(secondary.MR,[(1-secondary.MR)+secondary.R_n,0,0],[0,0,0], primary.R/rNorm, secondary.R_n, primary.J2, 0);
[JC_SmxJ2] = JacobiConstantCalculator_J2(secondary.MR,[(1-secondary.MR)-secondary.R_n,0,0],[0,0,0], primary.R/rNorm, secondary.R_n, primary.J2, 0);
[JC_SpyJ2] = JacobiConstantCalculator_J2(secondary.MR,[(1-secondary.MR),secondary.R_n,0],[0,0,0], primary.R/rNorm, secondary.R_n, primary.J2, 0);
[JC_SmyJ2] = JacobiConstantCalculator_J2(secondary.MR,[(1-secondary.MR),-secondary.R_n,0],[0,0,0], primary.R/rNorm, secondary.R_n, primary.J2, 0);
[JC_SpzJ2] = JacobiConstantCalculator_J2(secondary.MR,[(1-secondary.MR),0,secondary.R_n],[0,0,0], primary.R/rNorm, secondary.R_n, primary.J2, 0);
[JC_SmzJ2] = JacobiConstantCalculator_J2(secondary.MR,[(1-secondary.MR),0,-secondary.R_n],[0,0,0], primary.R/rNorm, secondary.R_n, primary.J2, 0);

%%% J2, Lagrange
[JC_L1J2] = JacobiConstantCalculator_J2(secondary.MR,L123_J2n(1,:),[0,0,0], primary.R/rNorm, secondary.R_n, primary.J2, 0);
[JC_L2J2] = JacobiConstantCalculator_J2(secondary.MR,L123_J2n(2,:),[0,0,0], primary.R/rNorm, secondary.R_n, primary.J2, 0);
[JC_L3J2] = JacobiConstantCalculator_J2(secondary.MR,L123_J2n(3,:),[0,0,0], primary.R/rNorm, secondary.R_n, primary.J2, 0);

figure
subplot(2,1,1)
plot([JC_Spx, JC_Smx,JC_Spy, JC_Smy,JC_Spz, JC_Smz])
subplot(2,1,2)
plot([JC_SpxJ2, JC_SmxJ2,JC_SpyJ2, JC_SmyJ2,JC_SpzJ2, JC_SmzJ2])

%%% Theoretical dV to get from L1 to L2
[dv_L1L2] = dJC_2_dv(JC_L1-JC_L2);
dv_m_per_sec = (dv_L1L2*vNorm)*1000
[dv_L1L2J2] = dJC_2_dv(JC_L1J2-JC_L2J2);
dvJ2_m_per_sec = (dv_L1L2J2*vNorm)*1000
delta_m_per_sec = (dvJ2_m_per_sec-dv_m_per_sec)


%%% Looking at theoretical dv to get from L1/L2 to stationary on surface
% [dv] = dJC_2_dv(abs(JC_L1-JC_Smx));
% dv_L1 = dv*vNorm*1000
% [dvJ2] = dJC_2_dv(abs(JC_L1J2-JC_SmxJ2));
% dvJ2_L1 = dvJ2*vNorm*1000
% 
% dvDiff_L1 = dvJ2_L1-dv_L1
% 
% 989898989898
% 
% [dv] = dJC_2_dv(abs(JC_L2-JC_Spx));
% dv_L2 = dv*vNorm*1000
% [dvJ2] = dJC_2_dv(abs(JC_L2J2-JC_SpxJ2));
% dvJ2_L2 = dvJ2*vNorm*1000
% 
% dvDiff_L2 = dvJ2_L2-dv_L2


end % run_SurfaceJCAnalysis
