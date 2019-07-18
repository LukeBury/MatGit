% ========================================================================
%%% Description
% ========================================================================
% 

% Created: 
% Author : Luke Bury, luke.bury@colorado.edu
% ========================================================================
%%% Initialization
% ========================================================================
clear
clc
close all
mbinPath = '~/CU_Google_Drive/Documents/MatGit/mbin';
addpath(genpath(mbinPath))
ticWhole = tic;

% ========================================================================
%%% Importing Data
% ========================================================================
%%% General data on solar system bodies
bodies = getBodyData(mbinPath);

%%% Color options/schemes
colors = get_colors();

%%% Periodic orbit ICs
PO_ICs = get_PO_ICs();

% ========================================================================
%%% Run Switches
% ========================================================================
storeAndPlotFullTrajectories = 1;

% ========================================================================
%%% Setup
% ========================================================================
% -------------------------------------------------
%%% System
% -------------------------------------------------
%%% Primary
uranus   = bodies.uranus;

%%% Secondaries
ophelia  = bodies.ophelia;
cordelia = bodies.cordelia;

%%% Normalizing constants
rNorm_ophelia  = ophelia.a;         % n <-> km
tNorm_ophelia  = ophelia.T_sec / 2*pi; % n <-> sec
vNorm_ophelia  = rNorm_ophelia / tNorm_ophelia; % n <-> km/sec

rNorm_cordelia = cordelia.a;         % n <-> km
tNorm_cordelia  = cordelia.T_sec / 2*pi; % n <-> sec
vNorm_cordelia  = rNorm_cordelia / tNorm_cordelia; % n <-> km/sec

%%% Equillibrium Points
rLPs_ophelia_n  = EquilibriumPoints(ophelia.MR);
rLPs_cordelia_n = EquilibriumPoints(cordelia.MR);

% ------------------------------------------------- 
%%% Integration Options
% -------------------------------------------------
tol = 1e-13;
options = odeset('RelTol',tol,'AbsTol',tol);
prmsOphelia.u  = ophelia.MR;
prmsCordelia.u = cordelia.MR;

% ========================================================================
%%% Perturb Lagrange points and compute manifolds
% ========================================================================
% -------------------------------------------------
%%% Perturb the lagrange points
% -------------------------------------------------
%%% Set perturbation vector
pert_n = [1,1].*1e-7;

%%% Add perturbation vector to L1 of Ophelia and L2 of Cordelia
X0n_ophelia_m  = [rLPs_ophelia_n(1,1:2) - pert_n, zeros(1,4)]';
X0n_cordelia_p = [rLPs_cordelia_n(2,1:2) + pert_n, zeros(1,4)]';

% ========================================================================
%%% Create initial conditions for ring particles
% ========================================================================
% Mean distance: 51149.3 km
% velocity: 10.650 km/s
% Period: 8.3823 hr
% width: 20 to 96 km
% e: 0.00794
% i: 0
a_epsilonRing_km = 51149.3; % km
minWidth_km = 26; % km
maxWidth_km = 96; % km
eccentricity_epsilonRing = 0.00794; % km

rp_epsilonRing_km = a_epsilonRing_km * (1 - eccentricity_epsilonRing);
ra_epsilonRing_km = a_epsilonRing_km * (1 + eccentricity_epsilonRing);

lowerBound_km = rp_epsilonRing_km - maxWidth_km;
upperBound_km = ra_epsilonRing_km + maxWidth_km;

n_particles = 5;
trueAnomalies_epsilonRing_deg = linspace(0,360 - 360/n_particles, n_particles);

X0s_inertial_epsilonRings = NaN(n_particles,6);
for kk = 1:n_particles
    [r, v] = COE2RV(a_epsilonRing_km, eccentricity_epsilonRing, 0, 0, 0, trueAnomalies_epsilonRing_deg(kk), uranus.u);
    X0s_inertial_epsilonRings(kk,:) = [r;v]';
end

% ========================================================================
%%% Prepare Epsilon Ring Particles for Integration
% ========================================================================
X0ns_rotating_epsilonRings_cordeliaFrame = NaN(size(X0s_inertial_epsilonRings));
X0ns_rotating_epsilonRings_opheliaFrame  = NaN(size(X0s_inertial_epsilonRings));

X0ns_rotating_epsilonRings_cordeliaFrame(:,1:3) = X0s_inertial_epsilonRings(:,1:3)./rNorm_cordelia;
X0ns_rotating_epsilonRings_opheliaFrame(:,1:3)  = X0s_inertial_epsilonRings(:,1:3)./rNorm_ophelia;

X0ns_rotating_epsilonRings_cordeliaFrame(:,4:6) = v_BodyCI2BodyCR(X0s_inertial_epsilonRings(:,4:6), X0s_inertial_epsilonRings(:,1:3), zeros(n_particles,1), cordelia.meanMot )./vNorm_cordelia;
X0ns_rotating_epsilonRings_opheliaFrame(:,4:6)  = v_BodyCI2BodyCR(X0s_inertial_epsilonRings(:,4:6), X0s_inertial_epsilonRings(:,1:3), zeros(n_particles,1), ophelia.meanMot )./vNorm_ophelia;


JC_ring_cordelia = JacobiConstantCalculator(cordelia.MR,X0ns_rotating_epsilonRings_cordeliaFrame(1,1:3),X0ns_rotating_epsilonRings_cordeliaFrame(1,4:6));
JC_ring_ophelia  = JacobiConstantCalculator(ophelia.MR,X0ns_rotating_epsilonRings_opheliaFrame(1,1:3),X0ns_rotating_epsilonRings_opheliaFrame(1,4:6));
JC_cordelia_L2   = JacobiConstantCalculator(cordelia.MR,rLPs_cordelia_n(2,:),[0,0,0]);
JC_ophelia_L1    = JacobiConstantCalculator(ophelia.MR,rLPs_ophelia_n(2,:),[0,0,0]);

% 989
% [Ttest, Xtest] = ode113(@Int_2BI, [0 86400*500], X0s_inertial_epsilonRings(3,:)', options, uranus.u);
% [testCordelia] = X_BCI2BCR(Xtest, Ttest, cordelia.meanMot );
% [testOphelia] = X_BCI2BCR(Xtest, Ttest, ophelia.meanMot );
% ========================================================================
%%% Integration
% ========================================================================
% -------------------------------------------------
%%% Integrate the perturbed states forwards and backworks in time
% -------------------------------------------------
%%% Create time vectors
t0 = 0;
% tf = 5000000*pi;
% tf = 500000*pi;
% tf = 50000*pi;
% tf = 5000*pi;
% tf = 500*pi;
tf = 10*pi;



time0_fwd = [t0,tf];
time0_bkd = [tf,t0];

X0s_parfor    = [X0n_ophelia_m'; X0n_ophelia_m'; X0n_cordelia_p'; X0n_cordelia_p'; X0ns_rotating_epsilonRings_cordeliaFrame; X0ns_rotating_epsilonRings_opheliaFrame];
time0s_parfor = [time0_bkd; time0_fwd; time0_bkd; time0_fwd; repmat(time0_fwd, n_particles*2,1)];
% prms_parfor   = [prmsOphelia; prmsOphelia; prmsCordelia; prmsCordelia; ];
rNorm_parfor = [rNorm_ophelia; rNorm_ophelia; rNorm_cordelia; rNorm_cordelia; repmat(rNorm_cordelia,n_particles,1); repmat(rNorm_ophelia,n_particles,1)];
vNorm_parfor = [vNorm_ophelia; vNorm_ophelia; vNorm_cordelia; vNorm_cordelia; repmat(vNorm_cordelia,n_particles,1); repmat(vNorm_ophelia,n_particles,1)];
tNorm_parfor = [tNorm_ophelia; tNorm_ophelia; tNorm_cordelia; tNorm_cordelia; repmat(tNorm_cordelia,n_particles,1); repmat(tNorm_ophelia,n_particles,1)];

options_yEqualsZero = odeset('Events',@event_yEqualsZero,'RelTol',tol,'AbsTol',tol);

trajs = {};
n_parfor = size(X0s_parfor,1);

parfor traj_i = 1:n_parfor
    
    
    X0_i = X0s_parfor(traj_i,:)';
    time0_i = time0s_parfor(traj_i,:);
%     prm_i = prms_parfor(traj_i);
    rNorm_i = rNorm_parfor(traj_i);
    vNorm_i = vNorm_parfor(traj_i);
    tNorm_i = tNorm_parfor(traj_i);
    
    prm_i= [];
    if rNorm_i == rNorm_ophelia
        prm_i = prmsOphelia;
    elseif rNorm_i == rNorm_cordelia
        prm_i = prmsCordelia;
    end
    
    [times_i, Xn_i, time_i_event, X_i_event, ~] = ode113(@Int_CR3Bn, time0_i, X0_i, options_yEqualsZero, prm_i);
    
    
    
    X_events_keep = [];
    t_events_keep = [];
    
    for kk = 1:size(X_i_event,1)
        if X_i_event(kk,1) < 0
            X_events_keep = [X_events_keep; X_i_event(kk,:)];
            t_events_keep = [t_events_keep; time_i_event(kk)];
        end
    end
    
    if storeAndPlotFullTrajectories == 1
        trajs{traj_i}.traj = [Xn_i(:,1:3).*rNorm_i, Xn_i(:,4:6).*vNorm_i];
        trajs{traj_i}.times = times_i.*tNorm_i;
    else
        trajs{traj_i}.traj = NaN;
        trajs{traj_i}.times = NaN;
    end
    trajs{traj_i}.X0   = [X0_i(1:3).*rNorm_i; X0_i(4:6).*vNorm_i];
    if isempty(X_events_keep) == 0
        trajs{traj_i}.X_events    = [X_events_keep(:,1:3).*rNorm_i, X_events_keep(:,4:6).*vNorm_i];
        trajs{traj_i}.time_events = t_events_keep.*tNorm_i;
    else
        trajs{traj_i}.X_events    = NaN(1,6);
        trajs{traj_i}.time_events = NaN;
    end
end

% -------------------------------------------------
%%% Seperate data for simplicity
% -------------------------------------------------
manifold_ophelia_m_stable    = trajs{1};
manifold_ophelia_m_unstable  = trajs{2};
manifold_cordelia_p_stable   = trajs{3};
manifold_cordelia_p_unstable = trajs{4};


cordeliaParticleEvents = [];
cordeliaParticleTraj   = [];
cordeliaParticleTimes  = [];
opheliaParticleEvents  = [];
opheliaParticleTraj    = [];
opheliaParticleTimes   = [];
% for kk = 5:length(trajs)
for kk = [5,10]
    if floor((kk-4-1)/n_particles) == 0
        %%% Cordelia set
        cordeliaParticleEvents = [cordeliaParticleEvents; trajs{kk}.X_events];
        cordeliaParticleTraj   = [cordeliaParticleTraj; trajs{kk}.traj; NaN(1,6)];
        cordeliaParticleTimes  = [cordeliaParticleTimes; trajs{kk}.times; NaN];
        
    elseif floor((kk-4-1)/n_particles) == 1
        %%% Ophelia set
        opheliaParticleEvents = [opheliaParticleEvents; trajs{kk}.X_events];
        opheliaParticleTraj   = [opheliaParticleTraj; trajs{kk}.traj; NaN(1,6)];
        opheliaParticleTimes  = [opheliaParticleTimes; trajs{kk}.times; NaN];
        
    else
        warning('algorithm error')
    end
end

[X_BCI_particle_cordelia] = X_BCR2BCI(cordeliaParticleTraj(1:end-1,:), cordeliaParticleTimes(1:end-1), cordelia.meanMot );
[X_BCI_particle_ophelia] = X_BCR2BCI(opheliaParticleTraj(1:end-1,:), opheliaParticleTimes(1:end-1), ophelia.meanMot );

X_ophelia_m_stable    = manifold_ophelia_m_stable.traj;
X_ophelia_m_unstable  = manifold_ophelia_m_unstable.traj;
X_cordelia_p_stable   = manifold_cordelia_p_stable.traj;
X_cordelia_p_unstable = manifold_cordelia_p_unstable.traj;

events_ophelia_m_stable    = manifold_ophelia_m_stable.X_events;
events_ophelia_m_unstable  = manifold_ophelia_m_unstable.X_events;
events_cordelia_p_stable   = manifold_cordelia_p_stable.X_events;
events_cordelia_p_unstable = manifold_cordelia_p_unstable.X_events;

if storeAndPlotFullTrajectories == 1
    %%% Ophelia
    r_ophelia_m_stable    = X_ophelia_m_stable(:,1:3);
    r_ophelia_m_unstable  = X_ophelia_m_unstable(:,1:3);

    %%% Cordelia
    r_cordelia_p_stable   = X_cordelia_p_stable(:,1:3);
    r_cordelia_p_unstable = X_cordelia_p_unstable(:,1:3);
end


% ========================================================================
%%% Create Fake Epsilon Ring
% ========================================================================
nPoints = 10000;
innerEdge = [ones(nPoints,1), zeros(nPoints,2)].*lowerBound_km;
outerEdge = [ones(nPoints,1), zeros(nPoints,2)].*upperBound_km;

thetas = linspace(0,2*pi,nPoints);
for kk = 2:nPoints
    innerEdge(kk,:) = R3([innerEdge(1,:)],thetas(kk));
    outerEdge(kk,:) = R3([outerEdge(1,:)],thetas(kk));
end

% ========================================================================
%%% Plot manifolds and epsilon ring
% ========================================================================
if storeAndPlotFullTrajectories == 1
lw = 2;


figure('position',[289 344 1061 392])
subplot(1,2,1); hold all

plot((1-ophelia.MR)*rNorm_ophelia,0,'o','markeredgecolor',colors.std.red,'markerfacecolor',colors.std.ltred,'markersize',8)
plot(rLPs_ophelia_n(1:2,1).*rNorm_ophelia,rLPs_ophelia_n(1:2,2).*rNorm_ophelia,'^','markeredgecolor',colors.std.red,'markerfacecolor',colors.std.ltred,'markersize',8)

plot((1-cordelia.MR)*rNorm_cordelia,0,'o','markeredgecolor',colors.std.blue,'markerfacecolor',colors.std.ltblue,'markersize',8)
plot(rLPs_cordelia_n(1:2,1).*rNorm_cordelia,rLPs_cordelia_n(1:2,2).*rNorm_cordelia,'^','markeredgecolor',colors.std.blue,'markerfacecolor',colors.std.ltblue,'markersize',8)

plot(r_ophelia_m_unstable(:,1), r_ophelia_m_unstable(:,2), 'r', 'linewidth', lw)
plot(r_ophelia_m_stable(:,1),   r_ophelia_m_stable(:,2),   'g', 'linewidth', lw)

plot(r_cordelia_p_unstable(:,1),r_cordelia_p_unstable(:,2),'r', 'linewidth', lw)
plot(r_cordelia_p_stable(:,1),  r_cordelia_p_stable(:,2),  'g', 'linewidth', lw)

%%% Plot Ring Coverage
% % fill([innerEdge(:,1); outerEdge(:,1)], [innerEdge(:,2); outerEdge(:,2)],colors.std.black)

PlotBoi2('$x$, $km$','$y$, $km$',20,'LaTex')
axis normal
axis equal


subplot(1,2,2); hold all

plot((1-ophelia.MR)*rNorm_ophelia,0,'o','markeredgecolor',colors.std.red,'markerfacecolor',colors.std.ltred,'markersize',8)
plot(rLPs_ophelia_n(1:2,1).*rNorm_ophelia,rLPs_ophelia_n(1:2,2).*rNorm_ophelia,'^','markeredgecolor',colors.std.red,'markerfacecolor',colors.std.ltred,'markersize',8)

plot((1-cordelia.MR)*rNorm_cordelia,0,'o','markeredgecolor',colors.std.blue,'markerfacecolor',colors.std.ltblue,'markersize',8)
plot(rLPs_cordelia_n(1:2,1).*rNorm_cordelia,rLPs_cordelia_n(1:2,2).*rNorm_cordelia,'^','markeredgecolor',colors.std.blue,'markerfacecolor',colors.std.ltblue,'markersize',8)


plot(r_ophelia_m_unstable(:,1), r_ophelia_m_unstable(:,2), 'r', 'linewidth', lw)
plot(r_ophelia_m_stable(:,1),   r_ophelia_m_stable(:,2),   'g', 'linewidth', lw)

p2 = plot(r_cordelia_p_unstable(:,1),r_cordelia_p_unstable(:,2),'r', 'linewidth', lw);
p1 = plot(r_cordelia_p_stable(:,1),  r_cordelia_p_stable(:,2),  'g', 'linewidth', lw);


%%% Plot Ring Coverage
% % ringPlot = fill([innerEdge(:,1); outerEdge(:,1)], [innerEdge(:,2); outerEdge(:,2)],colors.std.black);

PlotBoi2('$x$, $km$','$y$, $km$',20,'LaTex')
axis equal
% % legend([p1 p2 ringPlot],'Stable','Unstable','Epsilon Ring Range')

xlim([4.9 5.5].*1e4)
ylim([-1 1].*2000)

end % storeAndPlotFullTrajectories


% ========================================================================
%%% Plot Poincare Maps
% ========================================================================

minY_xdot = min([events_ophelia_m_stable(:,4); events_ophelia_m_unstable(:,4); events_cordelia_p_stable(:,4); events_cordelia_p_unstable(:,4)]);
maxY_xdot = max([events_ophelia_m_stable(:,4); events_ophelia_m_unstable(:,4); events_cordelia_p_stable(:,4); events_cordelia_p_unstable(:,4)]);
minY_ydot = min([events_ophelia_m_stable(:,5); events_ophelia_m_unstable(:,5); events_cordelia_p_stable(:,5); events_cordelia_p_unstable(:,5)]);
maxY_ydot = max([events_ophelia_m_stable(:,5); events_ophelia_m_unstable(:,5); events_cordelia_p_stable(:,5); events_cordelia_p_unstable(:,5)]);

fill_x = -[lowerBound_km lowerBound_km upperBound_km upperBound_km];
fill_y_xdot = [minY_xdot maxY_xdot maxY_xdot minY_xdot].*1.05;
fill_y_ydot = [minY_ydot maxY_ydot maxY_ydot minY_ydot].*1.05;

figure('position',[718 32 555 766])
subplot(2,1,1); hold all
plot(events_ophelia_m_stable(:,1),events_ophelia_m_stable(:,4),'g.')
plot(events_ophelia_m_unstable(:,1),events_ophelia_m_unstable(:,4),'r.')
plot(events_cordelia_p_stable(:,1),events_cordelia_p_stable(:,4),'g.')
plot(events_cordelia_p_unstable(:,1),events_cordelia_p_unstable(:,4),'r.')
fill(fill_x, fill_y_xdot, colors.std.black)
ylim([minY_xdot, maxY_xdot].*1.05)
PlotBoi2('$x$','$\dot{x}$',20,'LaTex')

subplot(2,1,2); hold all
plot(events_ophelia_m_stable(:,1),events_ophelia_m_stable(:,5),'g.')
plot(events_ophelia_m_unstable(:,1),events_ophelia_m_unstable(:,5),'r.')
plot(events_cordelia_p_stable(:,1),events_cordelia_p_stable(:,5),'g.')
plot(events_cordelia_p_unstable(:,1),events_cordelia_p_unstable(:,5),'r.')
fill(fill_x, fill_y_ydot, colors.std.black)
ylim([minY_ydot, maxY_ydot].*1.05)
PlotBoi2('$x$','$\dot{y}$',20,'LaTex')




% ========================================================================
%%% Formatting Structures
% ========================================================================
% -------------------------------------------------
%%% 
% -------------------------------------------------

% --------------------------
% 
% --------------------------

% ========================================================================
%%% Closeout
% ========================================================================
tocWhole = toc(ticWhole);
fprintf('Elapsed time: %1.4f seconds\n',tocWhole)
















