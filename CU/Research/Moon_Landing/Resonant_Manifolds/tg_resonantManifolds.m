% ========================================================================
%%% Description
% ========================================================================
% 

% Created: 07/12/19
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



% ========================================================================
%%% Setup
% ========================================================================
% -------------------------------------------------
%%% System
% -------------------------------------------------
primary   = bodies.jupiter;
secondary = bodies.europa;

%%% Normalizing constants
rNorm = secondary.a;         % n <-> km
tNorm = 1/secondary.meanMot; % n <-> sec
vNorm = rNorm / tNorm;       % n <-> km/sec

%%% Shortcut variables
mu   = secondary.MR;
R2_n = secondary.R_n;

%%% Equillibrium Points
rLPs_n = EquilibriumPoints(mu);

% ------------------------------------------------- 
%%% Integration Options
% -------------------------------------------------
tol = 5e-13;
options = odeset('RelTol',tol,'AbsTol',tol);
options_NegXAxisCrossing = odeset('Events',@event_negativeXAxisRetrogradeCrossing,'RelTol',tol,'AbsTol',tol);
prms.u    = mu;
prms.R2_n = R2_n;

% ========================================================================
%%% Simulation
% ========================================================================
% ------------------------------------------------- 
%%% Create initial conditions
% -------------------------------------------------
%%% Ranges for x0 and yDot0
nx0s    = 20;
nyDot0s = 20;
nTotal = nx0s * nyDot0s;
% x0_range    = linspace(-0.15,0,nx0s) - (1-secondary.MR);
x0_range    = linspace(-0.15,0,nx0s) - (1-secondary.MR);
% yDot0_range = linspace(0,0.07,nyDot0s);
yDot0_range = linspace(0,0.01,nyDot0s);

%%% Put all ICs into one matrix
X0s_mat = zeros(nTotal,6);
counter = 0;
for kk = 1:nx0s
    for jj = 1:nyDot0s
        counter = counter + 1;
        X0s_mat(counter,:) = [x0_range(kk), 0, 0, 0, yDot0_range(jj), 0];
    end
end
% 
% [x0_range_mesh,yDot0_range_mesh] = meshgrid(x0_range,yDot0_range);

tf = 20*pi;
timeVec = [0, tf];

%%% Preallocate
allCrossings = {};

parfor X0_index = 1:nTotal
    
    X0_i = X0s_mat(X0_index,:)';
    
    [time_i, X_i, time_i_event, X_i_event, index_i_event] =...
        ode113(@Int_CR3Bn, timeVec, X0_i, options_NegXAxisCrossing, prms);
    
    X_events_keep = [];
    t_events_keep = [];
    
    for kk = 1:size(X_i_event,1)
        if X_i_event(kk,1) < 0
            X_events_keep = [X_events_keep; X_i_event(kk,:)];
            t_events_keep = [t_events_keep; time_i_event(kk)];
        end
    end
    
    if isempty(X_events_keep) == 0
        allCrossings{X0_index}.X0s         = X0_i';
        allCrossings{X0_index}.X_events    = X_events_keep;
        allCrossings{X0_index}.time_events = t_events_keep;
    else
        allCrossings{X0_index}.X0s         = NaN(1,6);
        allCrossings{X0_index}.X_events    = NaN(1,6);
        allCrossings{X0_index}.time_events = NaN;
    end
    
end


allCrossings_structArray = [allCrossings{:}];

X0s = vertcat(allCrossings_structArray(:).X0s);
X_events = vertcat(allCrossings_structArray(:).X_events);
time_events = vertcat(allCrossings_structArray(:).time_events);

figure; hold all
plot(X_events(:,1),X_events(:,5),'k.')
PlotBoi2('$x_n$','$\dot{y}_n$',20,'LaTex')










% figure; hold all
% plotSecondary(secondary)
% plotPrimary(primary,secondary)
% axis equal
% PlotBoi3Norm
% 
% X0 = [-0.025399096010594 - (1-secondary.MR);0;0;0;0.041594674657154;0];
% tf = 50*pi;
% 
% [T_n, X_n] = ode113(@Int_CR3Bn, [0 tf], X0, options, prms);
% 
% plot3(X_n(:,1),X_n(:,2),X_n(:,3),'b')
% view(0,90)




% ========================================================================
%%% Closeout
% ========================================================================
tocWhole = toc(ticWhole);
fprintf('Elapsed time: %1.4f seconds\n',tocWhole)
















