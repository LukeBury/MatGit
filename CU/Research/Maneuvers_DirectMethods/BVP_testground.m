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
%%% Setup
% ========================================================================
%%% system
primary = bodies.jupiter;
secondary = bodies.europa;

%%% Collinear equillibria
L123 = EquilibriumPoints(secondary.MR,1:3); % [3x3] of L1, L2, L3 normalized BCR coordinates

%%% Normalizing constants
rNorm = secondary.a;         % n <-> km
tNorm = 1/secondary.meanMot; % n <-> sec
vNorm = rNorm / tNorm;       % n <-> km/sec

%%% Time
t_i = 0; % sec
t_f = 4*pi; % Long bc events are watching for impact or escape
n_dt = 10000;
time0_n = linspace(t_i,t_f,n_dt);

%%% Choosing ode tolerance
tol = 1e-13;

%%% Setting integrator options
options = odeset('Event',@event_Impact_CR3Bn, 'RelTol',tol,'AbsTol',tol);

%%% Setting necessary parameters for integration
prms.u = secondary.MR;
prms.R2_n = secondary.R_n;
prms.L1x = L123(1,1);
prms.L2x = L123(2,1);

%%% ICs (impact LatLon = [11.713313624820113, -17.884817025059917] )
X0_X_truth = [1.0204617015266166,-0.0019850725823437,0.0000000000000000,-0.0010976413130697,-0.0030157447223195,0.0098771743854318]';
introducedPositionError = [0; 1; 1].*(2.0)./rNorm;
X0_X_error = X0_X_truth + [introducedPositionError; 0; 0; 0];
% introducedVelocityError = -[0; 1; 1].*(1e-5)./vNorm;
% X0_X_error = X0_X_truth + [0; 0; 0; introducedVelocityError];
% introducedStateError = [[0; 0.001; 0.001]./rNorm; [0; 1e-9; 1e-9]./vNorm];
% X0_X_error = X0_X_truth + introducedStateError;


stm0 = reshape(eye(6),36,1);

X0 = [X0_X_error; stm0];


[time_n, X_BCR_n_nom] = ode113(@int_CR3BnSTM, time0_n, X0, options, prms);

[lat0, lon0] = BCR2latlon(X_BCR_n_nom(end,1:3), 'secondary', secondary.MR);

if (norm(X_BCR_n_nom(end,1:3) - [1-secondary.MR,0,0]) - secondary.R_n) > (1/rNorm)
    warning('no impact')
    return
end
targetLat = 11.713313624820113;
targetLon = -17.884817025059917;

timeLimit_n = time_n(end);
%%% Time
t_i = 0; % sec
t_f = timeLimit_n*1.05; % Long bc events are watching for impact or escape
n_dt = 10000;
time0_n = linspace(t_i,t_f,n_dt);


missLatLon0 = [lat0 - targetLat, lon0 - targetLon];

rTarget = [0.997807465412678; 0.000699369559596; 0.000472158462809];
% vTarget = [-0.022327519272868; -0.126350703871528; 0.041233558870785];
% stateTarget = [0.997807465412667; 0.000699369559550; 0.000472158462823; -0.022327519271346; -0.126350703871916; 0.041233558870487];

dXf = X_BCR_n_nom(end,1:3)' - rTarget;
% dXf = X_BCR_n_nom(end,4:6)' - vTarget;
% dXf = X_BCR_n_nom(end,1:6)' - stateTarget;

stm_tf_t0 = reshape(X_BCR_n_nom(end,7:42),6,6);
stm_tf_t0 = stm_tf_t0(1:3,1:3);
% stm_tf_t0 = stm_tf_t0(4:6,4:6);
stm_t0_tf = inv(stm_tf_t0);

dX0 = stm_t0_tf * dXf;


X0_new = [X0_X_error - [dX0; 0; 0; 0]; stm0];
% X0_new = [X0_X_error - [0; 0; 0; dX0]; stm0];
% X0_new = [X0_X_error - dX0; stm0];
% X0_new = [X0_X + dX0; stm0];
[time_n, X_BCR_n_new] = ode113(@int_CR3BnSTM, time0_n, X0_new, options, prms);
[lat_new, lon_new] = BCR2latlon(X_BCR_n_new(end,1:3), 'secondary', secondary.MR);


missLatLonf = [lat_new - targetLat, lon_new - targetLon];

dX0_scale = 1;
n_iterations = 30;

latLons_iter = zeros(n_iterations,2);

figure; hold all
for kk = 1:n_iterations

dXf = X_BCR_n_new(end,1:3)' - rTarget;
% dXf = X_BCR_n_new(end,4:6)' - vTarget;
% dXf = X_BCR_n_nom(end,1:6)' - stateTarget;

stm_tf_t0 = reshape(X_BCR_n_new(end,7:42),6,6);
stm_tf_t0 = stm_tf_t0(1:3,1:3);
% stm_tf_t0 = stm_tf_t0(4:6,4:6);
stm_t0_tf = inv(stm_tf_t0);

dX0 = stm_t0_tf * dXf;

X0_new = [X0_new(1:6) - [dX0; 0; 0; 0].*dX0_scale; stm0];
% X0_new = [X0_new(1:6) - [0; 0; 0; dX0].*dX0_scale; stm0];
% X0_new = [X0_new(1:6) - dX0.*dX0_scale; stm0];

% X0_new = [X0_X + dX0; stm0];
[time_n, X_BCR_n_new] = ode113(@int_CR3BnSTM, time0_n, X0_new, options, prms);
[lat_new, lon_new] = BCR2latlon(X_BCR_n_new(end,1:3), 'secondary', secondary.MR);
latLons_iter(kk,:) = [lat_new, lon_new];
missLatLonf = [lat_new - targetLat, lon_new - targetLon];

fprintf('-------------------------------------------- %1d\n',kk)
% fprintf('Latitude Miss: \t%1.4f\n',missLatLonf(1))
% fprintf('Longitude Miss:\t%1.4f\n',missLatLonf(2))
fprintf('Norm of Miss:  \t%1.4f\n',norm(missLatLonf))
fprintf('Norm of dX0:  \t%1.4e\n',norm(dXf))

plot(kk,norm(missLatLonf),'ro','markersize',10)
plot([0 kk],[norm(missLatLon0) norm(missLatLon0)],'b')
end

[ colorMatrix ] = colorScale([colors.std.mag; colors.std.cyan],n_iterations );
figure; hold all
% plot(latLons_iter(:,2),latLons_iter(:,1),'rx')
scatter(latLons_iter(:,2),latLons_iter(:,1),30,colorMatrix,'filled','marker','o')
plot(targetLon, targetLat, 'bo','markersize',15)


figure; hold all
plot3(X_BCR_n_nom(:,1),X_BCR_n_nom(:,2),X_BCR_n_nom(:,3),'b')
plot3(X_BCR_n_new(:,1),X_BCR_n_new(:,2),X_BCR_n_new(:,3),'r')
plotBodyTexture3(secondary.R_n, [1-secondary.MR, 0, 0], secondary.img)
% axis equal
PlotBoi3('x','y','z',18,'LaTex')






