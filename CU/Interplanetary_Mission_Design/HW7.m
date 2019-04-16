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
run_p1 = 0;
run_p2 = 0;
run_p3 = 1;

% ========================================================================
%%% Problem 1
% ========================================================================
if run_p1 == 1
% ------------------------------------
%%% System 
% ------------------------------------
AU = 1.49597870691e8; % km

primary.u = 1.32712440018e11; % km^3/s^2

secondary.u = 4.035032351966808e5; % km^3/s^2

r_earthMoonBC = [150000000.0; 6000.0; 1450.0]; % km
v_earthMoonBC = [0.00075; 0.08; 0.019]; % km/s

t_days = 450; % days
t_sec = 450*86400;

%%% Normalizing constants
rNorm = AU;         % n <-> km
tNorm = (bodies.earth.Tp)/(2*pi); % n <-> sec
vNorm = rNorm / tNorm;       % n <-> km/sec


normalizedPosition = r_earthMoonBC./rNorm
normalizedVelocity = v_earthMoonBC./vNorm
normalizedTime     = t_sec/tNorm


end % run_p1


% ========================================================================
%%% Problem 2
% ========================================================================
if run_p1 == 0 && run_p2 == 1
% ------------------------------------
%%% System 
% ------------------------------------
primary = bodies.earth;
secondary = bodies.moon;

secondary.a = 384747.962856037;
secondary.MR = 0.012150585609624;

%%% Normalizing constants
rNorm = secondary.a;         % n <-> km
tNorm = 1/secondary.meanMot; % n <-> sec
vNorm = rNorm / tNorm;       % n <-> km/sec

% ------------------------------------
%%% Setup
% ------------------------------------
%%% ICs given
X01 = [1.2; 0; 0; 0; -1.06110124; 0];
t_f1 = 6.20628;

X02 = [0.85; 0; 0.17546505; 0; 0.2628980369; 0];
t_f2 = 2.5543991;

%%% Integration times
t_i = 0; % sec
n_dt = 10000;
time01_n = linspace(t_i,t_f1,n_dt);
time02_n = linspace(t_i,t_f2,n_dt);

%%% Choosing ode tolerance
tol = 1e-13;

%%% Setting integrator options
options = odeset('RelTol',tol,'AbsTol',tol);

%%% Setting necessary parameters for integration
prms.u = secondary.MR;

% ------------------------------------
%%% Integrating
% ------------------------------------
[time1_n, X1_BCR_n] = ode113(@Int_CR3Bn, time01_n, X01, options, prms);
[time2_n, X2_BCR_n] = ode113(@Int_CR3Bn, time02_n, X02, options, prms);



% ------------------------------------
%%% Plotting
% ------------------------------------
figure; hold all
p1 = plot3(X1_BCR_n(:,1),X1_BCR_n(:,2),X1_BCR_n(:,3),'b','linewidth',2);
plot3(X1_BCR_n(1,1),X1_BCR_n(1,2),X1_BCR_n(1,3),'bo','markersize',10)
plot3(X1_BCR_n(end,1),X1_BCR_n(end,2),X1_BCR_n(end,3),'bx','markersize',10)
plotBodyTexture3(primary.R/rNorm, [-secondary.MR, 0, 0],primary.img)
plotBodyTexture3(secondary.R_n, [1-secondary.MR, 0, 0],secondary.img)
PlotBoi3('$X_n$','$Y_n$','$Z_n$',18,'LaTex')
axis equal
legend([p1],'Trajectory 1')

figure; hold all
p2 = plot3(X2_BCR_n(:,1),X2_BCR_n(:,2),X2_BCR_n(:,3),'r','linewidth',2);
plot3(X2_BCR_n(1,1),X2_BCR_n(1,2),X2_BCR_n(1,3),'ro','markersize',10)
plot3(X2_BCR_n(end,1),X2_BCR_n(end,2),X2_BCR_n(end,3),'rx','markersize',10)
plotBodyTexture3(secondary.R_n, [1-secondary.MR, 0, 0],secondary.img)
PlotBoi3('$X_n$','$Y_n$','$Z_n$',18,'LaTex')
axis equal
legend([p2],'Trajectory 2')

end % run_p1 == 0 && run_p2 == 1



% ========================================================================
%%% Problem 3
% ========================================================================
if run_p1 == 0 && run_p2 == 0 && run_p3 == 1
% ------------------------------------
%%% System 
% ------------------------------------
primary = bodies.earth;
secondary = bodies.moon;

secondary.a = 384747.962856037;
secondary.MR = 0.012150585609624;

%%% Normalizing constants
rNorm = secondary.a;         % n <-> km
tNorm = 1/secondary.meanMot; % n <-> sec
vNorm = rNorm / tNorm;       % n <-> km/sec

% ------------------------------------
%%% Setup
% ------------------------------------
%%% ICs given
X0_a = [-0.083; -0.03; 0.01; 3.53; -3.1; -0.1];
tf_a = 26;

X0_b = [0.05; -0.05; 0; 4.0; 2.6; 0];
tf_b = 15;

X0_c = [0.875; 0; 0.1914; 0; 0.23234; 0];
tf_c = 15;

X0_d = [-0.05; -0.02; 0; 4.09; -5.27; 0];
tf_d = 2.5;

%%% Integration times
t_i = 0; % sec
n_dt = 10000;
time0_a_n = linspace(t_i,tf_a,n_dt);
time0_b_n = linspace(t_i,tf_b,n_dt);
time0_c_n = linspace(t_i,tf_c,n_dt);
time0_d_n = linspace(t_i,tf_d,n_dt);

%%% Choosing ode tolerance
tol = 1e-13;

%%% Setting integrator options
options = odeset('RelTol',tol,'AbsTol',tol);

%%% Setting necessary parameters for integration
prms.u = secondary.MR;

% ------------------------------------
%%% Integrating
% ------------------------------------
[timea_n, Xa_BCR_n] = ode113(@Int_CR3Bn, time0_a_n, X0_a, options, prms);
[timeb_n, Xb_BCR_n] = ode113(@Int_CR3Bn, time0_b_n, X0_b, options, prms);
[timec_n, Xc_BCR_n] = ode113(@Int_CR3Bn, time0_c_n, X0_c, options, prms);
[timed_n, Xd_BCR_n] = ode113(@Int_CR3Bn, time0_d_n, X0_d, options, prms);


% ------------------------------------
%%% Rotating into inertial
% ------------------------------------
%%% Earth
[ rEartha_BCI ] = r_BCR2BCI( repmat([-secondary.MR, 0, 0],length(timea_n),1), timea_n.*tNorm, secondary.meanMot );
[ rEarthb_BCI ] = r_BCR2BCI( repmat([-secondary.MR, 0, 0],length(timeb_n),1), timeb_n.*tNorm, secondary.meanMot );
[ rEarthc_BCI ] = r_BCR2BCI( repmat([-secondary.MR, 0, 0],length(timec_n),1), timec_n.*tNorm, secondary.meanMot );
[ rEarthd_BCI ] = r_BCR2BCI( repmat([-secondary.MR, 0, 0],length(timed_n),1), timed_n.*tNorm, secondary.meanMot );

%%% Trajectory a
[ ra_BCI ] = r_BCR2BCI( Xa_BCR_n(:,1:3), timea_n.*tNorm, secondary.meanMot );
ra_ECI = ra_BCI - rEartha_BCI;
[ va_BCI ] = v_BCR2BCI( Xa_BCR_n(:,4:6), ra_BCI, timea_n, secondary.meanMot );
Xa_ECI_n = [ra_ECI, va_BCI];

%%% Trajectory b
[ rb_BCI ] = r_BCR2BCI( Xb_BCR_n(:,1:3), timeb_n.*tNorm, secondary.meanMot );
rb_ECI = rb_BCI - rEarthb_BCI;
[ vb_BCI ] = v_BCR2BCI( Xb_BCR_n(:,4:6), rb_BCI, timeb_n, secondary.meanMot );
Xb_ECI_n = [rb_ECI, vb_BCI];

%%% Trajectory c
[ rc_BCI ] = r_BCR2BCI( Xc_BCR_n(:,1:3), timec_n.*tNorm, secondary.meanMot );
rc_ECI = rc_BCI - rEarthc_BCI;
[ vc_BCI ] = v_BCR2BCI( Xc_BCR_n(:,4:6), rc_BCI, timec_n, secondary.meanMot );
Xc_ECI_n = [rc_ECI, vc_BCI];

%%% Trajectory d
[ rd_BCI ] = r_BCR2BCI( Xd_BCR_n(:,1:3), timed_n.*tNorm, secondary.meanMot );
rd_ECI = rd_BCI - rEarthd_BCI;
[ vd_BCI ] = v_BCR2BCI( Xd_BCR_n(:,4:6), rd_BCI, timed_n, secondary.meanMot );
Xd_ECI_n = [rd_ECI, vd_BCI];

% ------------------------------------
%%% Plotting
% ------------------------------------
%%% Figure (a)
figure('position',[440 520 1001 278])
subplot(1,2,1); hold all
title('Rotating')
pa1 = plot3(Xa_BCR_n(:,1),Xa_BCR_n(:,2),Xa_BCR_n(:,3),'b','linewidth',2);
plotBothBodiesHere(primary, secondary, rNorm)

subplot(1,2,2); hold all
title('Inertial')
pa2 = plot3(Xa_ECI_n(:,1),Xa_ECI_n(:,2),Xa_ECI_n(:,3),'b','linewidth',2);
plotPrimaryHere_ECI(primary, secondary, rNorm)

%%% Figure (b)
figure('position',[440 520 1001 278])
subplot(1,2,1); hold all
title('Rotating')
pb1 = plot3(Xb_BCR_n(:,1),Xb_BCR_n(:,2),Xb_BCR_n(:,3),'r','linewidth',2);
plotBothBodiesHere(primary, secondary, rNorm)

subplot(1,2,2); hold all
title('Inertial')
pb2 = plot3(Xb_ECI_n(:,1),Xb_ECI_n(:,2),Xb_ECI_n(:,3),'r','linewidth',2);
plotPrimaryHere_ECI(primary, secondary, rNorm)

%%% Figure (c)
figure('position',[440 520 1001 278])
subplot(1,2,1); hold all
title('Rotating')
pc1 = plot3(Xc_BCR_n(:,1),Xc_BCR_n(:,2),Xc_BCR_n(:,3),'g','linewidth',2);
plotSecondaryHere(secondary)

subplot(1,2,2); hold all
title('Inertial')
pc2 = plot3(Xc_ECI_n(:,1),Xc_ECI_n(:,2),Xc_ECI_n(:,3),'g','linewidth',2);
plotPrimaryHere_ECI(primary, secondary, rNorm)

%%% Figure (d)
figure('position',[440 520 1001 278])
subplot(1,2,1); hold all
title('Rotating')
pd1 = plot3(Xd_BCR_n(:,1),Xd_BCR_n(:,2),Xd_BCR_n(:,3),'c','linewidth',2);
plotBothBodiesHere(primary, secondary, rNorm)

subplot(1,2,2); hold all
title('Inertial')
pd2 = plot3(Xd_ECI_n(:,1),Xd_ECI_n(:,2),Xd_ECI_n(:,3),'c','linewidth',2);
plotPrimaryHere_ECI(primary, secondary, rNorm)













end % run_p1 == 0 && run_p2 == 0 && run_p3 == 1













function plotBothBodiesHere(primary, secondary, rNorm)
plotBodyTexture3(primary.R/rNorm, [-secondary.MR, 0, 0],primary.img)
plotBodyTexture3(secondary.R_n, [1-secondary.MR, 0, 0],secondary.img)
PlotBoi3('$X_n$','$Y_n$','$Z_n$',18,'LaTex')
axis equal
end


function plotPrimaryHere_ECI(primary, secondary, rNorm)
plotBodyTexture3(primary.R/rNorm, [0, 0, 0],primary.img)
PlotBoi3('$X_n$','$Y_n$','$Z_n$',18,'LaTex')
axis equal
end


function plotSecondaryHere(secondary)
plotBodyTexture3(secondary.R_n, [1-secondary.MR, 0, 0],secondary.img)
PlotBoi3('$X_n$','$Y_n$','$Z_n$',18,'LaTex')
axis equal
end






















