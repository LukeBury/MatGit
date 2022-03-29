% ========================================================================
%%% Description
% ========================================================================
% 

% Created: 05/25/20
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
primary = bodies.saturn;    secondary = bodies.enceladus;




%%% Normalizing constants
rNorm = secondary.a;         % n <-> km


tol = 1e-13;
options = odeset('RelTol',tol,'AbsTol',tol);


prms.u    = secondary.MR;
MR = secondary.MR;

m2 = secondary.mass;
m1 = primary.mass;
mTot = m1 + m2;
prms.up = m1*bodies.constants.G;
prms.us = m2*bodies.constants.G;
prms.m1 = m1;
prms.m2 = m2;

%%% Making mean motion a self-contained definition
989
secondary.meanMot = sqrt((m1+m2)*bodies.constants.G/(secondary.a^3));
tNorm = 1/secondary.meanMot; % n <-> sec
vNorm = rNorm / tNorm;       % n <-> km/sec


% prms.up   = primary.u;
prms.Rp   = primary.R;
prms.R1   = primary.R / rNorm;
prms.J2p  = primary.J2;

% prms.us      = secondary.u;
prms.Rs      = secondary.R;
prms.R2      = secondary.R_n;
prms.J2s     = secondary.J2;
prms.a       = secondary.a;
prms.meanMot = secondary.meanMot;

%%% Equillibrium Points
rLPs_n = EquilibriumPoints(prms.u);
rLPs_ZH_n = collinearEquilibriumPoints_ZH(prms);


t0 = 0;
n_t = 10000;
% tf = 3600*24*1;
tf_n = 3.0887514754980483;

Lp = 1;

%%% BaCI_n and PCI_n integration
% X0_BaCR_n = [1.0044576806787502,-0.0000000000000000,0.0004397507443672,-0.0000000000000013,-0.0036121436101641,0.0000000000000022];
% X0_BaCR_n = [rLPs_ZH_n(Lp,:)+[0.001/rNorm,0,0],zeros(1,3)];
X0_BaCR_n = [rLPs_ZH_n(Lp,:)+[0/rNorm,0,0],zeros(1,3)];
% X0_PCR_n  = X0_BaCR_n + [prms.u, zeros(1,5)];

% X0_BaCI_n = X_BCR2BCI(X0_BaCR_n, 0, 1);
% % X0_PCI_n  = X_BCR2BCI(X0_PCR_n, 0, 1);
% 
% X0_BaCI = [X0_BaCI_n(1:3).*rNorm, X0_BaCI_n(4:6).*vNorm];



X0_BaCR = [X0_BaCR_n(1:3).*rNorm, X0_BaCR_n(4:6).*vNorm];
X0_BaCI = X_BCR2BCI(X0_BaCR, 0, secondary.meanMot);



time0_n = linspace(t0,tf_n,n_t);
time0 = linspace(t0,tf_n*tNorm,n_t);

% ------------------------------------
%%% Integrate
% ------------------------------------
% [time_out, X_BaCI_n] = ode113(@Int_3BI_J2_testingMeanMotionDerivation,      time0_n, X0_BaCI_n', options, prms);
% [time_out, X_PCI_n] = ode113(@Int_3BI_J2_testingMeanMotionDerivation_PCIn, time0_n, X0_PCI_n',  options, prms);

[time_out, X_BaCI] = ode113(@Int_3BI_J2_testingMeanMotionDerivation, time0, X0_BaCI', options, prms);


% X_BaCR_n = X_BCI2BCR(X_BaCI_n, time_out, 1);
% 
% figure; hold all
% axis equal
% PlotBoi3_CR3Bn(20)
% plotSecondary(secondary)
% plot3(X_BaCR_n(:,1),X_BaCR_n(:,2),X_BaCR_n(:,3),'b')

% X_PCR_n = X_BCI2BCR(X_PCI_n, time_out, 1);
% 
% figure; hold all
% axis equal
% PlotBoi3_CR3Bn(20)
% % plotSecondary(secondary)
% plot3(X_PCR_n(:,1),X_PCR_n(:,2),X_PCR_n(:,3),'r')

X_BaCR = X_BCI2BCR(X_BaCI, time_out, prms.meanMot);

figure; hold all
axis equal
PlotBoi3('$x$','$y$','$z$',20,'LaTex')
plot3(X_BaCR(:,1),X_BaCR(:,2),X_BaCR(:,3),'b')


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









% My computed meanMot from n=sqrt(...) ... n = 5.305341284382199e-05
% The mean motion I have stored in bodyData: n = 5.307334465496030e-05
% computed / bodyData = 0.999624447803923

% Mean motion computed their way is: n = 5.309222050748919e-05
% from n = sqrt(bodies.constants.G * (m1+m2) * (1/(secondary.a^3) + 3*(prms.J2p*primary.R*primary.R + prms.J2s*secondary.R*secondary.R)/(2*(secondary.a^5))))
% theirs / bodyData = 1.000355655982332
% bodyData / theirs = 0.999644470463875




