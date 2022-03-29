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
%%% Integration options
tol = 1e-13;
options = odeset('RelTol',tol,'AbsTol',tol);


%%% Set primary & secondary
[primary, secondary] = assignPrimaryAndSecondary_CR3BP('Saturn_Enceladus', bodies);

%%% Factor for normalizing distances
rNorm = secondary.a; % n <-> km

%%% Setting parameters structure
prms.u = secondary.MR;
prms.R1 = primary.R / rNorm;
prms.R2 = secondary.R_n;
prms.J2p = primary.J2;
% prms.J4p = primary.J4;
% prms.J6p = primary.J6;
prms.J2s = secondary.J2;

%%% Getting normalized mean motion
tNorm = sqrt((rNorm^3)/(bodies.constants.G*(primary.mass + secondary.mass)));
prms.n = secondary.meanMot * tNorm;

% prms_theirMeanMot = prms;
% prms_theirMeanMot.n = 1 + 3*(prms.J2p*prms.R1*prms.R1 + prms.J2s*prms.R2*prms.R2)/2;

%%% Matching the review paper
% rNorm = 238042;
% prms.u = 0.1899309048e-6;
% prms.R1 = 60268.0/rNorm;
% prms.R2 = 252.1/rNorm;
% prms.J2p = 1.6298e-2;
% prms.J2s = 2.5e-3;
% prms.n = 1 + 3*(prms.J2p*prms.R1*prms.R1 + prms.J2s*prms.R2*prms.R2)/2;


% %%% Set primary & secondary
% [primary, secondary] = assignPrimaryAndSecondary_CR3BP('Jupiter_Europa', bodies);
% 
% %%% Factor for normalizing distances
% rNorm = secondary.a; % n <-> km
% 
% %%% Getting normalized mean motion
% tNorm = sqrt((rNorm^3)/(bodies.constants.G*(primary.mass + secondary.mass)));
% prms.n = secondary.meanMot * tNorm;
% 
% %%% Setting parameters structure
% prms.u = secondary.MR;
% prms.R1 = primary.R / rNorm;
% prms.R2 = secondary.R_n;
% prms.J2p = primary.J2;
% %     prms.J4p = primary.J4;
% %     prms.J6p = primary.J6;
% prms.J2s = secondary.J2;


%%% Equillibrium Points
% rLPs_n_van = EquilibriumPoints(prms.u, prms.n, 1:3);
rLPs_n_van = EquilibriumPoints(prms.u, 1, 1:3);
rLPs_n_ZH  = collinearEquilibriumPoints_ZH(prms);
rLPs_n_ZH_theirMeanMot  = collinearEquilibriumPoints_ZH(prms_theirMeanMot);
rLPs_n_ZH_meanMotEquals1  = collinearEquilibriumPoints_ZH(prms_n1);

dL123_n = rLPs_n_ZH(:,1) - rLPs_n_van(:,1);
dL123 = dL123_n .* rNorm;

fprintf('dL1n: %1.3e,\tdL2n: %1.3e,\tdL3n: %1.3e\n', dL123_n(1), dL123_n(2), dL123_n(3))
fprintf('dL1:  %1.3f,\t\tdL2:  %1.3f,\t\tdL3:  %1.3f km\n', dL123(1), dL123(2), dL123(3))



A1 = prms.J2p * (prms.R1^2);
A2 = prms.J2s * (prms.R2^2);


r1 = @(x,y,z) sqrt((x+prms.u)^2 + y^2 + z^2);
r2 = @(x,y,z) sqrt((x-1+prms.u)^2 + y^2 + z^2);

U_full_J2  = @(x,y,z) (1/2)*(prms.n^2)*(x^2+y^2) + (1-prms.u)/r1(x,y,z) + prms.u/r2(x,y,z)  - (1-prms.u)*A1*(3*z*z-r1(x,y,z)^2)/(2*(r1(x,y,z)^5)) - prms.u*A2*(3*z*z-r2(x,y,z)^2)/(2*(r2(x,y,z)^5));
U_full_van = @(x,y,z) (1/2)*(x^2+y^2) + (1-prms.u)/r1(x,y,z) + prms.u/r2(x,y,z);

% U_J2only = @(x) (1-prms.u)*A1/(2*(abs(x+prms.u)^3)) + prms.u*A2/(2*(abs(x-1+prms.u)^3));
% U_J2 = @(x)  .5*(prms.n^2)*(x^2)    + (1-prms.u)/abs(x+prms.u) + prms.u/abs(x-1+prms.u) + (1-prms.u)*A1/(2*(abs(x+prms.u)^3)) + prms.u*A2/(2*(abs(x-1+prms.u)^3));
% U_van = @(x) .5*(x^2) +               (1-prms.u)/abs(x+prms.u) + prms.u/abs(x-1+prms.u);



xs = [linspace(-1.1, 0.9, 1000000), linspace(0.9,1.1,1000000)];

Us_J2 = zeros(length(xs),1);
Us_van = zeros(length(xs),1);


% Us_J2only = zeros(length(xs),1);
for kk = 1:length(xs)
    Us_J2(kk) = U_full_J2(xs(kk),0,0);
    Us_van(kk) = U_full_van(xs(kk),0,0);
%     Us_J2only(kk) = U_J2only(xs(kk));
end



figure; hold all
plot(xs,Us_van,'b','linewidth',1.5)
ylim([-1000,1000])
PlotBoi2('$x_n$','$U_{van}$',20,'LaTex')

figure; hold all
plot(xs,Us_J2,'r','linewidth',1.5)
ylim([-1000,1000])
PlotBoi2('$x_n$','$U_{J_2}$',20,'LaTex')


% x_L22 = 1.000926133887796;
% ys = linspace(-.1,.1,1000000);
% Us_yaxis_J2 = zeros(length(ys),1);
% for kk = 1:length(ys)
%     Us_yaxis_J2(kk) = U_full_J2(x_L22,ys(kk),0);
% end
% figure; hold all
% plot(ys,Us_yaxis_J2, 'k', 'linewidth',1.5)
% PlotBoi2('$y_n$','$U_{J_2}$',20,'LaTex')

% figure; hold all
% plot(xs,Us_J2only,'m','linewidth',1.5)
% ylim([-1000,1000])
% PlotBoi2('$x_n$','$U_{J_2}$ Only',20,'LaTex')

% -------------------------------------------------
%%% Coding up their potential function
% -------------------------------------------------
U_them = @(x,y,z,mu) (1-mu)/r1(x,y,z) + (1-mu)*A1*(r1(x,y,z)^2 - 3*z*z)/(2*(r1(x,y,z)^5)) + mu/r2(x,y,z) + mu*A2*(r2(x,y,z)^2 - 3*z*z)/(2*(r2(x,y,z)^5)) + (prms.n^2)*((1-mu)*(r1(x,y,z)^2) + mu*(r2(x,y,z)^2) - z*z)/2;

Us_them = zeros(length(xs),1);
for kk = 1:length(xs)
    Us_them(kk) = U_them(xs(kk),0,0,prms.u);
end



figure; hold all
plot(xs,Us_J2,'r','linewidth',2)
plot(xs,Us_them - .5*prms.n*prms.n*(prms.u - prms.u^2),'k','linewidth',1)
ylim([-1000,1000])
PlotBoi2('$x_n$','$U_{J_2}$',20,'LaTex')
plot(rLPs_n_ZH(1,1), U_full_J2(rLPs_n_ZH(1,1),0,0), '^')
plot(rLPs_n_ZH(2,1), U_full_J2(rLPs_n_ZH(2,1),0,0), '^')
plot(rLPs_n_ZH(3,1), U_full_J2(rLPs_n_ZH(3,1),0,0), '^')


plot(rLPs_n_van(1,1), U_full_van(rLPs_n_van(1,1),0,0), 'b^')
plot(rLPs_n_van(2,1), U_full_van(rLPs_n_van(2,1),0,0), 'b^')
plot(rLPs_n_van(3,1), U_full_van(rLPs_n_van(3,1),0,0), 'b^')

% -------------------------------------------------
%%% Testing these points in integration
% -------------------------------------------------
time = [0, 2*pi];
X01 = [rLPs_n_ZH(1,1); 0; 0; 0; 0 ;0];
[~, X1] = ode113(@Int_CR3Bn_ZH, time, X01, options, prms);
X02 = [rLPs_n_ZH(2,1); 0; 0; 0; 0 ;0];
[~, X2] = ode113(@Int_CR3Bn_ZH, time, X02, options, prms);

figure
plot3(X1(:,1),X1(:,2),X1(:,3),'b')
PlotBoi3_CR3Bn(20)
% plotBody2(prms.R2,[1-prms.u,0,0],colors.white,colors.black,2,0.5)

figure
plot3(X2(:,1),X2(:,2),X2(:,3),'r')
PlotBoi3_CR3Bn(20)
% plotBody2(prms.R2,[1-prms.u,0,0],colors.white,colors.black,2,0.5)



% A1 = get_Amat_CR3BP_J2(prms.u, [1.003738; 0; 0]', prms.J2p, prms.J2s, prms.R1, prms.R2, prms.n);
% A2 = get_Amat_CR3BP_J2(prms.u, rLPs_n_ZH(2,:), prms.J2p, prms.J2s, prms.R1, prms.R2, prms.n);
% A3 = get_Amat_CR3BP_J2(prms.u, rLPs_n_ZH(1,:), prms.J2p, prms.J2s, prms.R1, prms.R2, prms.n);
% A4 = get_Amat_CR3BP_J2(prms.u, rLPs_n_ZH(1,:) + [1e-3, 0, 0], prms.J2p, prms.J2s, prms.R1, prms.R2, prms.n);
% 
% [v1,e1] = eig(A1);
% e1 = diag(e1)
% 
% [v2,e2] = eig(A2);
% e2 = diag(e2)
% 
% [v3,e3] = eig(A3);
% e3 = diag(e3)
% 
% [v4,e4] = eig(A4);
% e4 = diag(e4)
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





% %%% Manual results
% %L1_van - L1_J2
% dL1n_manual = 0.995911 - 0.996019
% %L2_van - L2_J2
% dL2n_manual = 1.003738 - 1.003991
% %L3_van - L3_J2
% dL3n_manual = (-0.99948) - (-1.000001)
% 
% dL1_manual = dL1n_manual*rNorm
% dL2_manual = dL2n_manual*rNorm
% dL3_manual = dL3n_manual*rNorm





