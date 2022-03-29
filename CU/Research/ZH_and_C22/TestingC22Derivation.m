% ========================================================================
%%% Description
% ========================================================================
% 

% Created: 3/25/22
% Author : Luke Bury, luke.bury@jpl.nasa.gov
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
%%% 
% -------------------------------------------------

% [rNorm, tNorm, vNorm] = cr3bp_norms(primary, secondary, bodies.constants.G);



syms mu r R J2 x y z C22 real

V_r_formB = (mu/r)*(1 - ((J2*R*R)/(r^2))*(3*z*z/(2*r*r)-1/2) + ((3*R*R*C22)/(r^2))*((x^2-y^2)/(x^2+y^2))*(1 - z*z/(r^2)));
% pretty(V_r)

V_r_formC = mu/r  - ((mu*R*R*J2)/(r^3))*(3*z*z/(2*r*r) - 1/2) + ((3*mu*R*R*C22)/(r^3))*((x^2-y^2)/(x^2+y^2) - (x^2-y^2)*(z*z)/((x^2+y^2)*r*r));

V_r_formD = mu/r - 3*mu*R*R*J2*z*z/(2*(r^5)) + mu*R*R*J2/(2*(r^3)) + 3*mu*R*R*C22*(x^2-y^2)/((x^2+y^2)*(r^3)) - 3*mu*R*R*C22*z*z*(x^2-y^2)/((x^2+y^2)*(r^5));

V_r_formE = mu*((x^2+y^2+z^2)^(-1/2)) - (3/2)*mu*R*R*J2*z*z*((x^2+y^2+z^2)^(-5/2)) + (1/2)*mu*R*R*J2*((x^2+y^2+z^2)^(-3/2)) + 3*mu*R*R*C22*(x^2-y^2)*((x^2+y^2)^(-1))*((x^2+y^2+z^2)^(-3/2)) - 3*mu*R*R*C22*z*z*(x^2-y^2)*((x^2+y^2)^(-1))*((x^2+y^2+z^2)^(-5/2));

V_r_formB_subr = subs(V_r_formB, r, (x^2 + y^2 + z^2)^(1/2));
V_r_formC_subr = subs(V_r_formC, r, (x^2 + y^2 + z^2)^(1/2));
V_r_formD_subr = subs(V_r_formD, r, (x^2 + y^2 + z^2)^(1/2));


V_r_formB_ddx = diff(V_r_formB_subr, x);
V_r_formC_ddx = diff(V_r_formC_subr, x);
V_r_formD_ddx = diff(V_r_formD_subr, x);
V_r_formE_ddx = diff(V_r_formE, x);

V_r_formB_ddy = diff(V_r_formB_subr, y);

V_r_formB_ddz = diff(V_r_formB_subr, z);

ddx_formF = -mu*x/(r^3) + 15*mu*R*R*J2*z*z*x/(2*(r^7)) - 3*mu*R*R*J2*x/(2*(r^5)) + 3*mu*R*R*C22*(2*x/((x^2+y^2)*(r^3)) + (x^2-y^2)*(-2*x/(((x^2+y^2)^2)*(r^3)) - 3*x/((x^2+y^2)*(r^5)))) - 3*mu*R*R*C22*z*z*(2*x/((x^2+y^2)*(r^5)) + (x^2-y^2)*(-2*x/(((x^2+y^2)^2)*(r^5)) - 5*x/((x^2+y^2)*(r^7))));


ddx_hand = (-mu/r^3 + ((3*mu*R*R*J2)/(2*r^5))*(5*z*z/(r^2) -1) + ((3*mu*R*R*C22)/((x^2+y^2)*r^3))*(2 - 2*(x^2-y^2)/(x^2+y^2) - (3*(x^2-y^2) + 2*z*z)/(r^2) + 2*z*z*(x^2-y^2)/((x^2+y^2)*r*r) + (5*z*z*(x^2-y^2))/(r^4)))*x;
% ddx = (-mu/r^3 + ((6*mu*R*R*J2)/(2*r^5))*(5*z*z/(r^2) -1) + ((3*mu*R*R*C22)/((x^2+y^2)*r^3))*(z - 2*(x^2-y^2)/(x^2+y^2) + (2*z*z*(x^2-y^2) - 3*(x^2-y^2) - 2*z*z)/(r^2) + (5*z*z*(x^2-y^2))/(r^4)))*x;
% pretty(ddx)

ddy_hand = (-mu/r^3 + ((3*mu*R*R*J2)/(2*r^5))*(5*z*z/(r^2) -1) + ((3*mu*R*R*C22)/((x^2+y^2)*r^3))*(-2 - 2*(x^2-y^2)/(x^2+y^2) - (3*(x^2-y^2) - 2*z*z)/(r^2) + 2*z*z*(x^2-y^2)/((x^2+y^2)*r*r) + (5*z*z*(x^2-y^2))/(r^4)))*y;

ddz_hand = -mu*z/(r^3) + (3*mu*R*R*J2*z/(2*(r^7)))*(5*z*z - 3*r*r) + (15*mu*R*R*C22*(x^2-y^2)*z/((x^2+y^2)*(r^7)))*(z^2 - r^2);


ddx_formF_subr = subs(ddx_formF, r, (x^2 + y^2 + z^2)^(1/2));
ddx_hand_subr  = subs(ddx_hand,  r, (x^2 + y^2 + z^2)^(1/2));

ddy_hand_subr  = subs(ddy_hand,  r, (x^2 + y^2 + z^2)^(1/2));

ddz_hand_subr  = subs(ddz_hand,  r, (x^2 + y^2 + z^2)^(1/2));

% isequal(V_r_formB_ddx, ddx)

mu = 1e-2;
x = 1.1;
y = 0.2;
z = 0.15;
J2 = 0.002;
C22 = 0.0005;
R = 0.012;

ans_ddx_formF = vpa(subs(ddx_formF_subr), 16);
ans_ddx_hand = vpa(subs(ddx_hand_subr), 16);
ans_B = vpa(subs(V_r_formB_ddx), 16);
ans_C = vpa(subs(V_r_formC_ddx), 16);
ans_D = vpa(subs(V_r_formD_ddx), 16);
ans_E = vpa(subs(V_r_formE_ddx), 16);


ans_ddy_hand = vpa(subs(ddy_hand_subr), 16);
ans_B_ddy = vpa(subs(V_r_formB_ddy), 16);

ans_ddz_hand = vpa(subs(ddz_hand_subr), 16);
ans_B_ddz = vpa(subs(V_r_formB_ddz), 16);


fprintf('Difference ddx_hand Vr_B: %1.1e\n', ans_ddx_hand - ans_B)
fprintf('Difference ddx_hand Vr_C: %1.1e\n', ans_ddx_hand - ans_C)
fprintf('Difference ddx_hand Vr_D: %1.1e\n', ans_ddx_hand - ans_D)
fprintf('Difference ddx_hand Vr_E: %1.1e\n', ans_ddx_hand - ans_E)

fprintf('Difference Vr_B Vr_C: %1.1e\n', ans_B - ans_C)
fprintf('Difference Vr_B Vr_D: %1.1e\n', ans_B - ans_D)
fprintf('Difference Vr_B Vr_E: %1.1e\n', ans_B - ans_E)

fprintf('Difference Vr_C Vr_D: %1.1e\n', ans_C - ans_D)
fprintf('Difference Vr_C Vr_E: %1.1e\n', ans_C - ans_E)

fprintf('Difference Vr_D Vr_E: %1.1e\n', ans_D - ans_E)



fprintf('Difference ddx_F ddx_hand: %1.1e\n', ans_ddx_formF - ans_ddx_hand)


fprintf('Difference ddx_F Vr_B: %1.1e\n', ans_ddx_formF - ans_B)
fprintf('Difference ddx_F Vr_C: %1.1e\n', ans_ddx_formF - ans_C)
fprintf('Difference ddx_F Vr_D: %1.1e\n', ans_ddx_formF - ans_D)
fprintf('Difference ddx_F Vr_E: %1.1e\n', ans_ddx_formF - ans_E)

fprintf('--------------------------\n')

fprintf('Difference ddy_hand Vr_B_ddy: %1.1e\n', ans_ddy_hand - ans_B_ddy)


fprintf('--------------------------\n')

fprintf('Difference ddz_hand Vr_B_ddz: %1.1e\n', ans_ddz_hand - ans_B_ddz)


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
fprintf('\nElapsed time: %1.4f seconds\n',tocWhole)
















