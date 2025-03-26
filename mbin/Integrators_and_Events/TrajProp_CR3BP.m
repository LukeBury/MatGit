% ========================================================================
%%% Description
% ========================================================================
% Propagate a trajectory with zonal harmonics and C22s

% Created: 3/27/22
% Author : Luke Bury, luke.bury@colorado.edu
% ========================================================================
%%% Initialization
% ========================================================================
clear
clc
close all
% mbinPath = '~/CU_Google_Drive/Documents/MatGit/mbin';
mbinPath = '~/Documents/MatGit/mbin';
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
use_scr_integration = false;

% ========================================================================
%%% Setup
% ========================================================================
% -------------------------------------------------
%%% PO IC
% -------------------------------------------------

% PO_IC = [0.9996165286614656;
%  0.0000000000000000;
%  -0.0009890792364341;
%  0.0000000000000000;
%  -0.0161973279788486;
%  0.0000000000000000;
%  4.8871289984803123];

% PO_IC = [0.9996165392231373;
%  0.0000000000000000;
%  -0.0009890876428918;
%  0.0000000000000000;
%  -0.0161901083758854;
%  0.0000000000000000;
%  4.8873805631699803];

% PO_IC = [1.0011102338259017;
%  0.0000000000000000;
%  -0.0002847905318245;
%  0.0000000000000000;
%  0.0138174040670202;
%  0.0000000000000000;
%  5.4397139997492401];
% 
% PO_IC = [1.0022677239299949;
%  0.0000000000000000;
%  -0.0031004095147306;
%  0.0000000000000000;
%  0.0031645135161638;
%  0.0000000000000000;
%  10.5724707789779533];
% 
% PO_IC = [1.0022480798114766;
%  0.0000000000000000;
%  -0.0030691798507714;
%  0.0000000000000000;
%  0.0033148412434564;
%  0.0000000000000000;
%  32.1080587744261834];
% 
% 
% PO_IC = [1.0005168282662340; % DPO_2P4_2T 10 km
%  0.0000000000000000;
%  -0.0009744895128070;
%  0.0000000000000020;
%  0.0158061372098970;
%  0.0000000000000020;
%  10.6103542894260929];

% PO_IC = [1.0022345263624222; % DPO_2P4_2T_1P3 10 km
%  0.0000000000000000;
%  -0.0030794633431889;
%  0.0000000000000000;
%  0.0033359547426397;
%  0.0000000000000000;
%  32.1839218574248633];

% PO_IC = [1.0020051348575794;
%  0.0000000000000000;
%  0.0037549529509182;
%  0.0000000000000000;
%  0.0007615935137484;
%  0.0000000000000000;
%  14.2080356367757474];
% 
% PO_IC = [1.0020363204721270;
%  0.0000000000000000;
%  0.0037632829773217;
%  0.0000000000000000;
%  0.0006276306774677;
%  0.0000000000000000;
%  14.2650184695806903];

% PO_IC = [1.0020133254143000;
%  0.0000000000017137;
%  0.0037568106645151;
%  0.0000000000060406;
%  0.0007292516743476;
%  -0.0000000000304203;
%  14.2238306725398758];
% 
% PO_IC = [1.0000030479714523;
%  0.0002766826179322;
%  -0.0010357660860857;
%  -0.0163225022362263;
%  -0.0002137478615124;
%  -0.0001081231515967;
%  14.2238306725398758];
% 
% PO_IC = [1.0024543617692410; % Decent but unstable one from DRO_1P4_1P3
%  0.0000000000000000;
%  0.0000000000000001;
%  0.0000000000000000;
%  -0.0007926135196438;
%  0.0090394467011972;
%  17.7545986137019227];
% 
% PO_IC = [1.0026187383553777;
%  0.0000000000000000;
%  0.0037505782393589;
%  0.0000000000000000;
%  -0.0014121016688949;
%  0.0000000000000000;
%  22.5986649184966311];
% 
% 
% PO_IC = [0.9878819009166636;
%  0.0000000000000000;
%  -0.0166667481158792;
%  0.0000000000000000;
%  -0.0065998480969078;
%  0.0000000000000000;
%  4.5359646720470677];
% 
% PO_IC = [1.0020468619049392;
%  0.0000000000000000;
%  0.0037488696082717;
%  0.0000000000000000;
%  0.0006874546202919;
%  0.0000000000000000;
%  3.5199734394120314];

% X_BCI = [[0.254529859427125, -1.553560284878574, -2.340877184285415].*1e2, 0.108414229760154, 0.160662289823849, -0.094837907591893]
% [X_BCR] = X_BCI2BCR(X_BCI, 0, bodies.enceladus.meanMot );
% PO_IC =    1.0e+02 .*[0.254529859427125
%                   -1.553560284878574
%                   -2.340877184285415
%                    0.001084142297602
%                    0.001606622898238
%                   -0.000948379075919];
% PO_IC(1:3) = PO_IC(1:3)./rNorm;
% PO_IC(4:6) = PO_IC(4:6)./vNorm;
% 
% PO_IC(1) = PO_IC(1) + (1-prms.u)


% ========================================================================
%%% Integrate and plot the POs
% ========================================================================
% -------------------------------------------------
%%% Set up parameters
% -------------------------------------------------
%%% Set primary and secondary bodies
% [primary, secondary] = assignPrimaryAndSecondary_CR3BP('Jupiter_Europa.CR3BP', bodies)
% [primary, secondary] = assignPrimaryAndSecondary_CR3BP('Jupiter_Io.CR3BP', bodies)
[primary, secondary] = assignPrimaryAndSecondary_CR3BP('Saturn_Enceladus.CR3BP', bodies);

% warning('Manually overriding MR')
% secondary.MR = 1.898884589251784e-07;


%%% Normalizing constants
[rNorm, tNorm, vNorm] = cr3bp_norms(primary, secondary, bodies.constants.G);

%%% prms for integration
prms.u  = secondary.MR;
prms.R2 = secondary.R_n;
prms.R1 = primary.R / rNorm;








test = [1.0e+02 .*[0.572017196151696
0.000026841342967
-2.762397025517204
-0.000006654903285
0.002127140318234
-0.000001357380485]
1.478405366018706];


rNorm_temp = 2.384160591811723E+05;
Tp_temp = 2*pi*sqrt((rNorm_temp^3) / (7.210497553340731 + 3.793120723493890E+07));
tNorm_temp =  Tp_temp / (2*pi);
vNorm_temp = rNorm_temp / tNorm_temp;
test(1:3) = test(1:3)./ rNorm_temp;
test(4:6) = test(4:6)./ vNorm_temp;
test(7) = (test(7)*24*60*60)/tNorm_temp;

PO_IC = test;
PO_IC(1) = PO_IC(1) + (1-prms.u);
PO_IC(7) = PO_IC(7) + 0.05;
prettyColVec(PO_IC)
% prms.J2p  = primary.J2;
% prms.J4p  = primary.J4;
% prms.J6p  = primary.J6;
% prms.J2s  = secondary.J2;
% prms.C22s = secondary.C22;
% prms.J2p  = primary.J2;
% prms.J4p  = primary.J4;
% prms.J6p  = primary.J6;
% prms.J2s  = secondary.J2;
% prms.C22s = secondary.C22;
% prms.C22s = 0;


%%% Determine the mean motion via the ephemeris method
% tN_ephemeris = sqrt((secondary.a^3) / (bodies.constants.G*(primary.mass+secondary.mass)));
% prms.n = secondary.meanMot*tN_ephemeris;
prms.n = 1;
%%% Set integrator handle
integratorHandle = @Int_CR3BnSTM;
% integratorHandle = @Int_CR3BnSTM_J2pJ4pJ6pJ2s_C22s;

%%% Collinear equillibrium points
% rLPs_n = EquilibriumPoints(prms.u, prms.n);
rLPs_n = collinearEquilibriumPoints_ZH_C22s(prms);


if use_scr_integration
    integratorHandle = @Int_CR3BnSTM_SCR;
%     PO_IC(1) = PO_IC(1) - (1-prms.u);
%     rLPs_n(:,1) = rLPs_n(:,1) - (1-prms.u);
end

%%% Integration options
tol = 1e-13;
options = odeset('RelTol',tol,'AbsTol',tol);
options_apsis = odeset('Events', @event_Apsis_CR3BP, 'RelTol',tol,'AbsTol',tol);

%%% Initialize STM as column vector to be added with state
stm0_colVec = reshape(eye(6),36,1);

% [T_out, X_out, t_aps, x_aps, in_aps] = ode113(@Int_CR3BnSTM, [0, PO_IC(end)], [PO_IC(1:6); stm0_colVec], options_apsis, prms);
[T_out, X_out, t_aps, x_aps, in_aps] = ode113(integratorHandle, [0, PO_IC(end)], [PO_IC(1:6); stm0_colVec], options_apsis, prms);

if use_scr_integration
    X_out(:,1) = X_out(:,1) + (1-prms.u);
    x_aps(:,1) = x_aps(:,1) + (1-prms.u);
end

fprintf('Apses:\n')
for kk = 1:size(x_aps,1)
    fprintf('[%1.15f, %1.15f, %1.15f, %1.15f, %1.15f, %1.15f, %1.15f]\n', x_aps(kk,1:6), PO_IC(end))
end

figure(80); hold all
plot3(X_out(:,1),X_out(:,2),X_out(:,3),'linewidth', 2, 'color', colors.mag)
% plot3(x_aps(:,1),x_aps(:,2),x_aps(:,3),'.', 'markersize', 16, 'color', colors.black)
PlotBoi3_CR3Bn(28)
% plotSecondary(secondary)
axis equal

stm_tf_t0                           = reshape(X_out(end,7:42),6,6);
monodromy                           = stm_tf_t0;
[eigenVectors_new, eigenValues_new] = eig(monodromy);
[S1, S2]                            = getStabilityIndices(diag(eigenValues_new));


fprintf('Prop Error = %1.2e\n\n',norm(X_out(end,1:6)' - X_out(1,1:6)') / norm(X_out(1,1:6)));

% 
% 
aps_alts = rowNorm(x_aps(:,1:3) - [1-prms.u, 0, 0]) - prms.R2;
aps_alts_km = aps_alts.*rNorm;

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
















