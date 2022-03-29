% ========================================================================
%%% Description
% ========================================================================
% 

% Created: 12/14/20
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
shooter_allFree          = false;
shooter_TpFixed          = false;
shooter_x0y0Fixed        = false;
shooter_y0Fixed          = false;
shooter_y0xd0Fixed       = false;
shooter_y0zd0Fixed       = false;
shooter_y0xd0zd0Fixed    = false;
shooter_z0Fixed          = true; 
shooter_zd0Fixed         = false;
shooter_y0Fixed_targetTp = false;

shooter_targetJC         = false;
shooter_targetJC_y0Fixed = false;
%     JC_des               = 3.003595836097908; % Europa, VL2 = 50 mps
%     JC_des               = 3.003556098395709; % Europa, VL2 = 100 mps
%     JC_des               = 3.003489868892044; % Europa, VL2 = 150 mps
    JC_des               = 3.003132229572251; % Europa, VL2 = 300 mps


% ========================================================================
%%% Setup
% ========================================================================
% -------------------------------------------------
%%% Initial Guess
% -------------------------------------------------
% famName = 'Jupiter_Europa';
% famName = 'Saturn_Titan';
famName = 'Saturn_Enceladus';
% famName = 'Earth_Moon';
% 989
[primary, secondary] = assignPrimaryAndSecondary_CR3BP(famName, bodies);
[rNorm, tNorm, vNorm] = cr3bp_norms(primary, secondary, bodies.constants.G);




%%% 
PO_guess = [0.9990598425656854;
 0.0000000000000000;
 -0.0028348610498525;
 -0.0000000000518546;
 -0.0071173508131177;
 0.0000000001688107;
 6.0800225227923965];

% -------------------------------------------------
%%% Shooter-specific setup
% -------------------------------------------------
if shooter_allFree
    
elseif shooter_y0Fixed

elseif shooter_x0y0Fixed

elseif shooter_y0xd0Fixed
%     if PO_guess(4) ~= 0
%         warning('Switching xd0 to 0')
%         PO_guess(4) = 0;
%     end
elseif shooter_y0zd0Fixed
    if PO_guess(6) ~= 0
%         warning('Switching zd0 to 0')
%         PO_guess(6) = 0;
    end
elseif shooter_y0xd0zd0Fixed
%     if PO_guess(4) ~= 0
%         warning('Switching xd0 to 0')
%         PO_guess(4) = 0;
%     end
%     if PO_guess(6) ~= 0
%         warning('Switching zd0 to 0')
%         PO_guess(6) = 0;
%     end
elseif shooter_z0Fixed
    
elseif shooter_zd0Fixed
    
elseif shooter_y0Fixed_targetTp
    Tdes = PO_guess(end);
end

% -------------------------------------------------
%%% Shooter options
% -------------------------------------------------
%%% Error tolerance for constraint vector norm in multiple shooter
% error_tol = 1e-8; 
error_tol = 1e-9;  % Usually start here when working from a guess
% error_tol = 1e-10; 
error_tol = 1e-11; 
% error_tol = 1e-12;  %
% error_tol = 5e-13; 
% error_tol = 1e-13; 

%%% Number of nodes for multiple shooter. Generally, higher for bigger POs,
%%% but also be aware that more nodes increase inherent error, so you may
%%% need to lower error tolerances for the shooter
% n_Nodes = 1; 
% n_Nodes = 2; 
% n_Nodes = 3;
% n_Nodes = 4;
n_Nodes = 5;
% n_Nodes = 6;
% n_Nodes = 7; 
% n_Nodes = 8;
% n_Nodes = 9;
% n_Nodes = 10;
% n_Nodes = 11;
% n_Nodes = 12;
% n_Nodes = 13;
% n_Nodes = 14;
% n_Nodes = 16;
% n_Nodes = 20;

%%% Maximum number of multiple-shooter iterations for any given family
%%% member before kicking out of the loop and adjusting the tuning
%%% parameters
ms_iterMax = 1200;

%%% Artificial scalar for increasing or decreasing the step size for 
% guessing at the state of the 2nd family member. 1 has no effect, 2 and 
% 3 can be helpful in certain regions
stepSize = 1;

% -------------------------------------------------
%%% System
% -------------------------------------------------


%%% Setting parameters
prms.u  = secondary.MR;
prms.n  = 1;
prms.R2 = secondary.R_n;

%%% Collinear equillibrium points
rLPs_n = EquilibriumPoints(prms.u, prms.n);

% ------------------------------------------------- 
%%% Integration Options
% -------------------------------------------------
tol = 2.22045e-14;
options = odeset('RelTol',tol,'AbsTol',tol);


% ========================================================================
%%% Correct guess with shooter
% ========================================================================
% -------------------------------------------------
%%% Correct the guess using desired shooting method
% -------------------------------------------------
if shooter_allFree
    [PO_result, counter, constraint_error] = correctPO_multShooter_stateContinuity(PO_guess, n_Nodes, @Int_CR3BnSTM, options, prms, error_tol, ms_iterMax, stepSize);
elseif shooter_TpFixed
    [PO_result, counter, constraint_error] = correctPO_multShooter_stateContinuity_TpFixed(PO_guess, n_Nodes, @Int_CR3BnSTM, options, prms, error_tol, ms_iterMax, stepSize, PO_guess(7));
elseif shooter_x0y0Fixed
    [PO_result, counter, constraint_error] = correctPO_multShooter_stateContinuity_x0y0Fixed(PO_guess, n_Nodes, @Int_CR3BnSTM, options, prms, error_tol, ms_iterMax, stepSize, PO_guess(1), PO_guess(2));
elseif shooter_y0Fixed
    [PO_result, counter, constraint_error] = correctPO_multShooter_stateContinuity_y0Fixed(PO_guess, n_Nodes, @Int_CR3BnSTM, options, prms, error_tol, ms_iterMax, stepSize, PO_guess(2));
elseif shooter_y0xd0Fixed
    [PO_result, counter, constraint_error] = correctPO_multShooter_stateContinuity_y0xd0Fixed(PO_guess, n_Nodes, @Int_CR3BnSTM, options, prms, error_tol, ms_iterMax, stepSize, PO_guess(2), PO_guess(4));
elseif shooter_y0zd0Fixed
    [PO_result, counter, constraint_error] = correctPO_multShooter_stateContinuity_y0zd0Fixed(PO_guess, n_Nodes, @Int_CR3BnSTM, options, prms, error_tol, ms_iterMax, stepSize, PO_guess(2), PO_guess(6));
elseif shooter_y0xd0zd0Fixed
    [PO_result, counter, constraint_error] = correctPO_multShooter_stateContinuity_y0xd0zd0Fixed(PO_guess, n_Nodes, @Int_CR3BnSTM, options, prms, error_tol, ms_iterMax, stepSize, PO_guess(2), PO_guess(4), PO_guess(6));
elseif shooter_z0Fixed
    [PO_result, counter, constraint_error] = correctPO_multShooter_stateContinuity_z0Fixed(PO_guess, n_Nodes, @Int_CR3BnSTM, options, prms, error_tol, ms_iterMax, stepSize);
elseif shooter_zd0Fixed
    [PO_result, counter, constraint_error] = correctPO_multShooter_stateContinuity_zd0Fixed(PO_guess, n_Nodes, @Int_CR3BnSTM, options, prms, error_tol, ms_iterMax, stepSize);
elseif shooter_y0Fixed_targetTp
    [PO_result, counter, constraint_error] = correctPO_multShooter_stateContinuity_y0Fixed_targetTp(PO_guess, n_Nodes, @Int_CR3BnSTM, options, prms, error_tol, ms_iterMax, stepSize, Tdes);
elseif shooter_targetJC
    [PO_result, counter, constraint_error] = correctPO_multShooter_stateContinuity_JCFixed(PO_guess, n_Nodes, @Int_CR3BnSTM, options, prms, error_tol, ms_iterMax, stepSize, JC_des);
elseif shooter_targetJC_y0Fixed
    [PO_result, counter, constraint_error] = correctPO_mS_sC_JCy0Fixed(PO_guess, n_Nodes, @Int_CR3BnSTM, options, prms, error_tol, ms_iterMax, stepSize, JC_des, PO_guess(2));
end
% -------------------------------------------------
%%% Integrate and plot the results
% -------------------------------------------------
[T_guess, X_guess] = ode113(@Int_CR3BnSTM, linspace(0, PO_guess(7), 50000), [PO_guess(1:6); reshape(eye(6),36,1)], options, prms);
[T_result, X_result] = ode113(@Int_CR3BnSTM, linspace(0, PO_result(7), 50000), [PO_result(1:6); reshape(eye(6),36,1)], options, prms);


figure('position',[43 423 560 420]); hold all
p_guess  = plot3(X_guess(:,1),X_guess(:,2),X_guess(:,3),'b');
PlotBoi3_CR3Bn(26)
legend([p_guess], 'guess');

figure('position',[604 423 560 420]); hold all
p_guess  = plot3(X_guess(:,1),X_guess(:,2),X_guess(:,3),'b');
p_result = plot3(X_result(:,1),X_result(:,2),X_result(:,3),'m');
PlotBoi3_CR3Bn(26)
legend([p_guess p_result], 'guess', 'result');

figure('position', [1165 423 560 420]); hold all
p_result = plot3(X_result(:,1),X_result(:,2),X_result(:,3),'m');
plot3(X_result(1,1),X_result(1,2),X_result(1,3),'ko')
plot3(X_result(end,1),X_result(end,2),X_result(end,3),'kx')
PlotBoi3_CR3Bn(26)
legend([p_result],'result');

%%% Find stability info
stm_tf_t0                           = reshape(X_result(end,7:42),6,6);
monodromy                           = stm_tf_t0;
[eigenVectors_new, eigenValues_new] = eig(monodromy);
[S1, S2]                            = getStabilityIndices(diag(eigenValues_new));

fprintf('PO_0:\n')
prettyColVec(PO_result)

fprintf('Stability Indices: [%1.1f, %1.1f]\n', S1, S2)

propError = norm(X_result(end,1:6)' - X_result(1,1:6)') / norm(X_result(1,1:6));

fprintf('iterations: %1.3d\nconstraint error = %1.2e\nprop error = %1.2e\n', counter, constraint_error, propError)

%%% Use this to choose a y=0 crossing. Ideally looking for one
%%% where xd and zd are also 0
options_yEquals0_nonTerminal = odeset('Event',@event_yEqualsZero_nonTerminal,'RelTol',tol,'AbsTol',tol);
[Tfix_fwd_non, Xfix_fwd_non, tEv_fwd_non, xEv_fwd_non, ~] = ode113(@Int_CR3BnSTM, [0, PO_result(7)], [PO_result(1:6); reshape(eye(6),36,1)], options_yEquals0_nonTerminal, prms);
fprintf('y-axis crossings:\n')
for kk = 1:size(xEv_fwd_non,1)
    fprintf('[%1.15f, %1.15f, %1.15f, %1.15f, %1.15f, %1.15f]''\n', xEv_fwd_non(kk,1:6))
end

% % X_test = [1.0206462604042468; % From Hg1 - the nearest to going inside L1/L2
% %  0.0000000000000000;
% %  -0.0000000000000000;
% %  0.0000000000000000;
% %  -0.0258834005296429;
% %  0.0000000000000000;
% %  5.4390036058067066];
% % 
% % 
% % X_test = [1.0189157565954854; % from Hg2
% %  0.0000000000000000;
% %  0.0000000000000000;
% %  -0.0000000000000;
% %  -0.0026288580808985;
% %  0.0000000000000000;
% %  8.2086350010682185];
% % 
% % [1.0196292235860616; % from Hg4
% %  0.0000000000000000;
% %  -0.0000000000000000;
% %  0.0000000000000000;
% %  0.0014894682946456;
% %  0.0000000000000000;
% %  13.4305140203815032];
% % 
% % % % % % % 
% % % % % X_test = [1.0209731634679116;
% % % % %  0.0000000000000000;
% % % % %  0.0000000000000000;
% % % % %  0.0000000000000013;
% % % % %  -0.0054460597115365;
% % % % %  0.0000000000000000; 
% % % % %  11.4890325492471597];

% % 
% X_test = [1.0093690704525875;
%  -0.0;
%  -0.0000000000000000;
%  -0.000;
%  -0.0626;
%  0.0000000000000000;
%  1.01];
% 
% % % % 
% [T_test_out, X_test_out] = ode113(@Int_CR3BnSTM, linspace(0, X_test(7), 10000), [X_test(1:6); reshape(eye(6),36,1)], options, prms);
% figure; hold all
% plot3(X_test_out(:,1)-(1-prms.u),X_test_out(:,2),X_test_out(:,3),'k')
% plot3(rLPs_n(1:2,1)-(1-prms.u), [0,0], [0,0], 'b^')
% PlotBoi3_CR3Bn(26)


% difference of 1.8e-4


%%% 
% L2 = 1.020461385981616;
% X_test = [1.01;
%  0;
%  0;
%  0;
%  0.932;
%  0;
%  2*pi];
% X_test = [1.01;
%  0;
%  0;
%  0;
%  -1.18;
%  0;
%  2*pi];
% X_test = [0.993;
%  0;
%  0;
%  0;
%  -1.8;
%  0;
%  4*pi];
% [T_test_out, X_test_out] = ode113(@Int_CR3BnSTM, linspace(0, X_test(7), 10000), [X_test(1:6); reshape(eye(6),36,1)], options, prms);
% figure; hold all
% plot3(X_test_out(:,1),X_test_out(:,2),X_test_out(:,3),'k')
% plot3(rLPs_n(1:2,1), [0,0], [0,0], 'b^')
% PlotBoi3_CR3Bn(26)
% plotSecondary(secondary)
% plotPrimary(primary, secondary)
% view(0, 90)



% X_test = [1.0189157565954854; %% Looks like a weird Hg1
%  0.0000000000000000;
%  0.0000000000000000;
%  -0.0000000000000;
%  -0.04;
%  0.0000000000000000;
% 5.5];



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















