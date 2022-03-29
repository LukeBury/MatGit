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

% ========================================================================
%%% Run Switches
% ========================================================================
chooseCustomPO = true;

% ========================================================================
%%% Setup
% ========================================================================
% -------------------------------------------------
%%% Shooter options
% -------------------------------------------------
% -----------------------
%%% Which shooter to run? true or false
% -----------------------
shooter_allFree         = false;
shooter_y0Fixed         = false;
shooter_y0xd0Fixed      = false;    
shooter_y0zd0Fixed      = false;
shooter_y0xd0zd0Fixed   = false;
shooter_z0Fixed         = false;
shooter_zd0Fixed        = false;

shooter_TpFixed         = true;
shooter_y0TpFixed       = false; 
shooter_y0xd0TpFixed    = false; 
shooter_y0xd0zd0TpFixed = false;
%     extraTpScale = 1;
    extraTpScale = 1 + 5e-5;

% -----------------------
%%% Choose how to select the PO to continue from a family. 
% -----------------------
choosePO_middle         = true;  % Middle member on Tp-vs-JC curve
choosePO_mostStable     = false; % Most stable member

if choosePO_middle && choosePO_mostStable
    warning('Only choose one')
    return
end

% -----------------------
%%% General shooting options
% -----------------------
%%% Number of nodes for shooter
% n_Nodes = 3;
% n_Nodes = 4;
% n_Nodes = 5;
% n_Nodes = 6;
% n_Nodes = 7;
n_Nodes = 8;
% n_Nodes = 9; %
% n_Nodes = 10;
% n_Nodes = 11;
% n_Nodes = 12;
% n_Nodes = 13;
% n_Nodes = 14;
% n_Nodes = 15;
% n_Nodes = 16;

%%% Error tolerance for constraint vector norm in multiple shooter
% error_tol = 1e-9; 989
% error_tol = 1e-10; 989
error_tol = 1e-11; 
% error_tol = 1e-12; 
% error_tol = 5e-13; 
% error_tol = 1e-13; %
% 
%%% Maximum number of iterations for shooter
ms_iterMax = 1000;

%%% I should really get rid of this as a parameter
stepSize = 1;

% -------------------------------------------------
%%% Options for continuation
% -------------------------------------------------
%%% Choose number of steps taken to get from unperturbed CR3BP to perturbed
% n_steps = 50; 
% n_steps = 100; 
% n_steps = 200; 
% n_steps = 300; 989
n_steps = 400; 989
% n_steps = 600; 989
% n_steps = 1000; 989
% n_steps = 100000; 989


% -------------------------------------------------
%%% Choose PO data for continuation
% -------------------------------------------------
if chooseCustomPO
    %%% Set family name
    family = 'Jupiter_Europa.CR3BP'; 
    
%     %%% Set [X0; Tp] of initial unperturbed periodic orbit 
%     PO_IC = [1.0162778160411139;
%          0.0000000000000000;
%          -0.0045663681005709;
%          -0.0058921216560423;
%          -0.0518126082370704;
%          0.0110421781965304;
%          25.1422570610823612];
     
     
PO_IC = [1.0177368259875725;
 0.0000000000000000;
 -0.0000000000002130;
 0.0000000000003213;
 0.0049616210115872;
 -0.0057972982040939;
 23.0711171494392211];

    
% PO_IC = [1.0179169125744147;
%  0.0012288820527487;
%  -0.0000748630089328;
%  -0.0011099415331474;
%  0.0040409868699154;
%  -0.0057058898863069;
%  23.0830178398857484];
     
    
else
%     family = 'Jupiter_Europa.CR3BP.L2_Lyapunov.txt';
%     family = 'Jupiter_Europa.CR3BP.L2_NHalo.txt';
%     family = 'Jupiter_Europa.CR3BP.L2_SHalo.txt';
%     family = 'Jupiter_Europa.CR3BP.L2_Vertical.txt';
%     family = 'Jupiter_Europa.CR3BP.L2_L_1T_1P4.txt'; 
%     family = 'Jupiter_Europa.CR3BP.L2_L_1T_1P3.txt'; 
%     family = 'Jupiter_Europa.CR3BP.L2_L_1T_1P2.txt'; 
%     family = 'Jupiter_Europa.CR3BP.L2_L_1T_2P2.txt'; 
%     family = 'Jupiter_Europa.CR3BP.L2_L_1T_3P2.txt'; 
%     family = 'Jupiter_Europa.CR3BP.L2_L_1T_4P2.txt'; % Butterly
%     family = 'Jupiter_Europa.CR3BP.L2_L_1T_4P2_1P3.txt'; % the end of this is the 2P2 bifurcation from LoPO_2P2_1T (so this is also LoPO_2P2_1T_2P2)
%     family = 'Jupiter_Europa.CR3BP.L2_L_1T_4P2_1P2.txt'; % 
%     family = 'Jupiter_Europa.CR3BP.L2_L_1T_4P2_2P2.txt'; % 
%     family = 'Jupiter_Europa.CR3BP.L2_L_1T_4P2_2P3.txt'; % 
%     family = 'Jupiter_Europa.CR3BP.L2_L_1T_4P2_2P4.txt'; % 
%     family = 'Jupiter_Europa.CR3BP.L2_L_2T.txt'; % Full Axial family (Also L2_V_T)
%     family = 'Jupiter_Europa.CR3BP.L2_L_1P3.txt'; 
%     family = 'Jupiter_Europa.CR3BP.L2_L_1P2.txt'; 

%     family = 'Jupiter_Europa.CR3BP.DRO.txt'; 
%     family = 'Jupiter_Europa.CR3BP.DRO_1P4.txt'; 
%     family = 'Jupiter_Europa.CR3BP.DRO_1P4_1P4.txt'; 
%     family = 'Jupiter_Europa.CR3BP.DRO_1P4_1P4_1T.txt'; 
%     family = 'Jupiter_Europa.CR3BP.DRO_1P3.txt'; 
%     family = 'Jupiter_Europa.CR3BP.DRO_1P3_1P3.txt'; 
%     family = 'Jupiter_Europa.CR3BP.DRO_2P3.txt'; 
%     family = 'Jupiter_Europa.CR3BP.DRO_2P3_1P2.txt'; 
%     family = 'Jupiter_Europa.CR3BP.DRO_2P3_2P2.txt'; 
%     family = 'Jupiter_Europa.CR3BP.DRO_2P4.txt'; 
%     family = 'Jupiter_Europa.CR3BP.DRO_2P4_1T.txt'; 
%     family = 'Jupiter_Europa.CR3BP.DRO_2P4_2T.txt'; 
% 
%     family = 'Jupiter_Europa.CR3BP.DPO.txt'; 
%     family = 'Jupiter_Europa.CR3BP.DPO_1P4.txt'; 
%     family = 'Jupiter_Europa.CR3BP.DPO_1P4_1T.txt'; 
%     family = 'Jupiter_Europa.CR3BP.DPO_1P3.txt'; 
%     family = 'Jupiter_Europa.CR3BP.DPO_1P2.txt'; 
%     family = 'Jupiter_Europa.CR3BP.DPO_2P3.txt'; 
%     family = 'Jupiter_Europa.CR3BP.DPO_2P2.txt'; 
%     family = 'Jupiter_Europa.CR3BP.DPO_3P2.txt'; 
%     family = 'Jupiter_Europa.CR3BP.DPO_2P4.txt'; 
%     family = 'Jupiter_Europa.CR3BP.DPO_3P3.txt'; 
%     family = 'Jupiter_Europa.CR3BP.DPO_3P4.txt'; 
% % %     family = 'Jupiter_Europa.CR3BP.DPO_3P4_symmetricAboutXb2.txt'; 
%     family = 'Jupiter_Europa.CR3BP.DPO_4P4.txt'; 
%     family = 'Jupiter_Europa.CR3BP.DPO_4P3.txt'; 
%     family = 'Jupiter_Europa.CR3BP.DPO_5P2.txt'; 
%     family = 'Jupiter_Europa.CR3BP.DPO_2T.txt'; 
% 
%     family = 'Jupiter_Europa.CR3BP.LoPO.txt'; 
%     family = 'Jupiter_Europa.CR3BP.LoPO_1P4.txt'; 
%     family = 'Jupiter_Europa.CR3BP.LoPO_1P3.txt'; 
%     family = 'Jupiter_Europa.CR3BP.LoPO_1P2.txt'; 
%     family = 'Jupiter_Europa.CR3BP.LoPO_2P2.txt'; 
    family = 'Jupiter_Europa.CR3BP.LoPO_2P2_1T.txt'; 
%     family = 'Jupiter_Europa.CR3BP.LoPO_3P2.txt'; 
%     family = 'Jupiter_Europa.CR3BP.LoPO_2P4.txt'; 
%     family = 'Jupiter_Europa.CR3BP.LoPO_3P3.txt'; 
%     family = 'Jupiter_Europa.CR3BP.LoPO_3P4.txt'; 

%     family = 'Jupiter_Europa.CR3BP.Hg2.txt'; 
%     family = 'Jupiter_Europa.CR3BP.Hg2_1P2.txt'; 
%     family = 'Jupiter_Europa.CR3BP.Hg2_2P2.txt'; 
%     family = 'Jupiter_Europa.CR3BP.Hg2_1P3.txt'; 
%     family = 'Jupiter_Europa.CR3BP.Hg2_2T.txt'; 
%     family = 'Jupiter_Europa.CR3BP.Hg2_3T.txt'; 
%     family = 'Jupiter_Europa.CR3BP.Hg2_5T.txt';
% 
%     family = 'Jupiter_Europa.CR3BP.Se7.txt'; 
%     family = 'Jupiter_Europa.CR3BP.Se7_2P2.txt'; 
%     family = 'Jupiter_Europa.CR3BP.Se7_2T.txt'; 
%     family = 'Jupiter_Europa.CR3BP.Se7_3T.txt'; 
%     family = 'Jupiter_Europa.CR3BP.Se7_5P2.txt'; 
%     family = 'Jupiter_Europa.CR3BP.Se7_6P2.txt'; 
%     family = 'Jupiter_Europa.CR3BP.Se7_5T.txt'; 
%     family = 'Jupiter_Europa.CR3BP.Se7_6T.txt'; 
% 
%     family = 'Jupiter_Europa.CR3BP.Hg1.txt'; 
%     family = 'Jupiter_Europa.CR3BP.Hg1_1P2.txt'; 
%     family = 'Jupiter_Europa.CR3BP.Hg1_2P2.txt'; 
%     family = 'Jupiter_Europa.CR3BP.Hg1_1P3.txt'; 
%     family = 'Jupiter_Europa.CR3BP.Hg1_1T.txt'; 
%     family = 'Jupiter_Europa.CR3BP.Hg1_2T.txt'; 

%     family = 'Jupiter_Europa.CR3BP.Hg3.txt'; 

    % --------------------------
    % Actually load family
    % --------------------------
    %%% Path from mbin to data
    dataPathFromMBin = '/Data/InitialConditions/PO_Families/';

    %%% PO data file
    PO_datafile = [mbinPath, dataPathFromMBin, family];
    
    %%% Load the data file
    PO_data = dlmread(PO_datafile,',',1,0);

    %%% Grab header line
    fid = fopen(PO_datafile, 'rt');  %the 't' is important!
    header = fgetl(fid);
    fclose(fid);
    
    % --------------------------
    % Create column specifiers
    % --------------------------
    PO_header_2020 = 'x0,y0,z0,xd0,yd0,zd0,Tp,JC,stabilityIndex1,stabilityIndex2,alpha,beta,impactFlag,error';

    if contains(header, PO_header_2020)
        c_x0 = 1;   c_y0 = 2;   c_z0 = 3;
        c_xd0 = 4;  x_yd0 = 5;  c_zd0 = 6;
        c_Tp = 7;   c_JC = 8;   c_S1 = 9;   c_S2 = 10;
        c_alpha = 11;   c_beta = 12;    c_impactFlag = 13;
        c_error = 14;
    end
    
    
    plot_PO_indices = getIndices_spacedByTpJcArclength(PO_data(:,c_Tp), PO_data(:,c_JC), 3);
    plot_PO_indices = plot_PO_indices(2);
    plot_PO_indices = 400; 989 % Overriding
    PO_IC = PO_data(plot_PO_indices, c_x0:c_Tp)';
    
%     newVec = PO_data(:,c_S1) + PO_data(:,c_S2);
%     plot_PO_indices = find(newVec == min(newVec));
%     plot_PO_indices = plot_PO_indices(1);
%     PO_IC = PO_data(plot_PO_indices, c_x0:c_Tp)';
    
    fprintf('Family: %s\n', family)
    fprintf('Initial PO Index: %1d\n', plot_PO_indices)
    fprintf('Stability Indices: [%1.3f, %1.3f]\n', PO_data(plot_PO_indices, c_S1), PO_data(plot_PO_indices, c_S2));
    prettyColVec(PO_IC)
           
end

% 989
% return

% -------------------------------------------------
%%% Parameter setup
% -------------------------------------------------
%%% Set primary and secondary bodies
[primary, secondary] = assignPrimaryAndSecondary_CR3BP(family, bodies);

%%% prms for integration
prms.u  = secondary.MR;
prms.R2 = secondary.R_n;
prms.R1 = primary.R / secondary.a;
prms.J2p = primary.J2;
prms.J4p = primary.J4;
prms.J6p = primary.J6;
prms.J2s = secondary.J2;

%%% Determine the mean motion via the ephemeris method
tN_ephemeris = sqrt((secondary.a^3) / (bodies.constants.G*(primary.mass+secondary.mass)));
prms.n_ephemeris = secondary.meanMot*tN_ephemeris;

%%% Integration options
% tol = 5e-13; 989
tol = 1e-13;
options = odeset('RelTol',tol,'AbsTol',tol);

%%% Initialize STM as column vector to be added with state
stm0_colVec = reshape(eye(6),36,1);

%%% Plot initial PO
prms.n = 1;
[~, X_PO_IC] = ode113(@Int_CR3BnSTM, [0, PO_IC(7)], [PO_IC(1:6); stm0_colVec], options, prms);
figure('position', [1179 597 560 420]); hold all;
plot3(X_PO_IC(:,1), X_PO_IC(:,2), X_PO_IC(:,3), 'linewidth', 2, 'color', colors.blue2)
PlotBoi3_CR3Bn(26)
if (~isequal(PO_IC(3), 0)) ||  (~isequal(PO_IC(6), 0))
    view(0,0)
end

%%% Stability info for PO_IC
stm_tf_t0                           = reshape(X_PO_IC(end,7:42),6,6);
monodromy                           = stm_tf_t0;
[eigenVectors_new, eigenValues_new] = eig(monodromy);
[S1, S2]                            = getStabilityIndices(diag(eigenValues_new));

%%% Perturbation scaling vector
perturbationScale_vec = linspace(0, 1, n_steps+1);
perturbationScale_vec = perturbationScale_vec(2:end); 

%%% Initialize guess
PO_guess = PO_IC;

%%% Get first fundamental Tp
rLPs_n = EquilibriumPoints(prms.u, 1);
A_mat_Lp = get_Amat_CR3BP(prms.u, rLPs_n(2,:), 1);
[Tp_fundamentalLyap_last] = get_fundamentalLyapunovTp(A_mat_Lp);

%%% optional color spectrum
color_spectrum = colorScale([colors.blue2; colors.mag], n_steps);
figure(10); hold all
if (~isequal(PO_IC(3), 0)) ||  (~isequal(PO_IC(6), 0))
    view(0,0)
end
PlotBoi3_CR3Bn(26)

%%% Preallocate space for PO solutions
PO_solutions = zeros(n_steps, 7);

% 989
% return
for step_i = 1:n_steps
% for step_i = 111:n_steps
    %%% Set current perturbation scaling factor
    prms.pertScale = perturbationScale_vec(step_i);
    
    %%% Set new mean motion for the current step
    prms.n = 1 + prms.pertScale*(prms.n_ephemeris-1);
    
    %%% If we're predicting the Tp
    if shooter_TpFixed || shooter_y0TpFixed || shooter_y0xd0TpFixed || shooter_y0xd0zd0TpFixed
        rLPs_n = collinearEquilibriumPoints_J2pJ4pJ6pJ2s_scaled(prms);
        A_mat_Lp = get_Amat_CR3BP(prms.u, rLPs_n(2,:), prms.n);
        [Tp_fundamentalLyap_current] = get_fundamentalLyapunovTp(A_mat_Lp);
        
        scale_Tp = Tp_fundamentalLyap_current / Tp_fundamentalLyap_last;
        
        PO_guess(7) = PO_guess(7)*scale_Tp*extraTpScale;

    end
    % -------------------------------------------------
    %%% Correct the guess using desired shooting method
    % -------------------------------------------------
    warning('off','MATLAB:nearlySingularMatrix')
    if shooter_allFree
        [PO_result, counter, constraint_error] = correctPO_multShooter_stateContinuity(PO_guess, n_Nodes, @Int_CR3BnSTM_J2pJ4pJ6pJ2s_scaled, options, prms, error_tol, ms_iterMax, stepSize);
    elseif shooter_y0Fixed
        [PO_result, counter, constraint_error] = correctPO_multShooter_stateContinuity_y0Fixed(PO_guess, n_Nodes, @Int_CR3BnSTM_J2pJ4pJ6pJ2s_scaled, options, prms, error_tol, ms_iterMax, stepSize, PO_guess(2));
    elseif shooter_y0xd0Fixed
        [PO_result, counter, constraint_error] = correctPO_multShooter_stateContinuity_y0xd0Fixed(PO_guess, n_Nodes, @Int_CR3BnSTM_J2pJ4pJ6pJ2s_scaled, options, prms, error_tol, ms_iterMax, stepSize, PO_guess(2), PO_guess(4));
    elseif shooter_y0zd0Fixed
        [PO_result, counter, constraint_error] = correctPO_multShooter_stateContinuity_y0zd0Fixed(PO_guess, n_Nodes, @Int_CR3BnSTM_J2pJ4pJ6pJ2s_scaled, options, prms, error_tol, ms_iterMax, stepSize, PO_guess(2), PO_guess(6));
    elseif shooter_y0xd0zd0Fixed
        [PO_result, counter, constraint_error] = correctPO_multShooter_stateContinuity_y0xd0zd0Fixed(PO_guess, n_Nodes, @Int_CR3BnSTM_J2pJ4pJ6pJ2s_scaled, options, prms, error_tol, ms_iterMax, stepSize, PO_guess(2), PO_guess(4), PO_guess(6));
    elseif shooter_z0Fixed
        [PO_result, counter, constraint_error] = correctPO_multShooter_stateContinuity_z0Fixed(PO_guess, n_Nodes, @Int_CR3BnSTM_J2pJ4pJ6pJ2s_scaled, options, prms, error_tol, ms_iterMax, stepSize);
    elseif shooter_zd0Fixed
        [PO_result, counter, constraint_error] = correctPO_multShooter_stateContinuity_zd0Fixed(PO_guess, n_Nodes, @Int_CR3BnSTM_J2pJ4pJ6pJ2s_scaled, options, prms, error_tol, ms_iterMax, stepSize);
    elseif shooter_TpFixed
        [PO_result, counter, constraint_error] = correctPO_multShooter_stateContinuity_TpFixed(PO_guess, n_Nodes, @Int_CR3BnSTM_J2pJ4pJ6pJ2s_scaled, options, prms, error_tol, ms_iterMax, stepSize, PO_guess(7));
    elseif shooter_y0TpFixed
        [PO_result, counter, constraint_error] = correctPO_multShooter_stateContinuity_y0TpFixed(PO_guess, n_Nodes, @Int_CR3BnSTM_J2pJ4pJ6pJ2s_scaled, options, prms, error_tol, ms_iterMax, stepSize, PO_guess(2), PO_guess(7));
    elseif shooter_y0xd0TpFixed
        [PO_result, counter, constraint_error] = correctPO_mS_sC_y0xd0TpFixed(PO_guess, n_Nodes, @Int_CR3BnSTM_J2pJ4pJ6pJ2s_scaled, options, prms, error_tol, ms_iterMax, stepSize, PO_guess(2), PO_guess(4), PO_guess(7));
    elseif shooter_y0xd0zd0TpFixed
        [PO_result, counter, constraint_error] = correctPO_mS_sC_y0xd0zd0TpFixed(PO_guess, n_Nodes, @Int_CR3BnSTM_J2pJ4pJ6pJ2s_scaled, options, prms, error_tol, ms_iterMax, stepSize, PO_guess(2), PO_guess(4), PO_guess(6), PO_guess(7));
    end
    warning('on')
    
    %%% Calculate propagation error
    [T_result, X_result] = ode113(@Int_CR3BnSTM_J2pJ4pJ6pJ2s_scaled, [0, PO_result(7)], [PO_result(1:6); stm0_colVec], options, prms);
    figure(10); plot3(X_result(:,1), X_result(:,2), X_result(:,3), 'linewidth', 1, 'color', color_spectrum(step_i,:))
    propError = norm(X_result(end,1:6)' - X_result(1,1:6)') / norm(X_result(1,1:6));
    
    %%% Store solution
    PO_solutions(step_i,:) = PO_result';

    fprintf('Step = %1d ... iterations: %1.3d ... prop error = %1.2e ... Elapsed time: %1.1f seconds\n',...
        step_i, counter, propError, toc(ticWhole))

    %%% Inialize for next step
    PO_guess = PO_result;
    if shooter_TpFixed || shooter_y0TpFixed || shooter_y0xd0TpFixed || shooter_y0xd0zd0TpFixed
        Tp_fundamentalLyap_last = Tp_fundamentalLyap_current;
    end
    
end

fprintf('PO Result, row:\n')
prettyColVec(PO_result, 'row')


fprintf('PO Result, column:\n')
prettyColVec(PO_result)


figure(1)
plot3(X_result(:,1),X_result(:,2),X_result(:,3), 'linewidth',2, 'color', color_spectrum(end,:))














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























% -------------------------------------------------
%%% Copy-paste script for restarting at some point along the way
% -------------------------------------------------
if 1+1 == 1
    
    close(figure(10))
    figure(10); hold all
    PlotBoi3_CR3Bn(26)
%     view(0,0)

    ticWhole = tic;
    startingGuessIndex = 36;
    n_Nodes = 6;


    PO_guess = PO_solutions(startingGuessIndex,:)';
    % for step_i = 1:n_steps
    for step_i = (startingGuessIndex+1):n_steps
    %%% Set current perturbation scaling factor
    prms.pertScale = perturbationScale_vec(step_i);
    %%% Set new mean motion for the current step
    prms.n = 1 + prms.pertScale*(prms.n_ephemeris-1);
    %%% If we're predicting the Tp
    if shooter_TpFixed || shooter_y0TpFixed || shooter_y0xd0TpFixed || shooter_y0xd0zd0TpFixed
        rLPs_n = collinearEquilibriumPoints_J2pJ4pJ6pJ2s_scaled(prms);
        A_mat_Lp = get_Amat_CR3BP(prms.u, rLPs_n(2,:), prms.n);
        [Tp_fundamentalLyap_current] = get_fundamentalLyapunovTp(A_mat_Lp);
        scale_Tp = Tp_fundamentalLyap_current / Tp_fundamentalLyap_last;
        PO_guess(7) = PO_guess(7)*scale_Tp*extraTpScale;
    end
    % -------------------------------------------------
    %%% Correct the guess using desired shooting method
    % -------------------------------------------------
    warning('off','MATLAB:nearlySingularMatrix')
    if shooter_allFree
        [PO_result, counter, constraint_error] = correctPO_multShooter_stateContinuity(PO_guess, n_Nodes, @Int_CR3BnSTM_J2pJ4pJ6pJ2s_scaled, options, prms, error_tol, ms_iterMax, stepSize);
    elseif shooter_y0Fixed
        [PO_result, counter, constraint_error] = correctPO_multShooter_stateContinuity_y0Fixed(PO_guess, n_Nodes, @Int_CR3BnSTM_J2pJ4pJ6pJ2s_scaled, options, prms, error_tol, ms_iterMax, stepSize, PO_guess(2));
    elseif shooter_y0xd0Fixed
        [PO_result, counter, constraint_error] = correctPO_multShooter_stateContinuity_y0xd0Fixed(PO_guess, n_Nodes, @Int_CR3BnSTM_J2pJ4pJ6pJ2s_scaled, options, prms, error_tol, ms_iterMax, stepSize, PO_guess(2), PO_guess(4));
    elseif shooter_y0zd0Fixed
        [PO_result, counter, constraint_error] = correctPO_multShooter_stateContinuity_y0zd0Fixed(PO_guess, n_Nodes, @Int_CR3BnSTM_J2pJ4pJ6pJ2s_scaled, options, prms, error_tol, ms_iterMax, stepSize, PO_guess(2), PO_guess(6));
    elseif shooter_y0xd0zd0Fixed
        [PO_result, counter, constraint_error] = correctPO_multShooter_stateContinuity_y0xd0zd0Fixed(PO_guess, n_Nodes, @Int_CR3BnSTM_J2pJ4pJ6pJ2s_scaled, options, prms, error_tol, ms_iterMax, stepSize, PO_guess(2), PO_guess(4), PO_guess(6));
    elseif shooter_z0Fixed
        [PO_result, counter, constraint_error] = correctPO_multShooter_stateContinuity_z0Fixed(PO_guess, n_Nodes, @Int_CR3BnSTM_J2pJ4pJ6pJ2s_scaled, options, prms, error_tol, ms_iterMax, stepSize);
    elseif shooter_zd0Fixed
        [PO_result, counter, constraint_error] = correctPO_multShooter_stateContinuity_zd0Fixed(PO_guess, n_Nodes, @Int_CR3BnSTM_J2pJ4pJ6pJ2s_scaled, options, prms, error_tol, ms_iterMax, stepSize);
    elseif shooter_TpFixed
        [PO_result, counter, constraint_error] = correctPO_multShooter_stateContinuity_TpFixed(PO_guess, n_Nodes, @Int_CR3BnSTM_J2pJ4pJ6pJ2s_scaled, options, prms, error_tol, ms_iterMax, stepSize, PO_guess(7));
    elseif shooter_y0TpFixed
        [PO_result, counter, constraint_error] = correctPO_multShooter_stateContinuity_y0TpFixed(PO_guess, n_Nodes, @Int_CR3BnSTM_J2pJ4pJ6pJ2s_scaled, options, prms, error_tol, ms_iterMax, stepSize, PO_guess(2), PO_guess(7));
    elseif shooter_y0xd0TpFixed
        [PO_result, counter, constraint_error] = correctPO_mS_sC_y0xd0TpFixed(PO_guess, n_Nodes, @Int_CR3BnSTM_J2pJ4pJ6pJ2s_scaled, options, prms, error_tol, ms_iterMax, stepSize, PO_guess(2), PO_guess(4), PO_guess(7));
    elseif shooter_y0xd0zd0TpFixed
        [PO_result, counter, constraint_error] = correctPO_mS_sC_y0xd0zd0TpFixed(PO_guess, n_Nodes, @Int_CR3BnSTM_J2pJ4pJ6pJ2s_scaled, options, prms, error_tol, ms_iterMax, stepSize, PO_guess(2), PO_guess(4), PO_guess(6), PO_guess(7));
    end
    warning('on')
    %%% Calculate propagation error
    [T_result, X_result] = ode113(@Int_CR3BnSTM_J2pJ4pJ6pJ2s_scaled, [0, PO_result(7)], [PO_result(1:6); stm0_colVec], options, prms);
    figure(10); plot3(X_result(:,1), X_result(:,2), X_result(:,3), 'linewidth', 1, 'color', color_spectrum(step_i,:))
    propError = norm(X_result(end,1:6)' - X_result(1,1:6)') / norm(X_result(1,1:6));
    %%% Store solution
    PO_solutions(step_i,:) = PO_result';
    fprintf('Step = %1d ... iterations: %1.3d ... prop error = %1.2e ... Elapsed time: %1.1f seconds\n',...
    step_i, counter, propError, toc(ticWhole))
    %%% Inialize for next step
    PO_guess = PO_result;
        if shooter_TpFixed || shooter_y0TpFixed || shooter_y0xd0TpFixed || shooter_y0xd0zd0TpFixed
            Tp_fundamentalLyap_last = Tp_fundamentalLyap_current;
        end
    end





    fprintf('PO Result, row:\n')
    prettyColVec(PO_result, 'row')
    fprintf('PO Result, column:\n')
    prettyColVec(PO_result)
    figure(1)
    plot3(X_result(:,1),X_result(:,2),X_result(:,3), 'linewidth',2, 'color', color_spectrum(end,:))











end