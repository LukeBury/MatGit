% ========================================================================
%%% Description
% ========================================================================
% Script for continuing families of periodic orbits in Jupiter-Europa
% system with pseudo-arclength continuation

% Created: 12/02/20
% Author : Luke Bury, luke.bury@colorado.edu
% ========================================================================
%%% Initialization
% ========================================================================
clear
clc
close all
mbinPath = '~/CU_Google_Drive/Documents/MatGit/mbin';
savePath = '~/CU_Google_Drive/Documents/MatGit/mbin/Data/InitialConditions/PO_Families/';
addpath(genpath(mbinPath))
ticWhole = tic;


% is it that the case that the null space increases dimension near a tangent bifurcation but
% reduces dimension near end of family? If so, I should change it so that my ds_PO increases 
% if it gains a dimension so that I can step past the bifurcation
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
%%% Which shooter to run? true or false
run_y0Fixed                         = true;

run_y0xd0Fixed                      = false;
    artificially_set_xd0_equal_0    = false; % ---
    
run_y0zd0Fixed                      = false;
    artificially_set_zd0_equal_0    = false; % ---
    
run_y0xd0zd0Fixed                   = false;
    artificially_set_xd0zd0_equal_0 = false; % ---
    
run_y0TpFixed                       = false; % NOT PSEUDO-ARCLENGTH
run_y0xd0TpFixed                    = false; % NOT PSEUDO-ARCLENGTH
run_y0xd0zd0TpFixed                 = false; % NOT PSEUDO-ARCLENGTH
%     artificially_set_xd0zd0_equal_0 = true;

%%% Correct initial guess before continuation?
correct_initial_guess_to_y0equals0 = true;

%%% Plotting switches
plot_reference_PO = true;
plot_current_PO = true;
    plotSkip = true; % Only plot every X POs

plot_stability  = true;
plot_JC         = true;

%%% Save switches
save_PO_database         = 0;
save_zMirror_PO_database = 0;
if (save_PO_database == 1) || (save_zMirror_PO_database == 1)
    warning('Save is on')
end


% ========================================================================
%%% Personal Setup
% ========================================================================
% ------------------------------------------------- 
%%% Choose PO family to find
% -------------------------------------------------
% --------------------------
% Jupiter-Europa
% --------------------------
% POName = 'L2_Lyapunov'; myPO_ICs = PO_ICs.Jupiter_Europa.CR3BP.L2_Lyapunov; 
% POName = 'L2_Vertical'; myPO_ICs = PO_ICs.Jupiter_Europa.CR3BP.L2_Vertical; 
% POName = 'L2_NHalo'; myPO_ICs = PO_ICs.Jupiter_Europa.CR3BP.L2_NHalo; 
% POName = 'L2_SHalo'; myPO_ICs = PO_ICs.Jupiter_Europa.CR3BP.L2_SHalo; 
% POName = 'L2_EasternAxial'; myPO_ICs = PO_ICs.Jupiter_Europa.CR3BP.L2_EasternAxial; 
% POName = 'L2_WesternAxial'; myPO_ICs = PO_ICs.Jupiter_Europa.CR3BP.L2_WesternAxial; 

% POName = 'L2_L_P2'; myPO_ICs = PO_ICs.Jupiter_Europa.CR3BP.L2_L_1P2;
% POName = 'L2_L_2T'; myPO_ICs = PO_ICs.Jupiter_Europa.CR3BP.L2_L_2T;
% POName = 'L2_L_1T_P2'; myPO_ICs = PO_ICs.Jupiter_Europa.CR3BP.L2_L_1T_1P2;
% POName = 'L2_L_1T_P3'; myPO_ICs = PO_ICs.Jupiter_Europa.CR3BP.L2_L_1T_1P3;
POName = 'L2_Temp';
systemName = 'Jupiter_Europa';
% systemName = 'Jupiter_Ganymede';
% systemName = 'Saturn_Enceladus';
% systemName = 'Saturn_Titan';
% systemName = 'Neptune_Triton';
% systemName = 'Earth_Moon';

modelName = 'CR3BP';
famName = [systemName, '.', modelName, '.', POName];


% ----------------------------------------

% Test
% % %%% Crazy one that's big
% % myPO_ICs = [1.493091065609089;
% %  0.000000000000000;
% %  0.000000000000000;
% %  -0.396739350276297;
% %  -0.750049037320416;
% %  -0.000000000000000;
% %  12.562542503358495];


myPO_ICs = [1.0058908442858501;
 0.0000000000000000;
 0.0058960075953716;
 0.0000000000000000;
 0.0620955148570595;
 0.0000000000000000;
 5.9508643351844785];

% ------------------------------------------------- 
%%% Continuation options
% -------------------------------------------------
%%% Maximum number of POs to find in family, regardless of step size
n_POs_max = 5000;

%%% Continuation direction for use with pseudo arclength continuation
familyDirection = -1; % 1 or -1

%%% Artificial scalar for increasing or decreasing the step size for 
% guessing at the state of the 2nd family member. 1 has no effect, 2 and 
% 3 can be helpful in certain regions
stepSize = 1;

%%% Maximum number of times the alorithm can try to relocate the family if
%%% it jumps to a different family/equilibrium 
lostFamilyCounter_max = 10;

%%% A value to scale down the continuation step size if a solution is not
%%% found within the maximum number of iterations. The step size will be
%%% updated and a solution will be tried for again.
stepSizePenalty = 0.9;

%%% Maximum percent change of position between successive PO_ICs. This is
%%% merely an imperfect way to try and detect if the algorithm took a step 
%%% to a different PO family.
max_percentChangeOfSecondaryCenteredXPosition = 5; % percent!
% max_percentChangeOfSecondaryCenteredXPosition = 10; % percent!
% max_percentChangeOfSecondaryCenteredXPosition = 25; % percent!

%%% If the continuation algorithm jumps to another PO family or
%%% equilibrium, this penalty will scale down the continaution step size so
%%% the algorithm can try again. 
lostFamilyPenalty = 0.9;

%%% If the multiple shooter required less iterations than this number, then
%%% scale up the continuation step size for the next PO
iterLimitForStepSizeIncrease = 100;

%%% Scale the continuation step size up by this amount if the iteration
%%% limit for step size increase was met
% iterLimitStepSizeReward = 1;
% iterLimitStepSizeReward = 1.01; % slow
iterLimitStepSizeReward = 1.02; % medium
% iterLimitStepSizeReward = 1.05; % fast
% iterLimitStepSizeReward = 1.1; % real fast

%%% If the null space changes dimensions, we'll reduce the step size and
%%% try again (which seems to work pretty often). This parameter is the
%%% maximum number of times we can start over before exiting the loop
nullSpaceDimensionCounter_max = 30;

%%% A value to scale down the continuation step size if the null space
%%% dimensions changes
% nullSpaceDimensionChangePenalty = 0.95;
nullSpaceDimensionChangePenalty = 0.98;

%%% Certain penalty flags will gradually lower the continuation step size,
%%% espcially in sensitive regions or toward the end of a family. This
%%% value sets a lower limit, so when the step size dips below this, the
%%% loop is terminated
ds_PO_minimumValue = 1e-8;

%%% Initial size to next PO in family. This value will dynamically change
%%% based on performance as the algorithm runs
% ds_PO = 2e-1;
% % ds_PO = 1e-1;
% ds_PO = 5e-2;
% ds_PO = 2e-2;
% ds_PO = 1e-2;
% ds_PO = 5e-3; % 
ds_PO = 1e-3; % Good guess for Europa
% ds_PO = 5e-4;
% ds_PO = 2e-4;
% ds_PO = 1e-4; 
% ds_PO = 5e-5;
% ds_PO = 1e-5;
% ds_PO = 4e-6;
% ds_PO = 1e-6;
% ds_PO = 1e-7;
% ds_PO = 1e-8;

%%% If continuing over time period, choose time period range
if run_y0TpFixed || run_y0xd0zd0TpFixed || run_y0xd0TpFixed
%     Tp_range = linspace(myPO_ICs(7),   14.674007955286328, 100); 
    
    Tp_range = myPO_ICs(7):1e-3:50; 
%     Tp_range = myPO_ICs(7):-1e-6:0; 


    n_POs_max = length(Tp_range);
    if n_POs_max > 5000
        n_POs_max = 5000;
        Tp_range = Tp_range(1:n_POs_max);
    end
end

warning('on')
% ------------------------------------------------- 
%%% Shooter options
% -------------------------------------------------

%%% Error tolerance for constraint vector norm in multiple shooter
% error_tol = 1e-8; 
% error_tol = 1e-9; 
% error_tol = 1e-10; 
% error_tol = 5e-11; 
% error_tol = 1e-11; 
% error_tol = 5e-12; 
error_tol = 1e-12; 
% error_tol = 5e-13; 
% error_tol = 1e-13; %
% error_tol = 5e-14; 
% error_tol = 1e-14; 

%%% Number of nodes for multiple shooter. Generally, higher for bigger POs,
%%% but also be aware that more nodes increase inherent error, so you may
%%% need to lower error tolerances for the shooter
% n_Nodes = 1; 
n_Nodes = 2;
n_Nodes = 3;
n_Nodes = 4;
n_Nodes = 5; %
% n_Nodes = 6; % %
% n_Nodes = 7;
% n_Nodes = 8; 
% n_Nodes = 9;
% n_Nodes = 10;
% n_Nodes = 11;
% n_Nodes = 12;
% n_Nodes = 13;
% n_Nodes = 14;
% n_Nodes = 15;
% n_Nodes = 16;


%%% Maximum number of multiple-shooter iterations for any given family
%%% member before kicking out of the loop and adjusting the tuning
%%% parameters
% ms_iterMax = 500;
ms_iterMax = 750; % standard
% ms_iterMax = 1000;
% ms_iterMax = 1500;

% ------------------------------------------------- 
%%% Integration Options
% -------------------------------------------------
tol = 1e-13;
options = odeset('RelTol',tol,'AbsTol',tol);

% -------------------------------------------------
%%% Parameter setup
% -------------------------------------------------
% --------------------------
% Set primary & secondary
% --------------------------
[primary, secondary] = assignPrimaryAndSecondary_CR3BP(famName, bodies);

% --------------------------
%%% System
% --------------------------
%%% Normalizing constants
% rNorm:  n <-> km
% tNorm:  n <-> sec
% vNorm:  n <-> km/sec
[rNorm, tNorm, vNorm] = cr3bp_norms(primary, secondary, bodies.constants.G);

%%% Setting parameters
prms.u  = secondary.MR;
prms.n  = 1;
prms.R2 = secondary.R_n;

%%% Collinear equillibrium points
rLPs_n = EquilibriumPoints(prms.u, prms.n);

%%% Identifying Lagrange point 
if contains(POName,'L1')
    LP = 1;
elseif contains(POName,'L2')
    LP = 2;
end

% ========================================================================
%%% Correct initial guess, preallocate data structures, and initialize
% ========================================================================
%%% Integrate and stop at y=0 with event function
options_yEquals0_terminal = odeset('Event',@event_yEqualsZero,'RelTol',tol,'AbsTol',tol);
options_yEquals0_positiveTerminal = odeset('Event',@event_yEqualsZero_positiveTerminal,'RelTol',tol,'AbsTol',tol);
options_yEquals0_nonTerminal = odeset('Event',@event_yEqualsZero_nonTerminal,'RelTol',tol,'AbsTol',tol);

if correct_initial_guess_to_y0equals0
    % ------------------------------------------------- 
    %%% Integrate forward and backward to nearest y=0 crossings so the shortest
    %%% of the two can be used - this ensures a guess closest to the initial
    % ------------------------------------------------- 
    if myPO_ICs(2) ~= 0
        %%% Integrate forward and backward
        [Tfix_fwd, Xfix_fwd, ~, ~, ~] = ode113(@Int_CR3BnSTM, [0, myPO_ICs(7)], [myPO_ICs(1:6); reshape(eye(6),36,1)], options_yEquals0_terminal, prms);
        [Tfix_bkwd, Xfix_bkwd, ~, ~, ~] = ode113(@Int_CR3BnSTM, [myPO_ICs(7), 0], [myPO_ICs(1:6); reshape(eye(6),36,1)], options_yEquals0_terminal, prms);
        
        %%% Use this to choose a y=0 crossing. Ideally looking for one
        %%% where xd and zd are also 0
        if 1+1==1
            [Tfix_fwd_non, Xfix_fwd_non, tEv_fwd_non, xEv_fwd_non, ~] = ode113(@Int_CR3BnSTM, [0, myPO_ICs(7)], [myPO_ICs(1:6); reshape(eye(6),36,1)], options_yEquals0_nonTerminal, prms);
            xEv_fwd_non(:,1:6)
        end
        %%% Times to y=0
        dt_fwd = abs(Tfix_fwd(end) - Tfix_fwd(1));
        dt_bkwd = abs(Tfix_bkwd(end) - Tfix_bkwd(1));

        %%% Find minmum of two times and use that state as new guess
        if dt_fwd < dt_bkwd
            PO_0_guess = [Xfix_fwd(end,1:6)'; myPO_ICs(7)];
        elseif dt_fwd >= dt_fwd
            PO_0_guess = [Xfix_bkwd(end,1:6)'; myPO_ICs(7)];
        end
    else
        PO_0_guess = myPO_ICs;
    end

    %%% Plot initial guess
    [T_guess, X_guess] = ode113(@Int_CR3BnSTM, [0, myPO_ICs(7)], [myPO_ICs(1:6); reshape(eye(6),36,1)], options, prms);
    figure; hold all
    title('guess')
    plot3(X_guess(:,1),X_guess(:,2),X_guess(:,3),'m')
    PlotBoi3_CR3Bn(26)

    %%% Correct initial guess into a satisfactory PO
    if run_y0Fixed
        [PO_0, counter, constraint_error] = correctPO_multShooter_stateContinuity_y0Fixed(PO_0_guess, n_Nodes, @Int_CR3BnSTM, options, prms, error_tol, ms_iterMax, stepSize, 0);
    elseif run_y0xd0Fixed
        [PO_0, counter, constraint_error] = correctPO_multShooter_stateContinuity_y0xd0Fixed(PO_0_guess, n_Nodes, @Int_CR3BnSTM, options, prms, error_tol, ms_iterMax, stepSize, 0, 0);
    elseif run_y0zd0Fixed
        [PO_0, counter, constraint_error] = correctPO_multShooter_stateContinuity_y0zd0Fixed(PO_0_guess, n_Nodes, @Int_CR3BnSTM, options, prms, error_tol, ms_iterMax, stepSize, 0, 0);
    elseif run_y0xd0zd0Fixed
        [PO_0, counter, constraint_error] = correctPO_multShooter_stateContinuity_y0xd0zd0Fixed(PO_0_guess, n_Nodes, @Int_CR3BnSTM, options, prms, error_tol, ms_iterMax, stepSize, 0, 0, 0);
    elseif run_y0TpFixed
        [PO_0, counter, constraint_error] = correctPO_multShooter_stateContinuity_y0TpFixed(PO_0_guess, n_Nodes, @Int_CR3BnSTM, options, prms, error_tol, ms_iterMax, stepSize, 0, PO_0_guess(7));
    elseif run_y0xd0TpFixed
        [PO_0, counter, constraint_error] = correctPO_mS_sC_y0xd0TpFixed(PO_0_guess, n_Nodes, @Int_CR3BnSTM, options, prms, error_tol, ms_iterMax, stepSize, 0, 0, PO_0_guess(7));
    elseif run_y0xd0zd0TpFixed
        [PO_0, counter, constraint_error] = correctPO_mS_sC_y0xd0zd0TpFixed(PO_0_guess, n_Nodes, @Int_CR3BnSTM, options, prms, error_tol, ms_iterMax, stepSize, 0, 0, 0, PO_0_guess(7));
    end
    
    fprintf('Initial guess has converged\n')

    % [T, X] = ode113(@Int_CR3BnSTM, [0, PO_0(7)], [PO_0(1:6); reshape(eye(6),36,1)], options, prms);
    % [T_corrected, X_corrected] = ode113(@Int_CR3BnSTM, linspace(0, PO_0(7), 10000), [PO_0(1:6); reshape(eye(6),36,1)], options, prms);
    [T_corrected, X_corrected] = ode113(@Int_CR3BnSTM, linspace(0, PO_0(7), 10000), [PO_0(1:6); reshape(eye(6),36,1)], options, prms);
    figure; hold all
    title('corrected')
    plot3(X_corrected(:,1),X_corrected(:,2),X_corrected(:,3),'b')
    PlotBoi3_CR3Bn(26)
else
    fprintf('Not correcting initial guess\n')
    PO_0 = myPO_ICs;
    end


% ------------------------------------------------- 
%%% Artificially change initial values if desired
% -------------------------------------------------
if artificially_set_xd0_equal_0
    PO_0(4) = 0;
elseif artificially_set_zd0_equal_0
    PO_0(6) = 0;
elseif artificially_set_xd0zd0_equal_0
    PO_0(4) = 0;
    PO_0(6) = 0;
end
% ------------------------------------------------- 
%%% Preallocating and Initializing
% -------------------------------------------------
%%% Preallocating space for outputs
if run_y0Fixed
    Fs                  = NaN(n_POs_max, 6*n_Nodes);
elseif run_y0xd0Fixed
    Fs                  = NaN(n_POs_max, 6*n_Nodes-1);
elseif run_y0zd0Fixed
    Fs                  = NaN(n_POs_max, 6*n_Nodes-1);
elseif run_y0xd0zd0Fixed
    Fs                  = NaN(n_POs_max, 6*n_Nodes-2);
end
POs                     = NaN(n_POs_max,7);
actualErrors            = NaN(n_POs_max,1);
impactFlags             = NaN(n_POs_max,1);
stabilityIndices        = NaN(n_POs_max,2);
alphas_betas            = NaN(n_POs_max,2);
jacobiConstants         = NaN(n_POs_max,1);
landingVelocities_mps   = NaN(n_POs_max,1);

%%% Storing first PO
POs(1,:) = PO_0(1:7)';

%%% Initialize STM as column vector to be added with state
stm0_colVec = reshape(eye(6),36,1);

%%% Initializing a counter for helping to identify when a family has been
%%% lost track of
lostFamilyCounter = 0;

%%% Initializing a counter for tracking when the null space changes
%%% dimension
nullSpaceDimensionChangeCounter = 0;

% ------------------------------------------------- 
%%% Integrate reference trajectory and create nodes for multiple shooter
% -------------------------------------------------
%%% Integrating and plotting reference trajectory
[Tref_n, Xref_n] = ode113(@Int_CR3BnSTM, [0, PO_0(7)], [PO_0(1:6); stm0_colVec], options, prms);

if plot_reference_PO
%     figure(100); hold all
    hFig = figure(100); set( hFig, 'Position', [2785 154 560 420]); hold all

    plot3(Xref_n(:,1),Xref_n(:,2),Xref_n(:,3),'m','linewidth',2)
    PlotBoi3_CR3Bn(26)
%     view(0, 90) % x-y
%     view(90, 0) % y-z
    view(0, 0) % x-z
end

%%% Storing first actual error
actualErrors(1) = norm(Xref_n(end,1:6)' - Xref_n(1,1:6)') / norm(Xref_n(1,1:6));

%%% Storing the first impact flag (0 if it never impacts, 1 if it does)
impactFlags(1) = min(rowNorm(Xref_n(:,1:3) - [1-prms.u,0,0])) < prms.R2;

%%% Integrating reference trajectory and determining initial set of nodes
%%% which are evenly spaced in time
[~, X_nodes] = get_nodes([PO_0(1:6); stm0_colVec], [0, PO_0(7)], n_Nodes+1, @Int_CR3BnSTM, options, prms);
X_nodes = X_nodes(1:n_Nodes,1:6);

% ------------------------------------------------- 
%%% Store data pertaining to the reference PO
% -------------------------------------------------
%%% Jacobi constant
jacobiConstants(1) = getJacobiConstant_ZH(POs(1,1:6), prms);

%%% Approximate landing velocity
landingVelocities_mps(1) = JC_2_approxLandingVelocity(jacobiConstants(1), prms, vNorm);

%%% Stability indicies
stm_tf_t0                           = reshape(Xref_n(end,7:42),6,6);
monodromy                           = stm_tf_t0;
[eigenVectors_new, eigenValues_new] = eig(monodromy);
[S1, S2]                            = getStabilityIndices(diag(eigenValues_new));
stabilityIndices(1,:)               = [S1, S2];
[alpha, beta]                       = getBrouckeStabilityParameters(diag(eigenValues_new), monodromy);
alphas_betas(1,:)                   = [alpha, beta];

% ------------------------------------------------- 
%%% Create first free-variable column vector
% -------------------------------------------------
%%% Start by including everything in the free-variable vector
F_new   = reshape(X_nodes',6*n_Nodes,1);
F_new   = [F_new; PO_0(7)];

%%% Take out necessary components from the full-free variable vector 
if run_y0Fixed
    F_new = [F_new(1); F_new(3:end)];
elseif run_y0xd0Fixed
    F_new = [F_new(1); F_new(3); F_new(5:end)];
elseif run_y0zd0Fixed
    F_new = [F_new(1); F_new(3:5); F_new(7:end)];
elseif run_y0xd0zd0Fixed
    F_new = [F_new(1); F_new(3); F_new(5); F_new(7:end)];
elseif run_y0TpFixed
    F_new = [F_new(1); F_new(3:end-1)];
elseif run_y0xd0TpFixed
    F_new = [F_new(1); F_new(3); F_new(5:end-1)];
elseif run_y0xd0zd0TpFixed
    F_new = [F_new(1); F_new(3); F_new(5); F_new(7:end-1)];
end
        
%%% Store first free variable vector
Fs(1,:) = F_new';
% ========================================================================
%%% Enter the continuation loop
% ========================================================================
%%% Initialize PO counter
PO_i = 1;

%%% Enter PO loop
while PO_i <= n_POs_max
    %%% Set free variable vector for current PO_i
    F_new = Fs(PO_i, :)'; 
    
    % -------------------------------------------------
    % Run pseudo-arclength continuation with multiple shooting to find the
    % next orbit
    % -------------------------------------------------
    %%% If this is the first PO, set null(DF) as zeros
    if PO_i == 1
        nullVecDF = zeros(size(F_new));
    end

    %%% Use pseudo-arclength continuation to obtain new F, c, and DF
    if run_y0Fixed
        [F_new, constraint_error, DF, counter] = pseudoArclengthContinuation_multShoot_stateContinuity_y0Fixed(...
            error_tol, ms_iterMax, n_Nodes, F_new, @Int_CR3BnSTM, options, prms, PO_i, stepSize, nullVecDF, ds_PO, 0);
    elseif run_y0xd0Fixed
        [F_new, constraint_error, DF, counter] = pAC_mS_sC_y0xd0Fixed(...
            error_tol, ms_iterMax, n_Nodes, F_new, @Int_CR3BnSTM, options, prms, PO_i, stepSize, nullVecDF, ds_PO, 0, 0);
    elseif run_y0zd0Fixed
        [F_new, constraint_error, DF, counter] = pAC_mS_sC_y0zd0Fixed(...
            error_tol, ms_iterMax, n_Nodes, F_new, @Int_CR3BnSTM, options, prms, PO_i, stepSize, nullVecDF, ds_PO, 0, 0);
    elseif run_y0xd0zd0Fixed
        [F_new, constraint_error, DF, counter] = pAC_mS_sC_y0xd0zd0Fixed(...
            error_tol, ms_iterMax, n_Nodes, F_new, @Int_CR3BnSTM, options, prms, PO_i, stepSize, nullVecDF, ds_PO, 0, 0, 0);
    elseif run_y0TpFixed
        Tp_i = Tp_range(PO_i);
        [PO_new, counter, constraint_error] = correctPO_multShooter_stateContinuity_y0TpFixed(...
            [F_new(1); 0; F_new(2:5)], n_Nodes, @Int_CR3BnSTM, options, prms, error_tol, ms_iterMax, stepSize, 0, Tp_i);
    elseif run_y0xd0TpFixed
        Tp_i = Tp_range(PO_i);
        PO_newGuess = POs(PO_i,:)';
        [PO_new, counter, constraint_error] = correctPO_mS_sC_y0xd0TpFixed(...
            PO_newGuess(1:6), n_Nodes, @Int_CR3BnSTM, options, prms, error_tol, ms_iterMax, stepSize, 0, 0, Tp_i);
        
    elseif run_y0xd0zd0TpFixed
        Tp_i = Tp_range(PO_i);
        PO_newGuess = POs(PO_i,:)';
        [PO_new, counter, constraint_error] = correctPO_mS_sC_y0xd0zd0TpFixed(...
            PO_newGuess(1:6), n_Nodes, @Int_CR3BnSTM, options, prms, error_tol, ms_iterMax, stepSize, 0, 0, 0, Tp_i);
    end
    
    % -------------------------------------------------
    % If error tolerance for the constraint vector is met, then compute the
    % null space of DF, adjust step-size if desired, and store the
    % converged free-variable vector. 
    %
    % If the tolerance was not met, lower the continuation step size
    % (ds_PO) and try again
    % -------------------------------------------------
    
    if constraint_error <= error_tol
        %%% Compute null space of DF
%         if (~run_y0TpFixed && ~run_y0xd0zd0TpFixed)
        if (~run_y0TpFixed && ~run_y0xd0zd0TpFixed && ~run_y0xd0TpFixed)
            if PO_i == 1
                nullVecDF = null(DF).*familyDirection;

            % -------------------------------------------------
            %%% If PO_i > 1, make sure the null space hasn't changed dimension
            %%% or flipped sign
            % -------------------------------------------------
            elseif PO_i > 1
                nullVecDFPrevious = nullVecDF;
                nullVecDF = null(DF).*familyDirection;

                %%% Check continuity of null space dimension and quit if it has
                %%% changed
                if ~isequal(size(nullVecDFPrevious),size(nullVecDF))

                    if nullSpaceDimensionChangeCounter < nullSpaceDimensionCounter_max
                        nullSpaceDimensionChangeCounter = nullSpaceDimensionChangeCounter + 1;
                        ds_PO = ds_PO * nullSpaceDimensionChangePenalty;

                        if ds_PO <= ds_PO_minimumValue
                            warning('Step size decreased below minimum set value')
                            break
                        end

                        if isempty(nullVecDF)
                            warning('Null space of DF collapsed ... Decreasing step size and trying again')
                        else
                            warning('Null space of DF gained a dimension ... Decreasing step size and trying again')
                        end
                        nullVecDF = nullVecDFPrevious;
                        continue
                    else
                        warning('Null space seems to have definitely changed dimensions')
                        break
                    end

                end

                %%% Check sign
                nullDot = nullVecDFPrevious' * nullVecDF;

                %%% Modify null vector so there's no sign change
                for kk = 1:size(nullVecDF,2)
                    if isequal(sign(nullDot(kk,kk)), -1)
                        nullVecDF(:,kk) = -1*nullVecDF(:,kk);
                    end
                end
            end
        end
        
        % -------------------------------------------------
        %%% Check for family continuity (no jumps to other equilibria)
        % -------------------------------------------------
%         percentChangeOfPosition = norm(F_new(1:2) - Fs(PO_i,1:2)')/norm(Fs(PO_i,1:2)');
%         if (~run_y0TpFixed && ~run_y0xd0zd0TpFixed)
        if (~run_y0TpFixed && ~run_y0xd0zd0TpFixed && ~run_y0xd0TpFixed)
            if PO_i >= 3
                oldXPosition = norm(Fs(PO_i, 1)' - (1-prms.u));
                newXPosition = norm(F_new(1) - (1-prms.u));
                percentChangeOfSecondaryCenteredXPosition = norm(newXPosition - oldXPosition)*100 / oldXPosition;
                if (percentChangeOfSecondaryCenteredXPosition > max_percentChangeOfSecondaryCenteredXPosition) || (norm([F_new(1); 0; F_new(2)] - rLPs_n(LP,:)') < 1e-14)

                    %%% If it seems we've jumped to a new family or equilibria,
                    %%% lower the step size and try again. If this happens more
                    %%% than the maximum allowable times, end the loop.
                    if lostFamilyCounter < lostFamilyCounter_max
                        %%% Count the lost family iteration
                        lostFamilyCounter = lostFamilyCounter + 1;

                        %%% Decrease step size and try to find family again
                        ds_PO = ds_PO * lostFamilyPenalty;

                        if ds_PO <= ds_PO_minimumValue
                            warning('Step size decreased below minimum set value')
                            break
                        end

                        warning('Lost track of family ... Trying again')
                        continue

                    %%% If algorithm repeatedly fails to continue family, end
                    %%% the search
                    else
                        warning('Lost track of family ... It possibly ended')
                        break
                    end
                end
            end
        end
        
        % -------------------------------------------------
        %%% New solution is good - store it
        % -------------------------------------------------
        %%% Reset counter for when the null space changes dimension
        nullSpaceDimensionChangeCounter = 0;
        
        %%% Reset counter for when the algorithm loses the family
        lostFamilyCounter = 0;
        
        %%% Store converged F
        if run_y0TpFixed
            [~, X_nodes] = get_nodes([PO_new(1:6); stm0_colVec], [0, PO_new(7)], n_Nodes+1, @Int_CR3BnSTM, options, prms);
            X_nodes = X_nodes(1:n_Nodes,1:6);
            F_new   = reshape(X_nodes',6*n_Nodes,1);
            F_new = [F_new(1); F_new(3:end)];
        end
        
        if run_y0xd0TpFixed
            [~, X_nodes] = get_nodes([PO_new(1:6); stm0_colVec], [0, PO_new(7)], n_Nodes+1, @Int_CR3BnSTM, options, prms);
            X_nodes = X_nodes(1:n_Nodes,1:6);
            F_new   = reshape(X_nodes',6*n_Nodes,1);
            F_new = [F_new(1); F_new(3); F_new(5:end)];
        end
        
        if run_y0xd0zd0TpFixed
            [~, X_nodes] = get_nodes([PO_new(1:6); stm0_colVec], [0, PO_new(7)], n_Nodes+1, @Int_CR3BnSTM, options, prms);
            X_nodes = X_nodes(1:n_Nodes,1:6);
            F_new   = reshape(X_nodes',6*n_Nodes,1);
            F_new = [F_new(1); F_new(3); F_new(5); F_new(7:end)];
        end
        
        Fs(PO_i+1,:) = F_new';
        
        %%% Store converged PO
        if run_y0Fixed
            POs(PO_i+1,:) = [F_new(1); 0; F_new(2:5); F_new(end)]';
        elseif run_y0xd0Fixed
            POs(PO_i+1,:) = [F_new(1); 0; F_new(2); 0; F_new(3:4); F_new(end)]';
        elseif run_y0zd0Fixed
            POs(PO_i+1,:) = [F_new(1); 0; F_new(2:4); 0; F_new(end)]';
        elseif run_y0xd0zd0Fixed
            POs(PO_i+1,:) = [F_new(1); 0; F_new(2); 0; F_new(3); 0; F_new(end)]';
        elseif run_y0TpFixed
            POs(PO_i+1,:) = PO_new';
        elseif run_y0xd0TpFixed
            POs(PO_i+1,:) = PO_new';
        elseif run_y0xd0zd0TpFixed
            POs(PO_i+1,:) = PO_new';
        end

        % -------------------------------------------------
        %%% Integrate the converged trajectory and store the stability
        %%% indices, the energy, and the approximate landing velocity
        % -------------------------------------------------
        %%% Integrate new PO
        if run_y0Fixed
            [~, X_new] = ode113(@Int_CR3BnSTM, [0, F_new(end)], [F_new(1); 0; F_new(2:5); stm0_colVec], options, prms);
        elseif run_y0xd0Fixed
            [~, X_new] = ode113(@Int_CR3BnSTM, [0, F_new(end)], [F_new(1); 0; F_new(2); 0; F_new(3:4); stm0_colVec], options, prms);
        elseif run_y0zd0Fixed
            [~, X_new] = ode113(@Int_CR3BnSTM, [0, F_new(end)], [F_new(1); 0; F_new(2:4); 0; stm0_colVec], options, prms);
        elseif run_y0xd0zd0Fixed
            [~, X_new] = ode113(@Int_CR3BnSTM, [0, F_new(end)], [F_new(1); 0; F_new(2); 0; F_new(3); 0; stm0_colVec], options, prms);
        elseif run_y0TpFixed
            [~, X_new] = ode113(@Int_CR3BnSTM, [0, Tp_i], [PO_new(1:6); stm0_colVec], options, prms);
        elseif run_y0xd0TpFixed
            [~, X_new] = ode113(@Int_CR3BnSTM, [0, Tp_i], [PO_new(1:6); stm0_colVec], options, prms);
        elseif run_y0xd0zd0TpFixed
            [~, X_new] = ode113(@Int_CR3BnSTM, [0, Tp_i], [PO_new(1:6); stm0_colVec], options, prms);
        end
        
        %%% Store the actual error of the fully integrated PO
        actualErrors(PO_i+1) = norm(X_new(end,1:6)' - X_new(1,1:6)') / norm(X_new(1,1:6));
        
        %%% Storing the impact flag (0 if it never impacts, 1 if it does)
        impactFlags(PO_i+1) = min(rowNorm(X_new(:,1:3) - [1-prms.u,0,0])) < prms.R2;
        
        %%% Plotting PO
        if plot_current_PO
            if mod(PO_i,plotSkip) == 0
                figure(100)
                plot3(X_new(:,1),X_new(:,2),X_new(:,3),'b')
                drawnow
            end
        end
        
        %%% Stability indices of new PO
        stm_tf_t0                           = reshape(X_new(end,7:42),6,6);
        monodromy                           = stm_tf_t0;
        [eigenVectors_new, eigenValues_new] = eig(monodromy);
        [S1, S2]                            = getStabilityIndices(diag(eigenValues_new));
        stabilityIndices(PO_i+1,:)          = [S1, S2];
        
        %%% Broucke stability parameters
        [alpha, beta]           = getBrouckeStabilityParameters(diag(eigenValues_new), monodromy);
        alphas_betas(PO_i+1, :) = [alpha, beta];
        
        %%% Energy of new PO
        jacobiConstants(PO_i+1)        = getJacobiConstant_ZH(X_new(1,1:6), prms);
        landingVelocities_mps(PO_i+1)  = JC_2_approxLandingVelocity(jacobiConstants(PO_i+1), prms, vNorm);
        
        % -------------------------------------------------
        %%% If the multiple shooter required few enough iterations, then
        %%% scale up the continuation step size for the next PO
        % -------------------------------------------------
        if counter <= iterLimitForStepSizeIncrease
            ds_PO = iterLimitStepSizeReward * ds_PO;
        end
        
        % -------------------------------------------------
        %%% Update PO_i for next while-loop iteration
        % -------------------------------------------------
        %%% Print info for current PO_i
%         if run_y0TpFixed || run_y0xd0zd0TpFixed
        if run_y0TpFixed || run_y0xd0zd0TpFixed || run_y0xd0TpFixed
            fprintf('PO_i = %1d ... iterations: %1.3d ... actual error = %1.2e ... Elapsed time: %1.1f seconds\n',...
                PO_i, counter, actualErrors(PO_i+1), toc(ticWhole))
        else
            fprintf('PO_i = %1d ... iterations: %1.3d ... ds_PO = %1.2e ... actual error = %1.2e ... Elapsed time: %1.1f seconds\n',...
                PO_i, counter, ds_PO, actualErrors(PO_i+1), toc(ticWhole))
        end
            
    
        %%% Update PO_i
        PO_i = PO_i + 1;
        
    else % (if constraint error is too large)
        %%% If this was the first PO, warn that a 2nd couldn't be found
        if PO_i == 1
            warning('Failed to converge on a 2nd PO')
        end
        
        %%% Lower the continuation step size
        ds_PO = stepSizePenalty * ds_PO; 
        fprintf('    Failed to converge in iteration limit ... lowering step size\n')
        
        if ds_PO <= ds_PO_minimumValue
            warning('Step size decreased below minimum set value')
            break
        end
    end % if constraint_error > error_tol
        
    
end % while PO_i <= n_POs_max


%%% Print the maximum actual error
fprintf('\nMaximum actual error: %1.2e\n\n', max(actualErrors))

% ========================================================================
%%% Plots/Studies
% ========================================================================
%%% Get rid of extra space
logIndices_POs = ~isnan(stabilityIndices(:,1));
% logIndices_POs = ~isnan(POs(1:174,1));
% logIndices_POs(1) = false;

POs              = POs(logIndices_POs, :);
actualErrors     = actualErrors(logIndices_POs);
impactFlags      = impactFlags(logIndices_POs);
stabilityIndices = stabilityIndices(logIndices_POs, :);
jacobiConstants  = jacobiConstants(logIndices_POs);
alphas_betas     = alphas_betas(logIndices_POs, :);
% -------------------------------------------------
%%% Plot Stability Indices
% -------------------------------------------------
if plot_stability == 1
    z_indices = linspace(1, size(POs,1), size(POs,1));
    
    figure('position',[209 322 948 302])
    subplot(1,2,1); hold all
    p1 = plot3(jacobiConstants, abs(stabilityIndices(:,1)), z_indices,'o','markeredgecolor',colors.blue,'markerfacecolor',colors.ltblue);
    p2 = plot3(jacobiConstants, abs(stabilityIndices(:,2)), z_indices,'o','markeredgecolor',colors.red,'markerfacecolor',colors.ltred);
    plot3(unique([min(jacobiConstants) max(jacobiConstants)]),[2 2],[1 1],'k','linewidth',1)
    PlotBoi2('Jacobi Constant','Stability Indices',26,'LaTex')
    legend([p1 p2],'S_1','S_2','FontSize',14)
    xlim(unique([min(jacobiConstants) max(jacobiConstants)]))
    view(0,90)
    
    subplot(1,2,2); hold all
    p1 = plot3(jacobiConstants, abs(stabilityIndices(:,1)), z_indices,'o','markeredgecolor',colors.blue,'markerfacecolor',colors.ltblue);
    p2 = plot3(jacobiConstants, abs(stabilityIndices(:,2)), z_indices,'o','markeredgecolor',colors.red,'markerfacecolor',colors.ltred);
    plot3(unique([min(jacobiConstants) max(jacobiConstants)]),[2 2],[1 1],'k','linewidth',1)
    PlotBoi2('Jacobi Constant','',26,'LaTex')
    ylim([1.9 2.1])
    xlim(unique([min(jacobiConstants) max(jacobiConstants)]))
    view(0,90)

end

% -------------------------------------------------
%%% Plot Jacobi constants
% -------------------------------------------------
if plot_JC
    figure; hold all
    yyaxis left
    plot(jacobiConstants,'o','markeredgecolor',colors.blue,'markerfacecolor',colors.ltblue)
    PlotBoi2('','Jacobi Constant',26,'LaTex')
    yyaxis right
    plot(landingVelocities_mps,'o','markeredgecolor',colors.red,'markerfacecolor',colors.ltred)
    PlotBoi2('PO Index','+$x$ Landing Velocity, $m/s$',26,'LaTex')
    
    figure; hold all
    plot3(POs(:,7), jacobiConstants, linspace(1,length(jacobiConstants), length(jacobiConstants)) ,'o','markeredgecolor',colors.blue,'markerfacecolor',colors.ltblue)
    PlotBoi2('$T_p$', 'Jacobi Constant', 26, 'LaTex')
end


% -------------------------------------------------
%%% Print bifurcation information from Broucke stability diagram
% -------------------------------------------------
[bifurcation_strings] = plot_BrouckeStabilityDiagram(alphas_betas(:,1), alphas_betas(:,2),false);

















% ========================================================================
%%% Saving Data
% ========================================================================
if save_PO_database == 1
    
    
    %%% Create unique filename and open file
    [uniqueFilename] = get_uniqueFilename(famName, savePath, 'txt');
    datafile = fopen(uniqueFilename,'wt');
    
    %%% Write header
    headerString = 'x0,y0,z0,xd0,yd0,zd0,Tp,JC,stabilityIndex1,stabilityIndex2,alpha,beta,impactFlag,error\n';
    fprintf(datafile,headerString);
    
    %%% Write data
    dataStart = 2;
    dataStop = size(POs,1);
    
    
    for kk = dataStart:dataStop % Order of generation
%     for kk = dataStop:-1:dataStart % Reverse-order of generation
        
        if POs(kk,2) == 0
            fprintf(datafile,'%1.16f,%1.1f,%1.16f,%1.16f,%1.16f,%1.16f,%1.16f,%1.16f,%1.6f,%1.6f,%1.13f,%1.13f,%1d,%1.1e\n',...
                POs(kk,1), POs(kk,2), POs(kk,3), POs(kk,4), POs(kk,5), POs(kk,6), POs(kk,end),...
                jacobiConstants(kk), stabilityIndices(kk,1), stabilityIndices(kk,2), alphas_betas(kk,1), alphas_betas(kk,2),...
                impactFlags(kk), actualErrors(kk));
        else
            fprintf(datafile,'%1.16f,%1.16f,%1.16f,%1.16f,%1.16f,%1.16f,%1.16f,%1.16f,%1.6f,%1.6f,%1.13f,%1.13f,%1d,%1.1e\n',...
                POs(kk,1), POs(kk,2), POs(kk,3), POs(kk,4), POs(kk,5), POs(kk,6), POs(kk,end),...
                jacobiConstants(kk), stabilityIndices(kk,1), stabilityIndices(kk,2), alphas_betas(kk,1), alphas_betas(kk,2),...
                impactFlags(kk), actualErrors(kk));
        end
    
        
    end
    
    %%% Close file
    fclose(datafile);




end


if save_zMirror_PO_database == 1
%     % -------------------------------------------------
%     %%% Preparing Save File
%     % -------------------------------------------------
%     %%% Create unique filename
%     [uniqueFilename] = get_uniqueFilename(famName, savePath, 'txt');
%     
%     %%% Open File
%     datafile = fopen(uniqueFilename,'wt');
%     
%     % -------------------------------------------------
%     %%% Writing data
%     % -------------------------------------------------
%     %%% Write header
%     headerString = ['x0,y0,z0,xd0,yd0,zd0,Tp,JC,stabilityIndex1,stabilityIndex2,alpha,beta,impactFlag,error\n'];
%     fprintf(datafile,headerString);
%     
%     %%% Write data
%     dataStart = 1;
%     dataStop = size(POs,1);
%     
%     
% %     for kk = dataStart:dataStop
%     for kk = dataStop:-1:dataStart
%         
%         if POs(kk,2) == 0
%             fprintf(datafile,'%1.16f,%1.1f,%1.16f,%1.16f,%1.16f,%1.16f,%1.16f,%1.16f,%1.6f,%1.6f,%1.13f,%1.13f,%1d,%1.1e\n',...
%                 POs(kk,1), POs(kk,2), -POs(kk,3), POs(kk,4), POs(kk,5), -POs(kk,6), POs(kk,end),...
%                 jacobiConstants(kk), stabilityIndices(kk,1), stabilityIndices(kk,2), alphas_betas(kk,1), alphas_betas(kk,2),...
%                 impactFlags(kk), actualErrors(kk));
%         else
%             fprintf(datafile,'%1.16f,%1.16f,%1.16f,%1.16f,%1.16f,%1.16f,%1.16f,%1.16f,%1.6f,%1.6f,%1.13f,%1.13f,%1d,%1.1e\n',...
%                 POs(kk,1), POs(kk,2), -POs(kk,3), POs(kk,4), POs(kk,5), -POs(kk,6), POs(kk,end),...
%                 jacobiConstants(kk), stabilityIndices(kk,1), stabilityIndices(kk,2), alphas_betas(kk,1), alphas_betas(kk,2),...
%                 impactFlags(kk), actualErrors(kk));
%         end
%     
%         
%     end
%     
%     %%% Close file
%     fclose(datafile);
%     
    
    
end



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
















