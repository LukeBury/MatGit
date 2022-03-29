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
% savePath = '~/CU_Google_Drive/Documents/MatGit/mbin/Data/InitialConditions/PO_Families/';
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

%%% Plotting switches
plot_reference_PO = true;
plot_current_PO = true;
    plotSkip = true; % Only plot every X POs

plot_stability  = true;
plot_JC     = true;

%%% Save switches
save_PO_X0_database = 0;
if save_PO_X0_database == 1
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
POName = 'L2_Lyapunov'; myPO_ICs = PO_ICs.Jupiter_Europa.CR3BP.L2_Lyapunov; 
% POName = 'L2_Vertical'; myPO_ICs = PO_ICs.Jupiter_Europa.CR3BP.L2_Vertical; 
% POName = 'L2_NHalo'; myPO_ICs = PO_ICs.Jupiter_Europa.CR3BP.L2_NHalo; 
% POName = 'L2_SHalo'; myPO_ICs = PO_ICs.Jupiter_Europa.CR3BP.L2_SHalo; 
% POName = 'L2_EasternAxial'; myPO_ICs = PO_ICs.Jupiter_Europa.CR3BP.L2_EasternAxial; 
% POName = 'L2_WesternAxial'; myPO_ICs = PO_ICs.Jupiter_Europa.CR3BP.L2_WesternAxial; 

systemName = 'Jupiter_Europa.CR3BP';
famName = [systemName, '.', POName];

myPO_ICs = [1.000352554165237 - 0.0000;
 0.000000000000000;
 0.006051942705126 + 0.000000;
 -0.000000000000004 + 0.000;
 0.085815067810111;
 0.000000000000001 + 0.0000;
 2.301345119115549];

bump = [-0.000000000011878;
 -0.120937348250111;
 -0.000000000026733;
 -0.183172205543633;
 0.000000000256463;
 0.975603151127984];

myPO_ICs(1:6) = myPO_ICs(1:6) - bump.*0.005;
% myPO_ICs(7) = myPO_ICs(7)*1.1;

% ------------------------------------------------- 
%%% Continuation options
% -------------------------------------------------
%%% Maximum number of POs to find in family, regardless of step size
n_POs_max = 100;

%%% Continuation direction
familyDirection = -1; % 1 or -1

%%% Artificial scalar for increasing or decreasing the step size for 
% guessing at the state of the 2nd family member. 1 has no effect, 2 and 
% 3 can be helpful in certain regions
stepSize = 1;

%%% Maximum number of times the alorithm can try to relocate the family if
%%% it jumps to a different family/equilibrium 
lostFamilyCounter_max = 5;

%%% A value to scale down the continuation step size if a solution is not
%%% found within the maximum number of iterations. The step size will be
%%% updated and a solution will be tried for again.
stepSizePenalty = 0.9;

%%% Maximum percent change of position between successive PO_ICs. This is
%%% merely an imperfect way to try and detect if the algorithm took a step 
%%% to a different PO family.
max_percentChangeOfPosition = 0.1;

%%% If the continuation algorithm jumps to another PO family or
%%% equilibrium, this penalty will scale down the continaution step size so
%%% the algorithm can try again. 
lostFamilyPenalty = 0.9;

%%% If the multiple shooter required less iterations than this number, then
%%% scale up the continuation step size for the next PO
iterLimitForStepSizeIncrease = 100;

%%% Scale the continuation step size up by this amount if the iteration
%%% limit for step size increase was met
% iterLimitStepSizeReward = 1.2;
iterLimitStepSizeReward = 1.01;

% ------------------------------------------------- 
%%% Shooter options
% -------------------------------------------------
%%% Initial size to next PO in family. This value will dynamically change
%%% based on performance as the algorithm runs
% ds_PO = 2e-2;
ds_PO = 1e-2;
% ds_PO = 1e-3; % Good guess for Europa
% ds_PO = 1e-4; 
% ds_PO = 1e-5;
% ds_PO = 1e-6;

%%% Error tolerance for constraint vector norm in multiple shooter
% error_tol = 1e-12; 
error_tol = 1e-10; 

%%% Number of nodes for multiple shooter
N_Nodes = 2; % Generally, higher for bigger POs

%%% Maximum number of multiple-shooter iterations for any given family
%%% member before kicking out of the loop and adjusting the tuning
%%% parameters
ms_iterMax = 500;


% ------------------------------------------------- 
%%% Integration Options
% -------------------------------------------------
tol = 1e-13;
options = odeset('RelTol',tol,'AbsTol',tol);





% ========================================================================
%%% Continuation Prep
% ========================================================================
% --------------------------
% Set primary & secondary
% --------------------------
[primary, secondary] = assignPrimaryAndSecondary_CR3BP(famName, bodies);

% -------------------------------------------------
%%% System
% -------------------------------------------------
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



% ------------------------------------------------- 
%%% Preallocating and setting first index results
% -------------------------------------------------
%%% Preallocating space for outputs
Fs                      = NaN(n_POs_max, 6*N_Nodes + 1);
POs                     = NaN(n_POs_max,7);
stabilityIndices        = NaN(n_POs_max,2);
alphas_betas            = NaN(n_POs_max,2);
jacobiConstants         = NaN(n_POs_max,1);
landingVelocities_mps   = NaN(n_POs_max,1);

%%% Storing first PO
POs(1,:) = myPO_ICs(1:7)';

%%% Initialize STM as column vector to be added with state
stm0_colVec = reshape(eye(6),36,1);

%%% Initializing a counter for helping to identify when a family has been
%%% lost track of
lostFamilyCounter = 0;



% ------------------------------------------------- 
%%% Integrate reference trajectory and create nodes for multiple shooter
% -------------------------------------------------
%%% Integrating and plotting reference trajectory
[Tref_n, Xref_n] = ode113(@Int_CR3BnSTM, [0, myPO_ICs(7)], [myPO_ICs(1:6); stm0_colVec], options, prms);

if plot_reference_PO
    figure(1); hold all
    plot3(Xref_n(:,1),Xref_n(:,2),Xref_n(:,3),'m','linewidth',2)
    PlotBoi3_CR3Bn(26)
end

%%% Integrating reference trajectory and determining initial set of nodes
%%% which are evenly spaced in time
[~, X_nodes] = get_nodes(myPO_ICs(1:6), [0, myPO_ICs(7)], N_Nodes+1, @Int_CR3Bn, options, prms);
X_nodes = X_nodes(1:N_Nodes,1:6);

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
warning('This is off')
% [S1, S2]                            = getStabilityIndices(diag(eigenValues_new));
% stabilityIndices(1,:)               = [S1, S2];
[alpha, beta]                       = getBrouckeStabilityParameters(diag(eigenValues_new), monodromy);
alphas_betas(1,:)                   = [alpha, beta];
% ------------------------------------------------- 
%%% Create first free-variable column vector
% -------------------------------------------------
F_new   = reshape(X_nodes',6*N_Nodes,1);
F_new   = [F_new; myPO_ICs(7)];
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
%     if PO_i == 2
%         return
%     end
    %%% Use pseudo-arclength continuation to obtain new F, c, and DF
    [F_new, constraint_error, DF, counter] = pseudoArclengthContinuation_multShoot_stateContinuity(...
            error_tol, ms_iterMax, N_Nodes, F_new, @Int_CR3BnSTM, options, prms, PO_i, stepSize, nullVecDF, ds_PO);
    
    %%% Print info for current PO_i
    fprintf('PO_i = %1d ... iterations: %1.3d ... ds_PO = %1.2e ... Elapsed time: %1.1f seconds\n',...
        PO_i, counter, ds_PO, toc(ticWhole))
    
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
                warning('Null space changed dimensions')
                break
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
        
        % -------------------------------------------------
        %%% Check for family continuity (no jumps to other equilibria)
        % -------------------------------------------------
        percentChangeOfPosition = norm(F_new(1:3) - Fs(PO_i,1:3)')/norm(Fs(PO_i,1:3)');
        if (percentChangeOfPosition > max_percentChangeOfPosition) || (norm(F_new(1:3) - rLPs_n(LP,:)') < 1e-14)
            
            %%% If it seems we've jumped to a new family or equilibria,
            %%% lower the step size and try again. If this happens more
            %%% than the maximum allowable times, end the loop.
            if lostFamilyCounter < lostFamilyCounter_max
                %%% Count the lost family iteration
                lostFamilyCounter = lostFamilyCounter + 1;

                %%% Decrease step size and try to find family again
                ds_PO = ds_PO * lostFamilyPenalty;
                continue
            
            %%% If algorithm repeatedly fails to continue family, end
            %%% the search
            elseif lostFamilyCounter >= 5
                warning('Lost track of family ... It possibly ended')
                break
            end
        end
        
        % -------------------------------------------------
        %%% New solution is good - store it
        % -------------------------------------------------
        %%% Reset counter for when the algorithm loses the family
        lostFamilyCounter = 0;
        
        %%% Store converged F
        Fs(PO_i+1,:) = F_new';
        POs(PO_i+1,:) = [F_new(1:6)', F_new(end)];

        % -------------------------------------------------
        %%% Integrate the converged trajectory and store the stability
        %%% indices, the energy, and the approximate landing velocity
        % -------------------------------------------------
        %%% Integrate new PO
        [~, X_new] = ode113(@Int_CR3BnSTM, [0, F_new(end)], [F_new(1:6); stm0_colVec], options, prms);
        
        %%% Plotting PO
        if plot_current_PO
            if mod(PO_i,plotSkip) == 0
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
        [alpha, beta] = getBrouckeStabilityParameters(diag(eigenValues_new), monodromy);
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
        %%% Update PO_i
        PO_i = PO_i + 1;
        
    else % (if constraint error is too large)
        %%% If this was the first PO, warn that a 2nd couldn't be found
        if PO_i == 1
            warning('Failed to converge on a 2nd PO')
        end
        
        %%% Lower the continuation step size
        ds_PO = stepSizePenalty * ds_PO; 
    end % if constraint_error > error_tol
        
    
end % while PO_i <= n_POs_max





% ========================================================================
%%% Plots/Studies
% ========================================================================
POs = POs(~isnan(POs(:,1)),:);


% -------------------------------------------------
%%% Plot Stability Indices
% -------------------------------------------------
%%% Find stability indices
s1_logicalInd = ~isnan(stabilityIndices(:,1));
s2_logicalInd = ~isnan(stabilityIndices(:,2));
% s1_logicalInd = ~isnan(stabilityIndices(1:1502,1));
% s2_logicalInd = ~isnan(stabilityIndices(1:1502,2));
S1 = stabilityIndices(s1_logicalInd,1);
S2 = stabilityIndices(s2_logicalInd,2);
    
if plot_stability == 1
    JCs_s1 = jacobiConstants(s1_logicalInd);
    JCs_s2 = jacobiConstants(s2_logicalInd);
    
    figure('position',[209 322 948 302])
    subplot(1,2,1); hold all
    p1 = plot(JCs_s1, abs(S1),'o','markeredgecolor',colors.blue,'markerfacecolor',colors.ltblue);
    p2 = plot(JCs_s2, abs(S2),'o','markeredgecolor',colors.red,'markerfacecolor',colors.ltred);
    plot(unique([min([JCs_s2, JCs_s1]) max([JCs_s2, JCs_s1])]),[2 2],'k','linewidth',1)
    PlotBoi2('Jacobi Constant','Stability Indices',26,'LaTex')
    legend([p1 p2],'S_1','S_2','FontSize',14)
    xlim(unique([min([JCs_s2, JCs_s1]) max([JCs_s2, JCs_s1])]))
    
    subplot(1,2,2); hold all
    p1 = plot(JCs_s1, abs(S1),'o','markeredgecolor',colors.blue,'markerfacecolor',colors.ltblue);
    p2 = plot(JCs_s2, abs(S2),'o','markeredgecolor',colors.red,'markerfacecolor',colors.ltred);
    plot(unique([min([JCs_s2, JCs_s1]) max([JCs_s2, JCs_s1])]),[2 2],'k','linewidth',1)
    PlotBoi2('Jacobi Constant','',26,'LaTex')
    ylim([1.9 2.1])
    xlim(unique([min([JCs_s2, JCs_s1]) max([JCs_s2, JCs_s1])]))
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
    
end






alphas_logicalInd = ~isnan(alphas_betas(:,1));
betas_logicalInd = ~isnan(alphas_betas(:,2));
% s1_logicalInd = ~isnan(stabilityIndices(1:1502,1));
% s2_logicalInd = ~isnan(stabilityIndices(1:1502,2));
alphas = alphas_betas(alphas_logicalInd,1);
betas = alphas_betas(betas_logicalInd,2);
plot_BrouckeStabilityDiagram(alphas, betas)





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
















