% ========================================================================
%%% Description
% ========================================================================
% Script for fleshing out a family of periodic orbits based on a set of
% initial conditions (X0 and Tp). The family is continued via
% pseudo-arclength continuation with a variable-node multiple shooter. 
% Generally, finding a new family requires tweaking some of the tuning
% parameters. Once satisfactory parameters are found and a full family is
% generated, turn on the option to save the X0 database to save all the
% ICs of the family, with Tp and JC included, to a csv in the savePath
% directory.

% Created: 09/09/16
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
print_PO_index  = 1;
plot_current_PO = 1;
    plotSkip = 1; % Only plot every X POs

plot_stability  = 1;
plot_energy     = 1;

save_PO_X0_database = 0;
if save_PO_X0_database == 1
    warning('Save is on')
end

% ========================================================================
%%% Setup
% ========================================================================
% ------------------------------------------------- 
%%% Choose PO family to find
% -------------------------------------------------
% --------------------------
% Earth-Moon
% --------------------------
% myPO_ICs = PO_ICs.Earth_Moon.CR3BP.L1_Lyapunov; famName = 'Earth_Moon.CR3BP.L1_Lyapunov';
% myPO_ICs = PO_ICs.Earth_Moon.CR3BP.L1_Vertical; famName = 'Earth_Moon.CR3BP.L1_Vertical';
% myPO_ICs = PO_ICs.Earth_Moon.CR3BP.L1_SHalo;    famName = 'Earth_Moon.CR3BP.L1_SHalo';

% myPO_ICs = PO_ICs.Earth_Moon.CR3BP.L2_Lyapunov; famName = 'Earth_Moon.CR3BP.L2_Lyapunov';


% --------------------------
% Jupiter-Europa
% --------------------------
myPO_ICs = PO_ICs.Jupiter_Europa.CR3BP.L2_Lyapunov; famName = 'Jupiter_Europa.CR3BP.L2_Lyapunov';
% myPO_ICs = PO_ICs.Jupiter_Europa.CR3BP.L2_Vertical; famName = 'Jupiter_Europa.CR3BP.L2_Vertical';
% myPO_ICs = PO_ICs.Jupiter_Europa.CR3BP.L2_NHalo; famName = 'Jupiter_Europa.CR3BP.L2_NHalo';
% myPO_ICs = PO_ICs.Jupiter_Europa.CR3BP.L2_SHalo; famName = 'Jupiter_Europa.CR3BP.L2_SHalo';
% myPO_ICs = PO_ICs.Jupiter_Europa.CR3BP.L2_EasternAxial; famName = 'Jupiter_Europa.CR3BP.L2_EasternAxial';
% myPO_ICs = PO_ICs.Jupiter_Europa.CR3BP.L2_WesternAxial; famName = 'Jupiter_Europa.CR3BP.L2_WesternAxial';

% myPO_ICs = PO_ICs.Jupiter_Europa.CR3BP_J2pJ4pJ6pJ2s.L2_Lyapunov; famName = 'Jupiter_Europa.CR3BP_J2pJ4pJ6pJ2s.L2_Lyapunov';
% myPO_ICs = PO_ICs.Jupiter_Europa.CR3BP_J2pJ4pJ6pJ2s.L2_SHalo; famName = 'Jupiter_Europa.CR3BP_J2pJ4pJ6pJ2s.L2_SHalo';
% myPO_ICs = PO_ICs.Jupiter_Europa.CR3BP_J2pJ4pJ6pJ2s.L2_Vertical; famName = 'Jupiter_Europa.CR3BP_J2pJ4pJ6pJ2s.L2_Vertical';

% --------------------------
% Jupiter-Ganymede
% --------------------------
% myPO_ICs = PO_ICs.Jupiter_Ganymede.CR3BP.L2_Lyapunov; famName = 'Jupiter_Ganymede.CR3BP.L2_Lyapunov';
% myPO_ICs = PO_ICs.Jupiter_Ganymede.CR3BP.L2_Vertical; famName = 'Jupiter_Ganymede.CR3BP.L2_Vertical';
% myPO_ICs = PO_ICs.Jupiter_Ganymede.CR3BP.L2_SHalo; famName = 'Jupiter_Ganymede.CR3BP.L2_SHalo';

% --------------------------
% Saturn-Enceladus
% --------------------------
% myPO_ICs = PO_ICs.Saturn_Enceladus.CR3BP.L2_Lyapunov; famName = 'Saturn_Enceladus.CR3BP.L2_Lyapunov';
% myPO_ICs = PO_ICs.Saturn_Enceladus.CR3BP.L1_Lyapunov; famName = 'Saturn_Enceladus.CR3BP.L1_Lyapunov';
% myPO_ICs = PO_ICs.Saturn_Enceladus.CR3BP.L2_SHalo; famName = 'Saturn_Enceladus.CR3BP.L2_SHalo';
% myPO_ICs = PO_ICs.Saturn_Enceladus.CR3BP.L2_NHalo; famName = 'Saturn_Enceladus.CR3BP.L2_NHalo';
% myPO_ICs = PO_ICs.Saturn_Enceladus.CR3BP.L1_SHalo; famName = 'Saturn_Enceladus.CR3BP.L1_SHalo';
% myPO_ICs = PO_ICs.Saturn_Enceladus.CR3BP.L1_NHalo; famName = 'Saturn_Enceladus.CR3BP.L1_NHalo';
% myPO_ICs = PO_ICs.Saturn_Enceladus.CR3BP.L2_Vertical; famName = 'Saturn_Enceladus.CR3BP.L2_Vertical';
% myPO_ICs = PO_ICs.Saturn_Enceladus.CR3BP.L2_WesternAxial; famName = 'Saturn_Enceladus.CR3BP.L2_WesternAxial';
% myPO_ICs = PO_ICs.Saturn_Enceladus.CR3BP.L2_EasternAxial; famName = 'Saturn_Enceladus.CR3BP.L2_EasternAxial';
% myPO_ICs = PO_ICs.Saturn_Enceladus.CR3BP.DRO; famName = 'Saturn_Enceladus.CR3BP.DRO';
% myPO_ICs = PO_ICs.Saturn_Enceladus.CR3BP.L1_NHalo_DoublePeriod; famName = 'Saturn_Enceladus.CR3BP.L1_NHalo_DoublePeriod';
% myPO_ICs = PO_ICs.Saturn_Enceladus.CR3BP.unlabeled1LeadingS; famName = 'Saturn_Enceladus.CR3BP.unlabeled1LeadingS';
% myPO_ICs = PO_ICs.Saturn_Enceladus.CR3BP.unlabeled1TrailingS; famName = 'Saturn_Enceladus.CR3BP.unlabeled1TrailingS';
% myPO_ICs = PO_ICs.Saturn_Enceladus.CR3BP.unlabeled1LeadingN; famName = 'Saturn_Enceladus.CR3BP.unlabeled1LeadingN';
% myPO_ICs = PO_ICs.Saturn_Enceladus.CR3BP.unlabeled1TrailingN; famName = 'Saturn_Enceladus.CR3BP.unlabeled1TrailingN';
% myPO_ICs = PO_ICs.Saturn_Enceladus.CR3BP.SLeadingSneakerToe; famName = 'Saturn_Enceladus.CR3BP.SLeadingSneakerToe';
% myPO_ICs = PO_ICs.Saturn_Enceladus.CR3BP.SLeadingSneakerHeel; famName = 'Saturn_Enceladus.CR3BP.SLeadingSneakerHeel';
% myPO_ICs = PO_ICs.Saturn_Enceladus.CR3BP.unknown1; famName = 'Saturn_Enceladus.CR3BP.unknown1';


% myPO_ICs = PO_ICs.Saturn_Enceladus.CR3BP_J2pJ4pJ6pJ2s.L2_Lyapunov; famName = 'Saturn_Enceladus.CR3BP_J2pJ4pJ6pJ2s.L2_Lyapunov';
% myPO_ICs = PO_ICs.Saturn_Enceladus.CR3BP_J2pJ4pJ6pJ2s.L2_Vertical; famName = 'Saturn_Enceladus.CR3BP_J2pJ4pJ6pJ2s.L2_Vertical';
% myPO_ICs = PO_ICs.Saturn_Enceladus.CR3BP_J2pJ4pJ6pJ2s.L2_SHalo; famName = 'Saturn_Enceladus.CR3BP_J2pJ4pJ6pJ2s.L2_SHalo';

% --------------------------
% Saturn-Titan
% --------------------------
% myPO_ICs = PO_ICs.Saturn_Titan.CR3BP.L2_Lyapunov; famName = 'Saturn_Titan.CR3BP.L2_Lyapunov';
% myPO_ICs = PO_ICs.Saturn_Titan.CR3BP.L2_Vertical; famName = 'Saturn_Titan.CR3BP.L2_Vertical';

% --------------------------
% Uranus-Cordelia
% --------------------------
% myPO_ICs = PO_ICs.Uranus_Cordelia.CR3BP.L2_Lyapunov; famName = 'Uranus_Cordelia.CR3BP.L2_Lyapunov';

% --------------------------
% Uranus-Ophelia
% --------------------------
% myPO_ICs = PO_ICs.Uranus_Ophelia.CR3BP.L1_Lyapunov; famName = 'Uranus_Ophelia.CR3BP.L1_Lyapunov';

% --------------------------
% Neptune - Triton
% --------------------------
% myPO_ICs = PO_ICs.Neptune_Triton.CR3BP.L2_Lyapunov; famName = 'Neptune_Triton.CR3BP.L2_Lyapunov';
% myPO_ICs = PO_ICs.Neptune_Triton.CR3BP.L2_Vertical; famName = 'Neptune_Triton.CR3BP.L2_Vertical';
% myPO_ICs = PO_ICs.Neptune_Triton.CR3BP.L2_SHalo; famName = 'Neptune_Triton.CR3BP.L2_SHalo';

% --------------------------
% Set primary & secondary
% --------------------------
[primary, secondary] = assignPrimaryAndSecondary_CR3BP(famName, bodies);

% -------------------------------------------------
%%% System
% -------------------------------------------------
%%% Normalizing constants
rNorm = secondary.a;         % n <-> km
tNorm = 1/secondary.meanMot; % n <-> sec
vNorm = rNorm / tNorm;       % n <-> km/sec

%%% Shortcut variables
mu   = secondary.MR;
R2_n = secondary.R_n;

%%% Setting parameters
prms.u    = mu;
prms.R2   = R2_n; 
prms.R1   = primary.R / rNorm;


%%% Equillibrium Points
if contains(famName,'.CR3BP.')
    prms.n = 1;
    rLPs_n = EquilibriumPoints(mu, prms.n);
elseif contains(famName,'.CR3BP_J2pJ4pJ6pJ2s.')
    prms.J2p = primary.J2; 
    prms.J4p = primary.J4; 
    prms.J6p = primary.J6; 
    prms.J2s = secondary.J2;
    rLPs_n = collinearEquilibriumPoints_ZH(prms);
    warning('Need prms.n')
end

% ------------------------------------------------- 
%%% Integration Options
% -------------------------------------------------
tol = 1e-13;
options = odeset('RelTol',tol,'AbsTol',tol);
options_XZStop = odeset('Event',@event_yEqualsZeroPastL2,'RelTol',tol,'AbsTol',tol);
t0 = 0;
% ------------------------------------------------- 
%%% Other options
% -------------------------------------------------
%%% Maximum number of POs to find in family, regardless of step size
n_POs_max = 10;

%%% Lagrange point? 
if contains(famName,'.L1_')
    LP = 1;
elseif contains(famName,'.L2_')
    LP = 2;
end


%%% Initial size to next PO in family. This value will dynamically change
%%% based on performance as the algorithm runs
% ds_PO = 1e-2;
% ds_PO = 7e-3; %  
ds_PO = 1e-3; % Europa 
% ds_PO = 1e-4; % Enceladus sneaker
% ds_PO = 1e-5;
% ds_PO = 1e-6; % Enceladus 
% ds_PO = 1e-7;
% ds_PO = 1e-8;
% ds_PO = 1e-9;

%%% Error tolerance for constraint vector in multiple shooter

error_tol = 1e-12; 
% error_tol = 1e-10;
% error_tol = 1e-9;

%%% Continue family in or out
familyDirection = 1; % 1 or -1

%%% Number of nodes for multiple shooter
n_Nodes = 4; % 3

%%% Artificial scalar for increasing or decreasing the step size for 
% guessing at the state of the 2nd family member. 1 has no effect, 2 and 
% 3 can be helpful in certain regions
stepSize = 1;

%%% Maximum number of multiple-shooter iterations for any given family
%%% member before kicking out of the loop and adjusting the tuning
%%% parameters
iterMax = 500;

%%% Maximum number of times the alorithm can try to relocate the family if
%%% it jumps to a different family/equilibrium 
lostFamilyCounter_max = 5;

% ------------------------------------------------- 
%%% Preparing to find family
% -------------------------------------------------
%%% Preallocating space for outputs
Fs                      = NaN(n_POs_max, 6*n_Nodes + 1);
POs                     = NaN(n_POs_max,7);
stabilityIndices        = NaN(n_POs_max,2);
jacobiConstants         = NaN(n_POs_max,1);
L2ExcessVelocities_mps  = NaN(n_POs_max,1);
landingVelocities_mps   = NaN(n_POs_max,1);

%%% Storing first PO
POs(1,:) = [myPO_ICs(1:6)', myPO_ICs(7)];

%%% Setting chosen ICs as first guess
X0_guess_n = myPO_ICs(1:6);
T_guess_n_new  = myPO_ICs(7);

%%% Storing first indices
JC_PO_i                        = getJacobiConstant_ZH(POs(1,1:6), prms);
L2_FlythroughVelocity_PO_i_mps = JC_2_L2FlyoverVelocity(JC_PO_i,prms,rLPs_n(2,:),vNorm);
jacobiConstants(1)             = JC_PO_i;
L2ExcessVelocities_mps(1)      = L2_FlythroughVelocity_PO_i_mps;
landingVelocities_mps(1)       = JC_2_approxLandingVelocity(JC_PO_i, prms, vNorm);

%%% Initialize STM
stm0        = eye(6);
stm0_colVec = reshape(stm0,36,1);

%%% Initializing a counter for helping to identify when a family has been
%%% lost track of
lostFamilyCounter = 0;

% ------------------------------------------------- 
%%% Integrate reference trajectory and divide into nodes for multiple
%%% shooter
% -------------------------------------------------
%%% Integrate reference PO
if contains(famName,'.CR3BP.')
    [Tref_n, Xref_n] = ode113(@Int_CR3BnSTM, [0, myPO_ICs(7)], [myPO_ICs(1:6); stm0_colVec], options, prms);
elseif contains(famName,'.CR3BP_J2pJ4pJ6pJ2s.')
    [Tref_n, Xref_n] = ode113(@Int_CR3BnSTM_J2pJ4pJ6pJ2s, [0, myPO_ICs(7)], [myPO_ICs(1:6); stm0_colVec], options, prms);
end

%%% Integrate again to discretize into nodes evenly spaced in time
if contains(famName,'.CR3BP.')
%     [T_nodes, X_nodes] = ode113(@Int_CR3BnSTM, linspace(0, myPO_ICs(7),n_Nodes+1), [myPO_ICs(1:6); stm0_colVec], options, prms);
    [T_nodes, X_nodes] = get_nodes(myPO_ICs(1:6), [0, myPO_ICs(7)], n_Nodes+1, @Int_CR3Bn, options, prms);
elseif contains(famName,'.CR3BP_J2pJ4pJ6pJ2s.')
    [T_nodes, X_nodes] = ode113(@Int_CR3BnSTM_J2pJ4pJ6pJ2s, linspace(0, myPO_ICs(7),n_Nodes+1), [myPO_ICs(1:6); stm0_colVec], options, prms);
end

%%% Setting initial nodes and node-time
X_nodes = X_nodes(1:n_Nodes,1:6);
T_PO_i  = myPO_ICs(7);

%%% Create first free-variable column vector
F_new   = X_nodes';
F_new   = F_new(:);
F_new   = [F_new; T_PO_i];
Fs(1,:) = F_new';

%%% Plotting reference PO
if plot_current_PO == 1
    figure(1); hold all
    if contains(famName,'.CR3BP.')
        plot3(Xref_n(:,1),Xref_n(:,2),Xref_n(:,3),'r','linewidth',2)
    elseif contains(famName,'.CR3BP_J2pJ4pJ6pJ2s.')
        plot3(Xref_n(:,1),Xref_n(:,2),Xref_n(:,3),'b','linewidth',2)
    end
%     PlotBoi3('$X_n$','$Y_n$','$Z_n$',20,'LaTex')
    PlotBoi3_CR3Bn(20)
    view(-10,50)
    
end % plot_current_PO


%%% Calculate stability indices of reference orbit
stm_tf_t0                           = reshape(Xref_n(end,7:42),6,6);
monodromy                           = stm_tf_t0;
[eigenVectors_new, eigenValues_new] = eig(monodromy);
[S1, S2]                            = getStabilityIndices(diag(eigenValues_new));
stabilityIndices(1,:)               = [S1, S2];

% ------------------------------------------------- 
%%% Enter loop to populate family of POs
% -------------------------------------------------
PO_i = 1;
while PO_i <= n_POs_max
    % --------------------------
    % Setup for new PO
    % --------------------------
    %%% Set free variable vector
    F_new = Fs(PO_i,:)';
    
%     if PO_i == 2
%         return
%     end
    
    % -------------------------------------------------
    % Run pseudo-arclength continuation with multiple shooting to find the
    % next orbit
    % -------------------------------------------------
    if PO_i == 1
        nullVecDF = zeros(size(F_new));
    end
    
    if contains(famName,'.CR3BP.')
        [F_new, constraint_error, DF, counter] = pseudoArclengthContinuation_multShoot_stateContinuity(...
            error_tol, iterMax, n_Nodes, F_new, @Int_CR3BnSTM, options, prms, PO_i, stepSize, nullVecDF, ds_PO);
    elseif contains(famName,'.CR3BP_J2pJ4pJ6pJ2s.')
        [F_new, constraint_error, DF, counter] = pseudoArclengthContinuation_multShoot_stateContinuity(...
            error_tol, iterMax, n_Nodes, F_new, @Int_CR3BnSTM_J2pJ4pJ6pJ2s, options, prms, PO_i, stepSize, nullVecDF, ds_PO);
    end
    
    
    if print_PO_index == 1
        fprintf('PO_i = %1d ... iterations: %1.3d ... ds_PO = %1.2e ... Elapsed time: %1.1f seconds\n',PO_i, counter, ds_PO, toc(ticWhole))
    end % print_PO_index
    
    % --------------------------
    % Checks
    % --------------------------
    
    
    % --------------------------
    % If checks are cleared and solution meets error tolerance, then
    % compute null space of DF, adjust step-size if necessary, and
    % store converged F
    % --------------------------
    if (constraint_error < error_tol)
        %%% Compute null space of DF
        if PO_i == 1
            nullVecDF = null(DF).*familyDirection;
        
        %%% If PO_i > 1, make sure null space hasn't flipped sign
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
   
        %%% Check for family continuity (no jumps to other equilibria)
        percentChangePosition = norm(F_new(1:3) - Fs(PO_i,1:3)')/norm(Fs(PO_i,1:3)');
        if (percentChangePosition > 0.1) || (norm(F_new(1:3) - rLPs_n(LP,:)') < 1e-14)
            if lostFamilyCounter < lostFamilyCounter_max
                %%% Count the lost family iteration
                lostFamilyCounter = lostFamilyCounter + 1;

                %%% Decrease step size and try to find family again
                ds_PO = ds_PO * 0.9;
                continue
            
            %%% If algorithm repeatedly fails to continue family, end
            %%% the search
            elseif lostFamilyCounter >= 5
                warning('Lost track of family ... It possibly ended')
                break
            end
        end
        
        %%% Reset counter for when the algorithm loses the family
        lostFamilyCounter = 0;
        
        %%% Store converged F
        Fs(PO_i+1,:) = F_new';
        POs(PO_i+1,:) = [F_new(1:6)', F_new(end)];
        
        % --------------------------
        % Integrate the converged trajectory and store the stability
        % indices, the energy, and the L2 flythrough velocity
        % --------------------------
        %%% Integrate new PO
        if contains(famName,'.CR3BP.')
            [~, X_new] = ode113(@Int_CR3BnSTM, [0, F_new(end)], [F_new(1:6); stm0_colVec], options, prms);
        elseif contains(famName,'.CR3BP_J2pJ4pJ6pJ2s.')
            [~, X_new] = ode113(@Int_CR3BnSTM_J2pJ4pJ6pJ2s, [0, F_new(end)], [F_new(1:6); stm0_colVec], options, prms);
        end
        
        %%% Plotting PO
        if plot_current_PO
            if mod(PO_i,plotSkip) == 0
                if contains(famName,'.CR3BP.')
                    plot3(X_new(:,1),X_new(:,2),X_new(:,3),'r')
                elseif contains(famName,'.CR3BP_J2pJ4pJ6pJ2s.')
                    plot3(X_new(:,1),X_new(:,2),X_new(:,3),'b')
                end
                drawnow
            end
        end
        
        %%% Stability indices of new PO
        stm_tf_t0                           = reshape(X_new(end,7:42),6,6);
        monodromy                           = stm_tf_t0;
        [eigenVectors_new, eigenValues_new] = eig(monodromy);
        [S1, S2]                            = getStabilityIndices(diag(eigenValues_new));
        stabilityIndices(PO_i+1,:)          = [S1, S2];
        
        %%% Energy of new PO
        JC_PO_new                     = getJacobiConstant_ZH(X_new(1,1:6), prms);
        L2_FlythroughVelocity_mps_new = JC_2_L2FlyoverVelocity(JC_PO_new, prms, rLPs_n(2,:), vNorm);
        
        jacobiConstants(PO_i+1)        = JC_PO_new;
        L2ExcessVelocities_mps(PO_i+1) = L2_FlythroughVelocity_mps_new;
        landingVelocities_mps(PO_i+1)  = JC_2_approxLandingVelocity(JC_PO_new, prms, vNorm);
        % --------------------------
        % Adjust step size
        % --------------------------
        %%% Make corrections to dS
%         if counter <= 5
%         if isequal(secondary.name,'europa')
%             if counter < 100
%                 ds_PO = ds_PO * 1.2;
%             end
%         end
%         if isequal(secondary.name,'enceladus')
%             if counter < 100
%                 ds_PO = ds_PO * 1.1; 
%             end
%         end
        
        if counter <= 100
%             if isequal(secondary.name,'enceladus')
%                 if PO_i == 66
%                     ds_PO = 9e-5;
% %                     if PO_i > 137
% %                         ds_PO = 0.00001;
% %                     end
%                 else
%                     ds_PO = 1.05 * ds_PO;
%                 end
%                 if PO_i == 192
%                     break
%                 end
            if isequal(secondary.name,'titan')
                ds_PO = 1.2*ds_PO;
%             elseif isequal(secondary.name,'cordelia')
%                 if PO_i > 138
%                     ds_PO = 0.1;
%                 end
            elseif isequal(secondary.name,'cordelia')
                if PO_i > 115
                    ds_PO = 5e-3;
                end
                if PO_i > 310
                    ds_PO = 3e-4;
                end
                
            elseif isequal(secondary.name,'ophelia')
                if PO_i > 118
                    ds_PO = 5e-3;
                end
                
            else % if none of these bodies
                ds_PO = 1.2*ds_PO;
            end
        end
        
        %%% Update PO_i
        PO_i = PO_i + 1;
    else
        if PO_i == 1
            warning('Failed to converge on a 2nd PO')
        end
        
        %%% Lower the dS
        if isequal(secondary.name,'titan')
            ds_PO = 0.7*ds_PO;
        else % if none of these bodies
            ds_PO = 0.5*ds_PO;
        end
    end
    
    
end % while PO_i < n_POs_max

%% =======================================================================
%%% Clean up the data
% ========================================================================
% -------------------------------------------------
%%% Get rid of extra 'NaN' states
% -------------------------------------------------
%%% Grab all non-NaN data
nonNaN_ind = ~isnan(POs(:,1));
POs_trimmed = POs(nonNaN_ind,:);

% -------------------------------------------------
%%% Integrating solution to x-y plane crossing for state simplicity
% -------------------------------------------------
POs_y0EqualsZero = NaN(size(POs_trimmed));

for kk = 1:size(POs_trimmed,1)
    %%% Integrate
    if contains(famName,'.CR3BP.')
        [~, X_BCR_n_XY, time_event, X_event, index_event] = ...
            ode113(@Int_CR3Bn, [0, POs_trimmed(kk,end)], POs_trimmed(kk,1:6)', options_XZStop, prms);
    elseif contains(famName,'.CR3BP_J2pJ4pJ6pJ2s.')
        [~, X_BCR_n_XY, time_event, X_event, index_event] = ...
            ode113(@Int_CR3Bn_J2pJ4pJ6pJ2s, [0, POs_trimmed(kk,end)], POs_trimmed(kk,1:6)', options_XZStop, prms);
    end
    
    %%% Attempt, but if trajectory doesn't cross x-y plane, then forget it
    if isempty(X_event) == 0
        POs_y0EqualsZero(kk,:) = [X_event(end,1:6),POs_trimmed(kk,end)];
    else
        POs_y0EqualsZero(kk,:) = [POs_trimmed(kk,1:6),POs_trimmed(kk,end)];
    end
end
    

% ========================================================================
%%% Plots/Studies
% ========================================================================
% -------------------------------------------------
%%% Plot Stability Indices
% -------------------------------------------------
%%% Find stability indices
s1_ind = ~isnan(stabilityIndices(:,1));
s2_ind = ~isnan(stabilityIndices(:,2));
S1 = stabilityIndices(s1_ind,1);
S2 = stabilityIndices(s2_ind,2);
    
if plot_stability == 1
    JC_temp1 = jacobiConstants(s1_ind);
    JC_temp2 = jacobiConstants(s2_ind);
    
    figure('position',[209 322 948 302])
    subplot(1,2,1); hold all
    p1 = plot(JC_temp1, abs(S1),'o','markeredgecolor',colors.blue,'markerfacecolor',colors.ltblue);
    p2 = plot(JC_temp2, abs(S2),'o','markeredgecolor',colors.red,'markerfacecolor',colors.ltred);
    plot(unique([min([JC_temp2, JC_temp1]) max([JC_temp2, JC_temp1])]),[2 2],'k','linewidth',1)
    PlotBoi2('Jacobi Constant','Stability Indices',18,'LaTex')
    legend([p1 p2],'S_1','S_2')
    xlim(unique([min([JC_temp2, JC_temp1]) max([JC_temp2, JC_temp1])]))
    
    subplot(1,2,2); hold all
    p1 = plot(JC_temp1, abs(S1),'o','markeredgecolor',colors.blue,'markerfacecolor',colors.ltblue);
    p2 = plot(JC_temp2, abs(S2),'o','markeredgecolor',colors.red,'markerfacecolor',colors.ltred);
    plot(unique([min([JC_temp2, JC_temp1]) max([JC_temp2, JC_temp1])]),[2 2],'k','linewidth',1)
    PlotBoi2('Jacobi Constant','',18,'LaTex')
    ylim([1.9 2.1])
    xlim(unique([min([JC_temp2, JC_temp1]) max([JC_temp2, JC_temp1])]))
end

% -------------------------------------------------
%%% Plot Energy
% -------------------------------------------------
if plot_energy
%     figure; hold all
%     yyaxis left
%     plot(jacobiConstants,'o','markeredgecolor',colors.blue,'markerfacecolor',colors.ltblue)
%     PlotBoi2('','Jacobi Constant',18,'LaTex')
%     yyaxis right
%     plot(L2ExcessVelocities_mps,'o','markeredgecolor',colors.red,'markerfacecolor',colors.ltred)
%     PlotBoi2('PO Index','$L_2$ Excess Velocity, $m/s$',18,'LaTex')
    
    figure; hold all
    yyaxis left
    plot(jacobiConstants,'o','markeredgecolor',colors.blue,'markerfacecolor',colors.ltblue)
    PlotBoi2('','Jacobi Constant',18,'LaTex')
    yyaxis right
    plot(landingVelocities_mps,'o','markeredgecolor',colors.red,'markerfacecolor',colors.ltred)
    PlotBoi2('PO Index','+$x$ Landing Velocity, $m/s$',18,'LaTex')
    
end

% ========================================================================
%%% Saving Data
% ========================================================================
if save_PO_X0_database == 1
    % -------------------------------------------------
    %%% Preparing Save File
    % -------------------------------------------------
    fileAlreadyExists = 1;
    fileVersion       = 0;
    famName_new       = famName;
    while fileAlreadyExists == 1
        if isfile([savePath,famName_new,'.txt']) == 1
            fileVersion = fileVersion + 1;
            famName_new = sprintf('%s_%1d',famName,fileVersion);
        else
            fileAlreadyExists = 0;
        end
    end
    fileName = [savePath,famName_new,'.txt'];
    
    %%% Open File
    datafile = fopen(fileName,'wt');
    
    % -------------------------------------------------
    %%% Writing data
    % -------------------------------------------------
    %%% Write header
    headerString = ['x0_n,y0_n,z0_n,xd0_n,yd0_n,zd0_n,Tp,JC,L2ExcessVelocity_mps,landingVelocity_mps,stabilityIndex1,stabilityIndex2...,',...
                    sprintf('error_tol=%1.0e\n',error_tol)];
    fprintf(datafile,headerString);
    
    %%% Write data
    dataStart = 1;
    dataStop = size(POs_y0EqualsZero,1);
    
    for kk = dataStart:dataStop
        fprintf(datafile,'%1.16f,%1.16f,%1.16f,%1.16f,%1.16f,%1.16f,%1.16f,%1.10f,%1.3f,%1.3f,%1.4f,%1.4f\n',...
            POs_y0EqualsZero(kk,1), POs_y0EqualsZero(kk,2), POs_y0EqualsZero(kk,3),...
            POs_y0EqualsZero(kk,4), POs_y0EqualsZero(kk,5), POs_y0EqualsZero(kk,6),...
            POs_y0EqualsZero(kk,end), jacobiConstants(kk), L2ExcessVelocities_mps(kk),...
            landingVelocities_mps(kk), stabilityIndices(kk,1), stabilityIndices(kk,2));
            
    end

    %%% Close file
    fclose(datafile);
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
fprintf('Elapsed time: %1.4f seconds\n',tocWhole)








