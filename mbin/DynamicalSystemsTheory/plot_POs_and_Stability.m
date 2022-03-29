% ========================================================================
%%% Description
% ========================================================================
% For plotting families of POs and describing the bifurcations that occur
% within them. For plotting, designed to integrate the trajectories in
% parallel, but the 'parfor' can be switched to 'for' if small, fast
% batches are desired

% Created: 12/16/20
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
plot_family        = true;
plot_energy        = true;
plot_stability     = true;
    check_for_false_Stability_Indices = false; % Warning, this takes a while
print_bifurcations = false;
    plot_BrouckeDiagram = false;

check_Tp_monotonic = false;

study_missionUtility = false; % Warning, this takes a while

% ========================================================================
%%% Setup
% ========================================================================
% -------------------------------------------------
%%% Personal options
% -------------------------------------------------


% -------------------------------------------------
%%% Choose data
% -------------------------------------------------
% family = 'Earth_Moon.CR3BP.L2_SHalo.txt';
% family = 'Earth_Moon.CR3BP.L2_L_1T_1P2.txt';
% family = 'Earth_Moon.CR3BP.L2_L_1T_2P2.txt';
% family = 'Earth_Moon.CR3BP.L2_L_1T_3P2.txt';
% family = 'Earth_Moon.CR3BP.L2_L_1T_4P2.txt';
% family = 'Earth_Moon.CR3BP.Apollo.txt';


% family = 'Jupiter_Europa.CR3BP.L2_Lyapunov.txt';
% family = 'Jupiter_Europa.CR3BP.L2_NHalo.txt';
% family = 'Jupiter_Europa.CR3BP.L2_SHalo.txt';
% family = 'Jupiter_Europa.CR3BP.L2_L_1T_1P4.txt'; 
% family = 'Jupiter_Europa.CR3BP.L2_L_1T_1P3.txt'; 
% family = 'Jupiter_Europa.CR3BP.L2_L_1T_1P2.txt'; 
% family = 'Jupiter_Europa.CR3BP.L2_L_1T_2P2.txt'; 
% family = 'Jupiter_Europa.CR3BP.L2_L_1T_3P2.txt'; 
% family = 'Jupiter_Europa.CR3BP.L2_L_1T_4P2.txt'; % Butterfly
% % % % % family = 'Jupiter_Europa.CR3BP.L2_L_1T_4P2_1P4.txt';
% family = 'Jupiter_Europa.CR3BP.L2_L_1T_4P2_1P3.txt'; % the end of this is the 2P2 bifurcation from LoPO_2P2_1T (so this is also LoPO_2P2_1T_2P2)
% family = 'Jupiter_Europa.CR3BP.L2_L_1T_4P2_1P2.txt';
% family = 'Jupiter_Europa.CR3BP.L2_L_1T_4P2_2P2.txt'; 
% % % family = 'Jupiter_Europa.CR3BP.L2_L_1T_4P2_2P3.txt'; 
% family = 'Jupiter_Europa.CR3BP.L2_L_1T_4P2_2P4.txt'; 
% family = 'Jupiter_Europa.CR3BP.L2_L_2T.txt'; % Full Axial family (Also L2_V_T)
% family = 'Jupiter_Europa.CR3BP.L2_L_1P3.txt'; 
% family = 'Jupiter_Europa.CR3BP.L2_L_1P2.txt'; 
% family = 'Jupiter_Europa.CR3BP.L2_Vertical.txt';

% family = 'Jupiter_Europa.CR3BP.DRO.txt'; 
% family = 'Jupiter_Europa.CR3BP.DRO_1P4.txt'; 
% family = 'Jupiter_Europa.CR3BP.DRO_1P4_1P4.txt'; 
% family = 'Jupiter_Europa.CR3BP.DRO_1P4_1P4_1T.txt'; 
% family = 'Jupiter_Europa.CR3BP.DRO_1P3.txt'; 
% family = 'Jupiter_Europa.CR3BP.DRO_1P3_1P3.txt'; 
% family = 'Jupiter_Europa.CR3BP.DRO_2P3.txt'; 
% family = 'Jupiter_Europa.CR3BP.DRO_2P3_1P2.txt'; 
% family = 'Jupiter_Europa.CR3BP.DRO_2P3_2P2.txt'; 
% family = 'Jupiter_Europa.CR3BP.DRO_2P4.txt'; 
% family = 'Jupiter_Europa.CR3BP.DRO_2P4_1T.txt'; 
% family = 'Jupiter_Europa.CR3BP.DRO_2P4_2T.txt'; 

family = 'Jupiter_Europa.CR3BP.DPO.txt'; 
% family = 'Jupiter_Europa.CR3BP.DPO_1P4.txt'; 
% family = 'Jupiter_Europa.CR3BP.DPO_1P4_1T.txt'; 
% family = 'Jupiter_Europa.CR3BP.DPO_1P3.txt'; 
% family = 'Jupiter_Europa.CR3BP.DPO_1P2.txt'; 
% family = 'Jupiter_Europa.CR3BP.DPO_2P3.txt'; 
% family = 'Jupiter_Europa.CR3BP.DPO_2P2.txt'; 
% family = 'Jupiter_Europa.CR3BP.DPO_3P2.txt'; 
% family = 'Jupiter_Europa.CR3BP.DPO_2P4.txt'; 
% family = 'Jupiter_Europa.CR3BP.DPO_3P3.txt'; 
% family = 'Jupiter_Europa.CR3BP.DPO_3P4.txt'; 
% family = 'Jupiter_Europa.CR3BP.DPO_3P4_symmetricAboutXb2.txt'; 
% family = 'Jupiter_Europa.CR3BP.DPO_4P4.txt'; 
% family = 'Jupiter_Europa.CR3BP.DPO_4P3.txt'; 
% family = 'Jupiter_Europa.CR3BP.DPO_4P2.txt';  % THIS IS DPO_3P2 FROM THE OTHER SIDE
% family = 'Jupiter_Europa.CR3BP.DPO_5P2.txt';  % THIS IS LoPO_3P2 FROM THE OTHER SIDE
% family = 'Jupiter_Europa.CR3BP.DPO_2T.txt'; 

% family = 'Jupiter_Europa.CR3BP.LoPO.txt'; 
% family = 'Jupiter_Europa.CR3BP.LoPO_1P4.txt'; 
% family = 'Jupiter_Europa.CR3BP.LoPO_1P3.txt'; 
% family = 'Jupiter_Europa.CR3BP.LoPO_1P2.txt'; 
% family = 'Jupiter_Europa.CR3BP.LoPO_2P2.txt'; 
% family = 'Jupiter_Europa.CR3BP.LoPO_2P2_1T.txt'; 
% family = 'Jupiter_Europa.CR3BP.LoPO_2P4.txt'; 
% family = 'Jupiter_Europa.CR3BP.LoPO_3P3.txt'; 
% family = 'Jupiter_Europa.CR3BP.LoPO_3P4.txt'; 

% family = 'Jupiter_Europa.CR3BP.Hg1.txt'; 
% family = 'Jupiter_Europa.CR3BP.Hg1_1P2.txt'; 
% family = 'Jupiter_Europa.CR3BP.Hg1_2P2.txt'; 
% family = 'Jupiter_Europa.CR3BP.Hg1_1P3.txt'; 
% family = 'Jupiter_Europa.CR3BP.Hg1_1T.txt'; 
% family = 'Jupiter_Europa.CR3BP.Hg1_2T.txt'; 

% family = 'Jupiter_Europa.CR3BP.Hg2.txt'; 
% family = 'Jupiter_Europa.CR3BP.Hg2_1P2.txt';
% family = 'Jupiter_Europa.CR3BP.Hg2_2P2.txt';
% family = 'Jupiter_Europa.CR3BP.Hg2_1P3.txt';
% family = 'Jupiter_Europa.CR3BP.Hg2_2T.txt'; 
% family = 'Jupiter_Europa.CR3BP.Hg2_3T.txt'; 
% family = 'Jupiter_Europa.CR3BP.Hg2_5T.txt';

% family = 'Jupiter_Europa.CR3BP.Se7.txt'; 
% family = 'Jupiter_Europa.CR3BP.Se7_2P2.txt'; 
% family = 'Jupiter_Europa.CR3BP.Se7_2T.txt'; 
% family = 'Jupiter_Europa.CR3BP.Se7_3T.txt'; 
% family = 'Jupiter_Europa.CR3BP.Se7_5P2.txt'; 
% family = 'Jupiter_Europa.CR3BP.Se7_6P2.txt'; 
% family = 'Jupiter_Europa.CR3BP.Se7_5T.txt'; 
% family = 'Jupiter_Europa.CR3BP.Se7_6T.txt'; 

% family = 'Jupiter_Europa.CR3BP.Hg1_full?.txt'; 
% family = 'Jupiter_Europa.CR3BP.Hg3.txt'; 



% 50 mps l2 flythrough velocity => 3.003595836097908
% --------------------------
% Actually load family
% --------------------------
%%% Path from mbin to data
dataPathFromMBin = '/Data/InitialConditions/PO_Families/';

%%% PO data file
PO_datafile = [mbinPath, dataPathFromMBin, family];
% -------------------------------------------------
%%% Set up parameters
% -------------------------------------------------
%%% Set primary and secondary bodies
[primary, secondary] = assignPrimaryAndSecondary_CR3BP(family, bodies);

%%% Normalizing constants
[rNorm, tNorm, vNorm] = cr3bp_norms(primary, secondary, bodies.constants.G);

%%% prms for integration
prms.u  = secondary.MR;
prms.R2 = secondary.R_n;

if contains(family,'.CR3BP.')
    prms.n    = 1;
    rLPs_n = EquilibriumPoints(prms.u, prms.n);
end

%%% Integration options
tol = 1e-13;
options = odeset('RelTol',tol,'AbsTol',tol);

% -------------------------------------------------
%%% Load data
% -------------------------------------------------
%%% Load the data file
PO_data = dlmread(PO_datafile,',',1,0);

%%% Grab header line
fid = fopen(PO_datafile, 'rt');  %the 't' is important!
header = fgetl(fid);
fclose(fid);

%%% Number of ICs
n_POs = size(PO_data,1);

%%% Indices for plotting purposes
PO_indices = linspace(1, n_POs, n_POs);
% --------------------------
% Create column specifiers
% --------------------------
PO_header_2020 = 'x0,y0,z0,xd0,yd0,zd0,Tp,JC,stabilityIndex1,stabilityIndex2,alpha,beta,impactFlag,error';

if contains(header, PO_header_2020)
    c_x0 = 1;   c_y0 = 2;   c_z0 = 3;
    c_xd0 = 4;  c_yd0 = 5;  c_zd0 = 6;
    c_Tp = 7;   c_JC = 8;   c_S1 = 9;   c_S2 = 10;
    c_alpha = 11;   c_beta = 12;    c_impactFlag = 13;
    c_error = 14;
end

% ========================================================================
%%% Stability Analysis
% ========================================================================
if plot_energy
% % %     %%% Calculate approximate landing velocities
% % %     [landingVelocities_mps] = JC_2_approxLandingVelocity(PO_data(:, c_JC), prms, vNorm);
% % %     
% % %     %%% Plot Jacobi constant and corresponding approximate landing velocity vs Tp
% % %     
% % %     
% % %     figure; hold all
% % %     yyaxis left
% % %     plot3(PO_data(:, c_Tp), PO_data(:, c_JC), PO_indices,'o','markeredgecolor',colors.blue,'markerfacecolor',colors.ltblue)
% % %     PlotBoi2('','Jacobi Constant',26,'LaTex')
% % %     yyaxis right
% % %     plot3(PO_data(:, c_Tp), landingVelocities_mps, PO_indices,'o','markeredgecolor',colors.red,'markerfacecolor',colors.ltred)
% % %     PlotBoi2('$T_P$','+$x$ Landing Velocity, $m/s$',26,'LaTex')
   
    
    %%% Plot Jacobi constant and Tp
    figure; hold all
    plot3(PO_data(:, c_Tp), PO_data(:, c_JC), PO_indices,'o','markeredgecolor',colors.blue,'markerfacecolor',colors.ltblue)
    PlotBoi2('$T_P$','Jacobi Constant',26,'LaTex')
end

if plot_stability    
    figure('position',[209 322 948 302])
    subplot(1,2,1); hold all
    p1 = plot3(PO_data(:, c_JC), abs(PO_data(:, c_S1)), PO_indices,'o','markeredgecolor',colors.blue,'markerfacecolor',colors.ltblue);
    p2 = plot3(PO_data(:, c_JC), abs(PO_data(:, c_S2)), PO_indices,'o','markeredgecolor',colors.red,'markerfacecolor',colors.ltred);
    plot3(unique([min(PO_data(:, c_JC)) max(PO_data(:, c_JC))]),[2 2], [1 1],'k','linewidth',1)
    PlotBoi2('Jacobi Constant','Stability Indices',26,'LaTex')
    legend([p1 p2],'S_1','S_2','FontSize',14)
    xlim(unique([min(PO_data(:, c_JC)) max(PO_data(:, c_JC))]))
    view(0,90)
    
    subplot(1,2,2); hold all
    p1 = plot3(PO_data(:, c_JC), abs(PO_data(:, c_S1)), PO_indices,'o','markeredgecolor',colors.blue,'markerfacecolor',colors.ltblue);
    p2 = plot3(PO_data(:, c_JC), abs(PO_data(:, c_S2)), PO_indices,'o','markeredgecolor',colors.red,'markerfacecolor',colors.ltred);
    plot3(unique([min(PO_data(:, c_JC)) max(PO_data(:, c_JC))]),[2 2], [1 1],'k','linewidth',1)
    PlotBoi2('Jacobi Constant','',26,'LaTex')
    ylim([0 4])
    xlim(unique([min(PO_data(:, c_JC)) max(PO_data(:, c_JC))]))
    view(0,90)
    
    %%% Recalculate stability
    if 1+1==1
%         stm0_colVec = reshape(eye(6),36,1);
%         newStabilityIndices = nan(n_POs,2);
%         for kk = 1:n_POs
%             [~, X_new] = ode113(@Int_CR3BnSTM, [0, PO_data(kk, c_Tp)], [PO_data(kk, c_x0:c_zd0)'; stm0_colVec], options, prms);
%             stm_tf_t0                           = reshape(X_new(end,7:42),6,6);
%             monodromy                           = stm_tf_t0;
%             [eigenVectors_new, eigenValues_new] = eig(monodromy);
%             [S1, S2]                            = getStabilityIndices(diag(eigenValues_new));
%             newStabilityIndices(kk,:)          = [S1, S2];
%         end
    end
    
    %%% Check for false positives of stability (where complex components
    %%% lie off the unit circle so magnitude is greater than one
    if check_for_false_Stability_Indices
        stm0_colVec = reshape(eye(6),36,1);
%         newStabilityIndices = nan(n_POs,2);
        for kk = 1:n_POs
            
            if max(abs([PO_data(kk, c_S1), PO_data(kk, c_S2)])) < 11
            
                [~, X_new] = ode113(@Int_CR3BnSTM, [0, PO_data(kk, c_Tp)], [PO_data(kk, c_x0:c_zd0)'; stm0_colVec], options, prms);
                stm_tf_t0                           = reshape(X_new(end,7:42),6,6);
                monodromy                           = stm_tf_t0;
                [eigenVectors_new, eigenValues_new] = eig(monodromy);
                [S1, S2]                            = getStabilityIndices(diag(eigenValues_new));
                if (max(rowNorm(diag(eigenValues_new))) > max(abs([S1,S2])))
                    warning('Stability indices may not be representative ... May have large imaginary components')
                    kk
                    break
                end
            end
%             newStabilityIndices(kk,:)          = [S1, S2];
        end
    end
        
end




















if print_bifurcations
    [bifurcation_strings] = plot_BrouckeStabilityDiagram(PO_data(:, c_alpha), PO_data(:, c_beta), plot_BrouckeDiagram);
end


% ========================================================================
%%% Check if the time period is monotonic
% ========================================================================
if check_Tp_monotonic
    [isMonotonic, Tp_turnPoints] = checkIfMonotonic(PO_data(:,c_Tp));
    
    fprintf('\n')
    if isMonotonic
        fprintf('Time period IS monotonic\n')
    else
        fprintf('Time period IS NOT monotonic\n')
        display(Tp_turnPoints)
    end
end



% ========================================================================
%%% Check if anything is 
% ========================================================================
if study_missionUtility
    %%% Linewidth of plotted POs
    lw = 2;
    
    color_nearStable = colors.sch.eighties.pink;
    color_stable     = colors.sch.eighties.blue;
    
    
    stability_sums   = zeros(size(PO_data,1),1);
    stable_flags     = zeros(size(PO_data,1),1);
    nearStable_flags = zeros(size(PO_data,1),1);
    
    maximumSumForNearStability = 10;
    
    for kk = 1:size(PO_data,1)
        stability_sums(kk) = abs(PO_data(kk,c_S1)) + abs(PO_data(kk,c_S2));
        
        if (abs(PO_data(kk,c_S1)) < 2) && (abs(PO_data(kk,c_S2)) < 2)
            stable_flags(kk) = 1;
        end
        
        if (stability_sums(kk) <= maximumSumForNearStability) && (max(abs(PO_data(kk,c_S1:c_S2))) >= 2)
            nearStable_flags(kk) = 1;
        end
    end
    stable_flags = logical(stable_flags);
    nearStable_flags = logical(nearStable_flags);
    unstable_flags = ~(stable_flags | nearStable_flags);
    
    fprintf('\n')
    if sum(stable_flags) > 0
        fprintf('Family exhibits stability\n')
    end
    
    if sum(nearStable_flags) > 0
        fprintf('Family exhibits near-stability\n')
    end
    
    if sum(stable_flags + nearStable_flags) == 0
        fprintf('Family does not exhibit stability\n')
    end
    
    thisShouldBeAll0s = stable_flags & nearStable_flags;
    if sum(thisShouldBeAll0s) > 0
        warning('stability flag vectors have a common member somewhere')
    end
    
    fprintf('\n')
    
    if sum(stable_flags + nearStable_flags) > 0
        
        %%% Grab indices of near-stable POs
        nearStable_PO_indices = find(nearStable_flags==1);
        stable_PO_indices = find(stable_flags==1);
        unstable_PO_indices = find(unstable_flags==1);

        %%% Preallocate
        traj_nearStable_POi = cell(length(nearStable_PO_indices),1);
        traj_stable_POi = cell(length(stable_PO_indices),1);

        %%% Initialize STM as column vector to be added with state
        stm0_colVec = reshape(eye(6),36,1);

        parfor PO_i = 1:length(nearStable_PO_indices)
            index = nearStable_PO_indices(PO_i);

            [T_PO, X_PO] = ode113(@Int_CR3BnSTM, [0, PO_data(index,c_Tp)], [PO_data(index,c_x0:c_zd0)'; stm0_colVec], options, prms);

            %%% Store data
            traj_nearStable_POi{PO_i}.X = X_PO;
            traj_nearStable_POi{PO_i}.T = T_PO;
        end

        parfor PO_i = 1:length(stable_PO_indices)
            index = stable_PO_indices(PO_i);

            [T_PO, X_PO] = ode113(@Int_CR3BnSTM, [0, PO_data(index,c_Tp)], [PO_data(index,c_x0:c_zd0)'; stm0_colVec], options, prms);

            %%% Store data
            traj_stable_POi{PO_i}.X = X_PO;
            traj_stable_POi{PO_i}.T = T_PO;
        end

        nearStable_PO_indices_nonImpacting = [];
        stable_PO_indices_nonImpacting = [];
        
        %%% Loop and plot
        figure; hold all
        nearStableImpactFlag = 0;
        for PO_i = 1:length(traj_nearStable_POi)
            r = rowNorm(traj_nearStable_POi{PO_i}.X(:,1:3) - [1-prms.u, 0, 0]);
            if min(r) < prms.R2
                if nearStableImpactFlag == 0
                    fprintf('Near-stable traj impacted body\n')
                    nearStableImpactFlag = 1;
                end
                plot3(traj_nearStable_POi{PO_i}.X(:,1), traj_nearStable_POi{PO_i}.X(:,2), traj_nearStable_POi{PO_i}.X(:,3), 'linewidth', lw, 'color', colors.sch.eighties.yellow);
                continue
            end
            nearStable_PO_indices_nonImpacting = [nearStable_PO_indices_nonImpacting; nearStable_PO_indices(PO_i)];
            p_nearStable = plot3(traj_nearStable_POi{PO_i}.X(:,1), traj_nearStable_POi{PO_i}.X(:,2), traj_nearStable_POi{PO_i}.X(:,3), 'linewidth', lw, 'color', colors.sch.eighties.yellow);
        end
        
        nearStable_PO_indices_impacting = nearStable_PO_indices(~ismember(nearStable_PO_indices, nearStable_PO_indices_nonImpacting));
        
        stableImpactFlag = 0;
        for PO_i = 1:length(traj_stable_POi)
            r = rowNorm(traj_stable_POi{PO_i}.X(:,1:3) - [1-prms.u, 0, 0]);
            if min(r) < prms.R2
                if stableImpactFlag == 0
                    fprintf('Stable traj impacted body\n')
                    stableImpactFlag = 1;
                end
                plot3(traj_stable_POi{PO_i}.X(:,1), traj_stable_POi{PO_i}.X(:,2), traj_stable_POi{PO_i}.X(:,3), 'linewidth', lw, 'color', colors.sch.eighties.purple);
                continue
            end
            stable_PO_indices_nonImpacting = [stable_PO_indices_nonImpacting; stable_PO_indices(PO_i)];
            p_stable = plot3(traj_stable_POi{PO_i}.X(:,1), traj_stable_POi{PO_i}.X(:,2), traj_stable_POi{PO_i}.X(:,3), 'linewidth', lw, 'color', colors.sch.eighties.blue);
        end
        
        stable_PO_indices_impacting = stable_PO_indices(~ismember(stable_PO_indices, stable_PO_indices_nonImpacting));
        
        PlotBoi3_CR3Bn(26)
        axis equal
        plot3(rLPs_n(1:2,1),[0,0],[0,0],'^','markeredgecolor',colors.black, 'markerfacecolor', colors.blue)
        if (isempty(nearStable_PO_indices_nonImpacting)) && (~isempty(stable_PO_indices_nonImpacting))
            legend([ p_stable], 'Stable')
        elseif (~isempty(nearStable_PO_indices_nonImpacting)) && (isempty(stable_PO_indices_nonImpacting))
            legend([p_nearStable], 'Near-Stable')
        elseif (~isempty(nearStable_PO_indices_nonImpacting)) && (~isempty(stable_PO_indices_nonImpacting))
            legend([p_nearStable, p_stable], 'Near-Stable', 'Stable')
        end
        
        
        
    elseif sum(stable_flags + nearStable_flags) == 0
        unstable_PO_indices = linspace(1,size(PO_data,1),size(PO_data,1));
        
        nearStable_PO_indices_impacting    = [];
        nearStable_PO_indices_nonImpacting = [];
        stable_PO_indices_impacting        = [];
        stable_PO_indices_nonImpacting     = [];
    end
    
    
    %%% Get indices of impacting & non-impacting unstable POs
    unstable_PO_indices_impacting = unstable_PO_indices((PO_data(unstable_PO_indices, c_impactFlag) == 1));
    unstable_PO_indices_nonImpacting = unstable_PO_indices((PO_data(unstable_PO_indices, c_impactFlag) == 0));
    
    
    
    %%% Plot Jacobi constant and Tp
    figure; hold all
    PlotBoi2('$T_P$','Jacobi Constant',26,'LaTex')
    
    %%% Plot impacting POs
    if ~isempty(unstable_PO_indices_impacting) % unstable, impacting POs
        plot3(PO_data(unstable_PO_indices_impacting, c_Tp), PO_data(unstable_PO_indices_impacting, c_JC), unstable_PO_indices_impacting,'x','markeredgecolor',colors.ltgrey,'markerfacecolor',colors.ltgrey, 'markersize', 7, 'DisplayName', 'Unstable, impacting')
    end
    if ~isempty(nearStable_PO_indices_impacting) % near-stable, impacting POs
        plot3(PO_data(nearStable_PO_indices_impacting, c_Tp), PO_data(nearStable_PO_indices_impacting, c_JC), nearStable_PO_indices_impacting,'x','markeredgecolor',color_nearStable,'markerfacecolor',color_nearStable, 'markersize', 7, 'DisplayName', 'Near-stable, impacting')
    end
    if ~isempty(stable_PO_indices_impacting) % stable, impacting POs
        plot3(PO_data(stable_PO_indices_impacting, c_Tp), PO_data(stable_PO_indices_impacting, c_JC), stable_PO_indices_impacting,'x','markeredgecolor',color_stable,'markerfacecolor',color_stable, 'markersize', 7, 'DisplayName', 'Stable, impacting')
    end
    
    
    if ~isempty(unstable_PO_indices_nonImpacting) % unstable, non-impacting POs
        plot3(PO_data(unstable_PO_indices_nonImpacting, c_Tp), PO_data(unstable_PO_indices_nonImpacting, c_JC), unstable_PO_indices_nonImpacting,'o','markeredgecolor',colors.ltgrey,'markerfacecolor',colors.ltgrey, 'markersize', 17, 'DisplayName', 'Unstable')
    end
    if ~isempty(nearStable_PO_indices_nonImpacting) % near-stable, non-impacting POs
        plot3(PO_data(nearStable_PO_indices_nonImpacting, c_Tp), PO_data(nearStable_PO_indices_nonImpacting, c_JC), nearStable_PO_indices_nonImpacting,'o','markeredgecolor',color_nearStable,'markerfacecolor',color_nearStable, 'markersize', 17, 'DisplayName', 'Near-stable')
    end
    if ~isempty(stable_PO_indices_nonImpacting) % stable, non-impacting POs
        plot3(PO_data(stable_PO_indices_nonImpacting, c_Tp), PO_data(stable_PO_indices_nonImpacting, c_JC), stable_PO_indices_nonImpacting,'o','markeredgecolor',color_stable,'markerfacecolor',color_stable, 'markersize', 17, 'DisplayName', 'Stable')
    end
    
    legend show 
    hLegend = findobj(gcf, 'Type', 'Legend');
    hLegend.FontSize = 14;
    hLegend.Location = 'best';
    
%     if sum(unstable_flags) > 0
%         plot3(PO_data(unstable_flags, c_Tp), PO_data(unstable_flags, c_JC), unstable_PO_indices,'o','markeredgecolor',colors.ltgrey,'markerfacecolor',colors.ltgrey, 'markersize', 8)
%     end
%     if sum(stable_flags) > 0
%         plot3(PO_data(stable_flags, c_Tp), PO_data(stable_flags, c_JC), stable_PO_indices,'o','markeredgecolor',color_stable,'markerfacecolor',color_stable, 'markersize', 8);
%     end
%     if sum(nearStable_flags) > 0
%         plot3(PO_data(nearStable_flags, c_Tp), PO_data(nearStable_flags, c_JC), nearStable_PO_indices,'o','markeredgecolor',color_nearStable,'markerfacecolor',color_nearStable, 'markersize', 8);
%     end
end


% ========================================================================
%%% Plot POs
% ========================================================================
if plot_family
    %%% Linewidth of plotted POs
    lw = 2;

    %%% Choose range
%     plot_PO_indices = 1:n_POs;
%     plot_PO_indices = 1;
%     plot_PO_indices = n_POs;
%     plot_PO_indices = [1, n_POs];
%     plot_PO_indices = getIndices_spacedByTpJcArclength(PO_data(:,c_Tp), PO_data(:,c_JC), 3);
%     plot_PO_indices = getIndices_spacedByTpJcArclength(PO_data(:,c_Tp), PO_data(:,c_JC), 5);
%     plot_PO_indices = getIndices_spacedByTpJcArclength(PO_data(:,c_Tp), PO_data(:,c_JC), 6);
%     plot_PO_indices = getIndices_spacedByTpJcArclength(PO_data(:,c_Tp), PO_data(:,c_JC), 10);
%     plot_PO_indices = getIndices_spacedByTpJcArclength(PO_data(:,c_Tp), PO_data(:,c_JC), 15);
    plot_PO_indices = getIndices_spacedByTpJcArclength(PO_data(:,c_Tp), PO_data(:,c_JC), 25);
%     plot_PO_indices = getIndices_spacedByTpJcArclength(PO_data(:,c_Tp), PO_data(:,c_JC), 30);
%     plot_PO_indices = getIndices_spacedByTpJcArclength(PO_data(:,c_Tp), PO_data(:,c_JC), 35);
%     plot_PO_indices = getIndices_spacedByTpJcArclength(PO_data(:,c_Tp), PO_data(:,c_JC), 50);
%     plot_PO_indices = getIndices_spacedByTpJcArclength(PO_data(:,c_Tp), PO_data(:,c_JC), 100);
%     plot_PO_indices = getIndices_spacedByTpJcArclength(PO_data(:,c_Tp), PO_data(:,c_JC), 150);
%     plot_PO_indices = getIndices_spacedByTpJcArclength(PO_data(:,c_Tp), PO_data(:,c_JC), 300);

    plot_PO_indices = 50;
%     plot_PO_indices = 405; %                
%     plot_PO_indices = 896:903; 
%     plot_PO_indices = [1, 20, 40, 60, 80, 100];
%     plot_PO_indices = [933];
%     plot_PO_indices = [1, floor(n_POs/2), n_POs];
%     plot_PO_indices = 336;
%     989
%     plot_PO_indices = 1665;
% 
% plot_PO_indices = [30; 60; plot_PO_indices];
% plot_PO_indices = plot_PO_indices(1:end-15);989
% plot_PO_indices = plot_PO_indices(2:end);989
% plot_PO_indices = sort(plot_PO_indices);
    
    %%% optional color spectrum
%     color_spectrum = colors.blue2;
    color_spectrum = colorScale([colors.blue2; colors.mag], length(plot_PO_indices));
%     color_spectrum = colorScale([colors.cyan; colors.mag], length(plot_PO_indices));
%     color_spectrum = colorScale([colors.blue2; colors.red2], length(plot_PO_indices));
%     color_spectrum = colorScale([colors.sch.d3_3(1,:); colors.sch.d3_3(2,:); colors.sch.d3_3(3,:)], length(plot_PO_indices));
%     color_spectrum = colorScale([colors.sch.d5_1(1,:); colors.sch.d5_1(3,:); colors.sch.d5_1(5,:)], length(plot_PO_indices));
%     color_spectrum = colorScale([colors.sch.d5_1(5,:); colors.mag], length(plot_PO_indices));
%     color_spectrum = colorScale([colors.sch.d5_1(5,:); colors.sch.eighties.yellow; colors.sch.eighties.pink], length(plot_PO_indices));
%     color_spectrum = colorScale([colors.sch.d5_1(1,:); colors.sch.d5_1(2,:); colors.sch.d5_1(3,:); colors.sch.d5_1(4,:); colors.sch.d5_1(5,:)], length(plot_PO_indices));
%     color_spectrum = colorScale([colors.sch.eighties.yellow; colors.sch.eighties.red], length(plot_PO_indices));
%     color_spectrum = colorScale([colors.sch.eighties.yellow; colors.sch.eighties.pink; colors.sch.eighties.cyan], length(plot_PO_indices));
%     color_spectrum = colorScale([colors.black; colors.cyan], length(plot_PO_indices));
%     color_spectrum = colorScale(colors.sch.eighties.yellow, length(plot_PO_indices));
%     color_spectrum = colorScale(colors.sch.d4_1, length(plot_PO_indices));
%     color_spectrum = colorScale([colors.red; colors.orange; colors.ylw; colors.grn; colors.blue; colors.purp2], length(plot_PO_indices));
    
    %%% Preallocate
    traj_POi = cell(length(plot_PO_indices),1);
    
    %%% Initialize STM as column vector to be added with state
    stm0_colVec = reshape(eye(6),36,1);
    
%     figure; hold all
    parfor PO_i = 1:length(plot_PO_indices)
        index = plot_PO_indices(PO_i);
        
        [T_PO, X_PO] = ode113(@Int_CR3BnSTM, [0, PO_data(index,c_Tp)], [PO_data(index,c_x0:c_zd0)'; stm0_colVec], options, prms);
        
        %%% Store data
        traj_POi{PO_i}.X = X_PO;
        traj_POi{PO_i}.T = T_PO;
    end
    
    %%% Loop and plot
    figure; hold all
    for PO_i = 1:length(traj_POi)
%         plot3(traj_POi{PO_i}.X(:,1), traj_POi{PO_i}.X(:,2), traj_POi{PO_i}.X(:,3), 'linewidth', lw, 'color', PO_color) % 0.1961    0.3922    0.7843
%         plot3(traj_POi{PO_i}.X(:,1), traj_POi{PO_i}.X(:,2), traj_POi{PO_i}.X(:,3), 'linewidth', lw, 'color', colors.sch.eighties.yellow) % 0.1961    0.3922    0.7843
        plot3(traj_POi{PO_i}.X(:,1), traj_POi{PO_i}.X(:,2), traj_POi{PO_i}.X(:,3), 'linewidth', lw, 'color', color_spectrum(PO_i,:)) % 0.1961    0.3922    0.7843
    end
    ax = gca;
    ax.FontSize = 14;
    PlotBoi3_CR3Bn(28)
    axis equal
%     plot3(rLPs_n(2,1),[0],[0],'^','markeredgecolor',colors.black, 'markerfacecolor', colors.blue)
    
    plot3(rLPs_n(1:2,1),[0,0],[0,0],'^','markeredgecolor',colors.black, 'markerfacecolor', colors.blue)
%     set(gca, 'Color', colors.orccaPPT)
%     set(gcf, 'Color', colors.orccaPPT)
%     plotSecondary(secondary)

    if sum(traj_POi{PO_i}.X(:,3) + traj_POi{PO_i}.X(:,6)) == 0
        view(0,90) % x-z view
    else
        view(0,0) % x-z view
    end


end

if 1+1==1
    
    
    newIndex = 100;
    
    [T_PO_new, X_PO_new] = ode113(@Int_CR3BnSTM, [0, PO_data(newIndex,c_Tp)], [PO_data(newIndex,c_x0:c_zd0)'; stm0_colVec], options, prms);
    
    stm_tf_t0                           = reshape(X_PO_new(end,7:42),6,6);
    monodromy                           = stm_tf_t0;
    [eigenVectors_new, eigenValues_new] = eig(monodromy);
    [S1, S2]                            = getStabilityIndices(diag(eigenValues_new));
    S1
    S2
    diag(eigenValues_new)
    

    figure;
    hold all
    PlotBoi3_CR3Bn(26)
    axis equal
%     plot3(rLPs_n(1:2,1),[0,0],[0,0],'^','markeredgecolor',colors.black, 'markerfacecolor', colors.blue)
    plotSecondary(secondary)
    
    plot3(X_PO_new(:,1), X_PO_new(:,2), X_PO_new(:,3), 'linewidth', lw, 'color', colors.sch.eighties.pink) % 0.1961    0.3922    0.7843


    latLons_PO = zeros(length(T_PO_new),2);
    for kk = 1:length(latLons_PO)
        [lat_deg, lon_deg] = BCR2latlon(X_PO_new(kk,1:3)', 'secondary', prms.u);
        latLons_PO(kk,:) = [lat_deg, lon_deg];
    end

    figure; hold all
    xlim([0 360])
    ylim([-90 90])
    PlotBoi2('Longitude, $^\circ$', 'Latitude, $^\circ$', 26, 'LaTex')
    h = image(xlim, -ylim, secondary.img);
    [lons_new] = convert_lon180_to_lon360(latLons_PO(:,2));
    plot(lons_new,latLons_PO(:,1),'.', 'color', colors.blue2)


end


%%% For saving trajectory
if 1+1 == 1
%     trajectory_X0 = PO_data(n_POs,c_x0:c_Tp)';
%     [trajectory_t, trajectory_X] = ode113(@Int_CR3BnSTM, [0, trajectory_X0(7)], [trajectory_X0(1:6); stm0_colVec], options, prms);
%     
%     newFileName = 'TrajFile';
%     savePath = '/Users/lukebury/Downloads/';
%     
%     % -------------------------------------------------
%     %%% Preparing Save File
%     % -------------------------------------------------
%     %%% Create unique filename
%     [uniqueFilename] = get_uniqueFilename(newFileName, savePath, 'txt');
%     
%     %%% Open File
%     datafile = fopen(uniqueFilename,'wt');
%     
%     % -------------------------------------------------
%     %%% Writing data
%     % -------------------------------------------------
%     %%% Write header
%     headerString = 'x,y,z,xd,yd,zd,t\n';
%     fprintf(datafile,headerString);
%     
%     %%% Write data
%     for kk = 1:length(trajectory_t)
%         
%         fprintf(datafile,'%1.16f,%1.16f,%1.16f,%1.16f,%1.16f,%1.16f,%1.16f\n',...
%             trajectory_X(kk,1), trajectory_X(kk,2), trajectory_X(kk,3),...
%             trajectory_X(kk,4), trajectory_X(kk,5), trajectory_X(kk,6),...
%             trajectory_t(kk));
%         
%     end
%     
%     %%% Close file
%     fclose(datafile);




end

%%% Plot PO with Shadows
if 1+1==1
    index = 1480;
    
    PO_0 = PO_data(index, c_x0:c_Tp)';
    [T_PO, X_PO] = ode113(@Int_CR3BnSTM, [0, PO_0(7)], [PO_0(1:6); stm0_colVec], options, prms);
    
    figure; hold all
    plot3(X_PO(:,1),X_PO(:,2),X_PO(:,3),'linewidth',2,'color',colors.mag)
%     plotTrajShadows(X_PO, 2, colors.grey, 'x', 0.965, 'y', 0.025, 'z', -4.5e-3)
    plotTrajShadows(X_PO, 2, colors.grey, 'x', 0.965, 'y', 0.035, 'z', -2.7e-2, 'bodyshadow', [1-prms.u, prms.R2])
    axis equal
    PlotBoi3_CR3Bn(26)
    plotSecondary(secondary)
    view(48,16)

end


%%% Propagate a perturbation
if 1+1==1
%     index = 524;
    index = 1000;
    
    pert = [10/rNorm; 0; 0; 0; 0; 0; 10*pi].*0;
    
    PO_0 = PO_data(index, c_x0:c_Tp)';
    PO_0 = PO_0 + pert;
    [T_PO, X_PO] = ode113(@Int_CR3BnSTM, [0, PO_0(7)], [PO_0(1:6); stm0_colVec], options, prms);
    
    figure; hold all
    plot3(X_PO(:,1),X_PO(:,2),X_PO(:,3),'linewidth',1.5,'color',colors.mag)
%     plotTrajShadows(X_PO, 2, colors.grey, 'x', 0.965, 'y', 0.025, 'z', -4.5e-3)
%     plotTrajShadows(X_PO, 2, colors.grey, 'x', 0.965, 'y', 0.035, 'z', -2.7e-2, 'bodyshadow', [1-prms.u, prms.R2])
    axis equal
    PlotBoi3_CR3Bn(26)
    plotSecondary(secondary)
    view(0,90)


end


%%% Progressively Draw
if 1+1==1
% %     plot_PO_indices = [270:285];
%     plot_PO_indices = getIndices_spacedByTpJcArclength(PO_data(:,c_Tp), PO_data(:,c_JC), 15);
%     color_spectrum = colorScale([colors.blue2; colors.mag], length(plot_PO_indices));
% 
%     
% %   %%% Loop and plot
%     figure; hold all
%     PlotBoi3_CR3Bn(26)
%     axis equal
%     plot3(rLPs_n(1:2,1),[0,0],[0,0],'^','markeredgecolor',colors.black, 'markerfacecolor', colors.blue)
% %     view(0, 90) % x-y
% %     view(90,0) % y-z
%     view(0,0) % x-z
%     
%     for PO_i = 1:length(plot_PO_indices)
%         index = plot_PO_indices(PO_i);
%         
%         [T_PO, X_PO] = ode113(@Int_CR3BnSTM, [0, PO_data(index,c_Tp)], [PO_data(index,c_x0:c_zd0)'; stm0_colVec], options, prms);
%         
%         plot3(X_PO(:,1), X_PO(:,2), X_PO(:,3), 'linewidth', lw, 'color', color_spectrum(PO_i,:)) % 0.1961    0.3922    0.7843
%         drawnow
%         pause(0.2)
%     end



end


%%% These were used in defense presentation
if 1+1==1
    %%% Chaotic 3B example
    X0 = [1.0176794902848845; 0; 0; 0; 0.0053; 0];
    [T_PO_new, X_PO_new] = ode113(@Int_CR3BnSTM, [0:.0001:3.2*pi], [X0; stm0_colVec], options, prms);
    figure;
    hold all
    PlotBoi3_CR3Bn(26)
    axis equal
    plot3(X_PO_new(:,1), X_PO_new(:,2), X_PO_new(:,3), '.', 'color', colors.cyan) % 0.1961    0.3922    0.7843
    set(gca, 'Color', colors.orccaPPT)
    set(gcf, 'Color', colors.orccaPPT)
    
    
    
    
    
    
    %%% Apollo orbit
    newIndex = 100;
    [T_PO_new, X_PO_new] = ode113(@Int_CR3BnSTM, [0, PO_data(newIndex,c_Tp)], [PO_data(newIndex,c_x0:c_zd0)'; stm0_colVec], options, prms);
    figure;
    hold all
    PlotBoi3_CR3Bn(26)
    axis equal
    plotSecondary(secondary)
    plot3(X_PO_new(:,1), X_PO_new(:,2), X_PO_new(:,3), 'linewidth', lw, 'color', colors.sch.eighties.pink) % 0.1961    0.3922    0.7843
    plotPrimary(primary, secondary)
    view(0,90)
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






