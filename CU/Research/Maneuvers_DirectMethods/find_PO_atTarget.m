% ========================================================================
%%% Description
% ========================================================================
% From a pre-computed database containing a family of periodic orbits,
% specify a set of desired Jacobi constants or time-periods at which to
% sample the family and create a new subset of POs at these Jacobi 
% constants

% Created: 10/02/19
% Author : Luke Bury, luke.bury@colorado.edu
% ========================================================================
%%% Initialization
% ========================================================================
clear
clc
close all
mbinPath = '~/CU_Google_Drive/Documents/MatGit/mbin';
POfamilyPath = '~/CU_Google_Drive/Documents/MatGit/mbin/Data/InitialConditions/PO_Families/';
savePath = '~/CU_Google_Drive/Documents/MatGit/mbin/Data/InitialConditions/PO_Families_Pretty/';
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
print_solution_index = 1;

spaceByJC = 0;
spaceByTp = 1;

savePrettyFamily = 0;

plot_surface_reconstruction = 0;

plot_secondary = 0;

if (spaceByJC + spaceByTp) == 2
    warning('These are mututally exclusive')
    return
end



if savePrettyFamily == 1
    warning('Save is on')
end

% ========================================================================
%%% Setup
% ========================================================================
% -------------------------------------------------
%%% Column identifiers for family data file
% -------------------------------------------------
c_x0_n                = 1;
c_y0_n                = 2;
c_z0_n                = 3;
c_xd0_n               = 4;
c_yd0_n               = 5;
c_zd0_n               = 6;
c_Tp_n                = 7;
c_JC                  = 8;
c_L2FlythroughVel_mps = 9;
c_landingVelocity_mps = 10;
c_stabilityIndex1     = 11;
c_stabilityIndex2     = 12;

if spaceByJC == 1
    c_target = c_JC;
elseif spaceByTp == 1
    c_target = c_Tp_n;
end

% -------------------------------------------------
%%% Choose PO family
% -------------------------------------------------
% for jj = 1:2
%     if jj == 1
%         famName = 'Saturn_Enceladus.CR3BP.L2_Vertical';
%     elseif jj == 2
%         famName = 'Saturn_Enceladus.CR3BP_J2pJ4pJ6pJ2s.L2_Vertical';
%     end

%%% Earth-Moon
% famName = 'Earth_Moon.CR3BP.L1_Lyapunov';
% famName = 'Earth_Moon.CR3BP.L1_Vertical';
% famName = 'Earth_Moon.CR3BP.L1_SHalo';

%%% Jupiter-Europa
% famName = 'Jupiter_Europa.CR3BP.L2_Lyapunov';
% famName = 'Jupiter_Europa.CR3BP.L2_Vertical';
% famName = 'Jupiter_Europa.CR3BP.L2_NHalo';
% famName = 'Jupiter_Europa.CR3BP.L2_SHalo';
% famName = 'Jupiter_Europa.CR3BP.L2_EasternAxial';
% famName = 'Jupiter_Europa.CR3BP.L2_WesternAxial';

% famName = 'Jupiter_Europa.CR3BP_J2pJ4pJ6pJ2s.L2_Lyapunov';
% famName = 'Jupiter_Europa.CR3BP_J2pJ4pJ6pJ2s.L2_Vertical';
% famName = 'Jupiter_Europa.CR3BP_J2pJ4pJ6pJ2s.L2_SHalo';


%%% Jupiter-Ganymede
% famName = 'Jupiter_Ganymede.CR3BP.L2_Lyapunov';
% famName = 'Jupiter_Ganymede.CR3BP.L2_Vertical';
% famName = 'Jupiter_Ganymede.CR3BP.L2_SHalo';


%%% Saturn-Enceladus
famName = 'Saturn_Enceladus.CR3BP.L2_Lyapunov';
% famName = 'Saturn_Enceladus.CR3BP.L2_Vertical';
% famName = 'Saturn_Enceladus.CR3BP.L2_SHalo';
% famName = 'Saturn_Enceladus.CR3BP.SLeadingSneakerToe';
% famName = 'Saturn_Enceladus.CR3BP_J2pJ4pJ6pJ2s.L2_Lyapunov';
% famName = 'Saturn_Enceladus.CR3BP_J2pJ4pJ6pJ2s.L2_Vertical';
% famName = 'Saturn_Enceladus.CR3BP_J2pJ4pJ6pJ2s.L2_SHalo';

%%% Uranus-Cordelia
% famName = 'Uranus_Cordelia.CR3BP.L2_Lyapunov';

%%% Uranus-Ophelia
% famName = 'Uranus_Ophelia.CR3BP.L1_Lyapunov';

%%% Neptune-Triton
% famName = 'Neptune_Triton.CR3BP.L2_Lyapunov';
% famName = 'Neptune_Triton.CR3BP.L2_Vertical';
% famName = 'Neptune_Triton.CR3BP.L2_SHalo';

% --------------------------
% Set primary & secondary
% --------------------------
[primary, secondary] = assignPrimaryAndSecondary_CR3BP(famName, bodies);

% --------------------------
%%% Load data and find JC range
% --------------------------
%%% Load PO family
FamilyDataFile = [POfamilyPath, famName, '.txt'];
POFamilyData = dlmread(FamilyDataFile,',',1,0);

%%% Find JC range
minJC = min(POFamilyData(:,c_JC));
maxJC = max(POFamilyData(:,c_JC));

%%% Find Tp range
minTp = min(POFamilyData(:,c_Tp_n));
maxTp = max(POFamilyData(:,c_Tp_n));

% -------------------------------------------------
%%% Set vector of target values to find POs at
% -------------------------------------------------
n_Solutions = 500;

if spaceByJC     == 1
    target_vec = linspace(minJC, maxJC, n_Solutions);
    % JC_des_vec = linspace(3.002200842600000, 3.003609009400000, n_Solutions); % For comparing vanilla to ZH Europa L2 Lyapunovs
    % JC_des_vec = linspace(3.000168585100000, 3.003609190100000, n_Solutions); % For comparing vanilla to ZH Europa L2 Verticals
    % JC_des_vec = linspace(3.000987456200000, 3.003326025500000, n_Solutions); % For comparing vanilla to ZH Europa L2 SHalos
    % JC_des_vec = linspace(, , n_Solutions); % For comparing vanilla to ZH Enceladus L2 Lyapunovs
    % JC_des_vec = linspace(, , n_Solutions); % For comparing vanilla to ZH Enceladus L2 Verticals
    % JC_des_vec = linspace(, , n_Solutions); % For comparing vanilla to ZH Enceladus L2 SHalos
    % JC_des_vec = [3.001];
    % JC_des_vec = [minJC];
    % % Cordelia particle
    % 3.000473172601698
    % % Ophelia particle
    % 3.001926680861342
elseif spaceByTp == 1
    target_vec = linspace(minTp, maxTp, n_Solutions);

%     target_vec = linspace(3.0815, 3.4747, n_Solutions); % For comparing vanilla to ZH Europa L2 Lyapunovs
%     target_vec = linspace(3.2, 6.0413, n_Solutions); % For comparing vanilla to ZH Europa L2 Verticals
%     target_vec = linspace(2.1952, 3.1242, n_Solutions); % For comparing vanilla to ZH Europa L2 SHalos
%     target_vec = linspace(2.1952, 3.128, n_Solutions); % For comparing vanilla to ZH Europa L2 SHalos
    
%     target_vec = linspace(3.185, 3.8947, n_Solutions); % For comparing vanilla to ZH Enceladus L2 Lyapunovs
%     target_vec = linspace(3.32, 6.2537, n_Solutions); % For comparing vanilla to ZH Enceladus L2 Verticals
%     target_vec = linspace(2.0559, 3.091, n_Solutions); % For comparing vanilla to ZH Enceladus L2 SHalos
end

% target_vec = fliplr(target_vec);

n_Solutions = 20;
% target_vec = linspace(target_vec(1), target_vec(3), n_Solutions);
% target_vec = linspace(3.1675396043512030, 3.1675271565021674, n_Solutions);

target_vec = linspace(3.0424294837044816, 3.0437044881280650, n_Solutions);
% n_Solutions = 20;
% target_vec = linspace(target_vec(1), target_vec(8), n_Solutions);


% ------------------------------------------------- 
%%% Other options
% -------------------------------------------------


%%% Error tolerance for constraint vector in multiple shooter
error_tol = 1e-13; 

%%% Number of nodes for multiple shooter
n_Nodes = 13; 

%%% Artificial scalar for increasing or decreasing the step size for 
%%% orbit correction. 1 has no effect, 2 and 3 are sometimes useful,
%%% sometimes not
stepSize = 1;

%%% Maximum number of multiple-shooter iterations for any given family
%%% member before kicking out of the loop and adjusting the tuning
%%% parameters
iterMax = 500;

% -------------------------------------------------
%%% System
% -------------------------------------------------
%%% Normalizing constants
rNorm = secondary.a;         % n <-> km
tNorm = 1/secondary.meanMot; % n <-> sec
vNorm = rNorm / tNorm;       % n <-> km/sec

%%% Setting parameter structure
prms.u  = secondary.MR;
prms.R1 = primary.R / rNorm;
prms.R2 = secondary.R_n;

if contains(famName,'.CR3BP_J2pJ4pJ6pJ2s.')
    prms.J2p = primary.J2; prms.J4p = primary.J4; prms.J6p = primary.J6; prms.J2s = secondary.J2;
end

%%% Equillibrium Points
if contains(famName,'.CR3BP.')
    rLPs_n = EquilibriumPoints(prms.u);
elseif contains(famName,'.CR3BP_J2pJ4pJ6pJ2s.')
    rLPs_n = collinearEquilibriumPoints_ZH(prms);
end

% ------------------------------------------------- 
%%% Integration Options
% -------------------------------------------------
tol     = 1e-13;
options = odeset('RelTol',tol,'AbsTol',tol);
options_XZStop = odeset('Event',@event_yEqualsZeroPastL2,'RelTol',tol,'AbsTol',tol);

% ------------------------------------------------- 
%%% Preparing to find POs
% -------------------------------------------------
%%% Preallocating space for outputs
POs                     = NaN(n_Solutions,7);
stabilityIndices        = NaN(n_Solutions,2);
jacobiConstants         = NaN(n_Solutions,1);
L2ExcessVelocities_mps  = NaN(n_Solutions,1);
landingVelocities_mps   = NaN(n_Solutions,1);

%%% Initialize STM
stm0        = eye(6);
stm0_colVec = reshape(stm0,36,1);

% -------------------------------------------------
%%% Loop through desired JC values and find PO at each
% -------------------------------------------------
solutionIndex = 0;
for target_i = target_vec
    if spaceByJC     == 1
        if target_i < minJC
            warning('Requested JC value is below minimum value from family')
            return
        elseif target_i > maxJC
            warning('Requested JC value is above maximum value from family')
            return
        end
    end
    
    %%% Count current PO
    solutionIndex = solutionIndex + 1;
    
    %%% Find index of closest solution in target-space
    indexOfClosest = find(abs(POFamilyData(:,c_target) - target_i) == min(abs(POFamilyData(:,c_target) - target_i)));
    
    if isequal(size(indexOfClosest),[1,1]) == 0
        indexOfClosest = indexOfClosest(1);
    end
    
    X0_guess_n = (POFamilyData(indexOfClosest, c_x0_n:c_zd0_n))';
    T_guess_n = POFamilyData(indexOfClosest, c_Tp_n);
    
    % ------------------------------------------------- 
    %%% Integrate guess and divide into nodes for multiple shooter
    % -------------------------------------------------
    %%% Integrate again to discretize into nodes evenly spaced in time
    if contains(famName,'.CR3BP.')
        [T_nodes, X_nodes] = ode113(@Int_CR3Bn, linspace(0, T_guess_n,n_Nodes+1), X0_guess_n, options, prms);
    elseif contains(famName,'.CR3BP_J2pJ4pJ6pJ2s.')
        [T_nodes, X_nodes] = ode113(@Int_CR3Bn_J2pJ4pJ6pJ2s, linspace(0, T_guess_n,n_Nodes+1), X0_guess_n, options, prms);
    end

    %%% Setting initial nodes (getting rid of repeat end state)
    X_nodes = X_nodes(1:n_Nodes,1:6);

    %%% Create first free-variable column vector
    F_new   = X_nodes';
    F_new   = F_new(:);
    if spaceByJC == 1
        F_new   = [F_new; T_guess_n];
    end

    %%% Set dynamics function
    if contains(famName,'.CR3BP.')
        dynamicsFunc = @Int_CR3BnSTM;
    elseif contains(famName,'.CR3BP_J2pJ4pJ6pJ2s.')
        dynamicsFunc = @Int_CR3BnSTM_J2pJ4pJ6pJ2s;
    end

    
    
    % -------------------------------------------------
    % Run  multiple shooting to find the correct PO
    % -------------------------------------------------
    % --------------------------
    %%% Pre-While-Loop work
    % --------------------------
    %%% Set while-loop/convergence variables
    iteration_counter = 0;
    constraint_error  = 100;

    % --------------------------
    % Enter while loop - attempt to converge on next PO in family
    % --------------------------
    while (constraint_error > error_tol) && (iteration_counter < iterMax)

        %%% Count iteration
        iteration_counter = iteration_counter + 1;

        % --------------------------
        %%% Using multiple shooting - loop through nodes, integrate, and 
        %%% populate constraint vector and create DF matrix
        % --------------------------
        if spaceByJC     == 1
            [DF_mat, constraints] = multShooter_stateContinuity_JCFixed(n_Nodes, F_new, dynamicsFunc, options, prms, target_i);
        elseif spaceByTp == 1
            [DF_mat, constraints] = multShooter_stateContinuity_TpFixed(n_Nodes, F_new, dynamicsFunc, options, prms, target_i);
        end

        % --------------------------
        %%% Determine current error and create next F if necessary
        % --------------------------
        %%% Compute error
        constraint_error = norm(constraints);

        %%% Compute new free-variable vector if error is not converged
        if (constraint_error > error_tol)
            warning('off','MATLAB:nearlySingularMatrix')
            F_new = F_new - stepSize * DF_mat'*((DF_mat*(DF_mat'))\constraints);
            warning('on','MATLAB:nearlySingularMatrix')
        end

    end % while constraint_error > error_tol

    if constraint_error > error_tol
        warning('Solution not converged on')
    end
    
    
    
    
    
    
    
    
    
    
    if print_solution_index == 1
        fprintf('Solution %1d ... error %1.1e ... iterations: %1d\n', solutionIndex, constraint_error, iteration_counter)
    end % print_PO_index
    
    % ------------------------------------------------- 
    %%% Store solution
    % -------------------------------------------------
    if spaceByJC     == 1
        POs(solutionIndex, :) = [F_new(1:6)', F_new(end)];
    elseif spaceByTp == 1
        POs(solutionIndex, :) = [F_new(1:6)', target_i];
    end
    
    % ------------------------------------------------- 
    %%% Find stability indices
    % -------------------------------------------------
    %%% Integrate new PO
    if contains(famName,'.CR3BP.')
        [~, X_new] = ode113(@Int_CR3BnSTM, [0, POs(solutionIndex, 7)], [POs(solutionIndex, 1:6)'; stm0_colVec], options, prms);
    elseif contains(famName,'.CR3BP_J2pJ4pJ6pJ2s.')
        [~, X_new] = ode113(@Int_CR3BnSTM_J2pJ4pJ6pJ2s, [0, POs(solutionIndex, 7)], [POs(solutionIndex, 1:6)'; stm0_colVec], options, prms);
    end

    %%% Stability indices of new PO
    stm_tf_t0                           = reshape(X_new(end,7:42),6,6);
    monodromy                           = stm_tf_t0;
    [eigenVectors_new, eigenValues_new] = eig(monodromy);
    [S1, S2]                            = getStabilityIndices(diag(eigenValues_new));
    stabilityIndices(solutionIndex,:)          = [S1, S2];
    
    % ------------------------------------------------- 
    %%% Store other important data
    % -------------------------------------------------
    jacobiConstants(solutionIndex)        = getJacobiConstant_ZH(F_new(1:6)', prms);
    L2ExcessVelocities_mps(solutionIndex) = JC_2_L2FlyoverVelocity(jacobiConstants(solutionIndex), prms, rLPs_n(2,:), vNorm);
    landingVelocities_mps(solutionIndex)  = JC_2_LandingVelocity(jacobiConstants(solutionIndex), prms, vNorm);
        
end

% ========================================================================
%%% Plotting Results
% ========================================================================
% -------------------------------------------------
%%% Plot new POs
% -------------------------------------------------
figure(100); hold all
axis equal
PlotBoi3('$X_n$','$Y_n$','$Z_n$',20,'LaTex')

if plot_secondary
    if isfield(secondary,'img') == 1
        plotSecondary(secondary, 0.8)
    elseif isfield(secondary,'color') == 1
        plotBody3(secondary.R_n,[1-secondary.MR,0,0],secondary.color,0.8)
    else
        plotBody3(secondary.R_n,[1-secondary.MR,0,0],colors.std.grn,0.8)
    end
end

allTrajs_CR3BP    = [];
allTrajs_CR3BP_ZH = [];

for kk = 1:n_Solutions
    
    %%% Integrate
    if contains(famName,'.CR3BP.')
        [~, X_PO] = ode113(@Int_CR3Bn, [0, POs(kk,7)], POs(kk,1:6)', options, prms);
        allTrajs_CR3BP = [allTrajs_CR3BP; X_PO];
        
        %%% Plot
        p_CR3BP = plot3(X_PO(:,1),X_PO(:,2),X_PO(:,3),'r','linewidth',1.5);
    elseif contains(famName,'.CR3BP_J2pJ4pJ6pJ2s.')
        [~, X_PO] = ode113(@Int_CR3Bn_J2pJ4pJ6pJ2s, [0, POs(kk,7)], POs(kk,1:6)', options, prms);
        allTrajs_CR3BP_ZH = [allTrajs_CR3BP_ZH; X_PO];
        
        %%% Plot
        p_CR3BP_ZH = plot3(X_PO(:,1),X_PO(:,2),X_PO(:,3),'b','linewidth',1.5);
    end    
    
end

% -------------------------------------------------
%%% Plot Stability Indices
% -------------------------------------------------
%%% Find stability indices
S1 = stabilityIndices(:,1);
S2 = stabilityIndices(:,2);

figure('position',[209 322 948 302])
subplot(1,2,1); hold all
p1 = plot(jacobiConstants, abs(S1),'o','markeredgecolor',colors.std.blue,'markerfacecolor',colors.std.ltblue);
p2 = plot(jacobiConstants, abs(S2),'o','markeredgecolor',colors.std.red,'markerfacecolor',colors.std.ltred);
% plot(unique([min([jacobiConstants, jacobiConstants]) max([jacobiConstants, jacobiConstants])]),[2 2],'k','linewidth',1)
plot([min(jacobiConstants), max(jacobiConstants)],[2 2],'k','linewidth',1)
PlotBoi2('Jacobi Constant','Stability Indices',18,'LaTex')
legend([p1 p2],'S_1','S_2')
% xlim(unique([min([jacobiConstants, jacobiConstants]) max([jacobiConstants, jacobiConstants])]))
xlim([min(jacobiConstants), max(jacobiConstants)])

subplot(1,2,2); hold all
p1 = plot(jacobiConstants, abs(S1),'o','markeredgecolor',colors.std.blue,'markerfacecolor',colors.std.ltblue);
p2 = plot(jacobiConstants, abs(S2),'o','markeredgecolor',colors.std.red,'markerfacecolor',colors.std.ltred);
% plot(unique([min([jacobiConstants, jacobiConstants]) max([jacobiConstants, jacobiConstants])]),[2 2],'k','linewidth',1)
plot([min(jacobiConstants), max(jacobiConstants)],[2 2],'k','linewidth',1)
PlotBoi2('Jacobi Constant','',18,'LaTex')
ylim([1.9 2.1])
xlim([min(jacobiConstants), max(jacobiConstants)])


% -------------------------------------------------
%%% Plot Surface Reconstruction
% -------------------------------------------------
if plot_surface_reconstruction
    
    oldTraj = [];
    oldTraj_ZH = [];
    PatchedShape_x = [];
    PatchedShape_y = [];
    PatchedShape_z = [];
    PatchedShape_ZH_x = [];
    PatchedShape_ZH_y = [];
    PatchedShape_ZH_z = [];
    for kk = 1:n_Solutions

        if kk > 1
            if contains(famName,'.CR3BP.')
                oldTraj = X_PO;
            elseif contains(famName,'.CR3BP_J2pJ4pJ6pJ2s.')
                oldTraj_ZH = X_PO_ZH;
            end
        end

        %%% Integrate
        if contains(famName,'.CR3BP.')
            [~, X_PO] = ode113(@Int_CR3Bn, linspace(0, POs(kk,7), 100), POs(kk,1:6)', options, prms);

        elseif contains(famName,'.CR3BP_J2pJ4pJ6pJ2s.')
            [~, X_PO_ZH] = ode113(@Int_CR3Bn_J2pJ4pJ6pJ2s, linspace(0, POs(kk,7), 100), POs(kk,1:6)', options, prms);

        end    

        if kk > 1
            if contains(famName,'.CR3BP.')
                [x_patch, y_patch, z_patch] = patchLines(oldTraj, X_PO);
                
                PatchedShape_x = [PatchedShape_x, x_patch]; 
                PatchedShape_y = [PatchedShape_y, y_patch]; 
                PatchedShape_z = [PatchedShape_z, z_patch]; 
            elseif contains(famName,'.CR3BP_J2pJ4pJ6pJ2s.')
                [x_patch_ZH, y_patch_ZH, z_patch_ZH] = patchLines(oldTraj_ZH, X_PO_ZH);

                PatchedShape_ZH_x = [PatchedShape_ZH_x, x_patch_ZH]; 
                PatchedShape_ZH_y = [PatchedShape_ZH_y, y_patch_ZH]; 
                PatchedShape_ZH_z = [PatchedShape_ZH_z, z_patch_ZH];
            end 
            
            

            
        end

    end
    
    figure(500); hold all
    if contains(famName,'.CR3BP.')
        patch(PatchedShape_x, PatchedShape_y, PatchedShape_z,'r')
    elseif contains(famName,'.CR3BP_J2pJ4pJ6pJ2s.')
        patch(PatchedShape_ZH_x, PatchedShape_ZH_y, PatchedShape_ZH_z,'b')
    end
    PlotBoi3('$x_n$','$y_n$','$z_n$',20,'LaTex')
    axis equal
    view(3)
    
    if plot_secondary
        if isfield(secondary,'img') == 1
            plotSecondary(secondary, 0.8)
        elseif isfield(secondary,'color') == 1
            plotBody3(secondary.R_n,[1-secondary.MR,0,0],secondary.color,0.8)
        else
            plotBody3(secondary.R_n,[1-secondary.MR,0,0],colors.std.grn,0.8)
        end
    end
end % plot_surface_reconstruction





%% =======================================================================
%%% Clean up the data
% ========================================================================
% -------------------------------------------------
%%% Integrating solution to x-y plane crossing for state simplicity
% -------------------------------------------------
POs_cleanedUp = NaN(size(POs));

for kk = 1:size(POs,1)
    %%% Integrate
    if contains(famName,'.CR3BP.')
        [~, X_BCR_n_XY, time_event, X_event, index_event] = ...
        ode113(@Int_CR3Bn, [0, POs(kk,end)], POs(kk,1:6)', options_XZStop, prms);
    elseif contains(famName,'.CR3BP_J2pJ4pJ6pJ2s.')
        [~, X_BCR_n_XY, time_event, X_event, index_event] = ...
        ode113(@Int_CR3Bn_J2pJ4pJ6pJ2s, [0, POs(kk,end)], POs(kk,1:6)', options_XZStop, prms);
    end
    
    
    
    %%% Attempt, but if trajectory doesn't cross x-y plane, then forget it
    if isempty(X_event) == 0
        POs_cleanedUp(kk,:) = [X_event(end,1:6),POs(kk,end)];
    else
        POs_cleanedUp(kk,:) = [POs(kk,1:6),POs(kk,end)];
    end
end
    

% ========================================================================
%%% Saving Data
% ========================================================================
if savePrettyFamily == 1
    % -------------------------------------------------
    %%% Preparing Save File
    % -------------------------------------------------
    fileAlreadyExists = 1;
    fileVersion       = 0;
    famName_pretty    = [famName, sprintf('.pretty_%1d',length(target_vec))];
    famName_new       = famName_pretty;
    while fileAlreadyExists == 1
        if isfile([savePath,famName_new,'.txt']) == 1
            fileVersion = fileVersion + 1;
            famName_new = sprintf('%s_%1d',famName_pretty,fileVersion);
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
    dataStop = size(POs_cleanedUp,1);
    
    
    for kk = dataStart:dataStop
        fprintf(datafile,'%1.16f,%1.16f,%1.16f,%1.16f,%1.16f,%1.16f,%1.16f,%1.10f,%1.3f,%1.3f,%1.4f,%1.4f\n',...
            POs_cleanedUp(kk,1), POs_cleanedUp(kk,2), POs_cleanedUp(kk,3),...
            POs_cleanedUp(kk,4), POs_cleanedUp(kk,5), POs_cleanedUp(kk,6),...
            POs_cleanedUp(kk,end), jacobiConstants(kk), L2ExcessVelocities_mps(kk),...
            landingVelocities_mps(kk), stabilityIndices(kk,1), stabilityIndices(kk,2));
            
    end

    %%% Close file
    fclose(datafile);






end


% ========================================================================
%%% Closeout
% ========================================================================
tocWhole = toc(ticWhole);
fprintf('Elapsed time: %1.4f seconds\n',tocWhole)




