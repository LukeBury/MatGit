% ========================================================================
%%% Description
% ========================================================================
% For continuing familes of periodic orbits perturbed by zonal harmonics
% with just a multiple shooter rather than pseudo-arclength continuation

% Created: 06/5/20
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
savePOs = false;

% ========================================================================
%%% Setup
% ========================================================================
% -------------------------------------------------
%%% Choose PO family
% -------------------------------------------------
% famName = 'Jupiter_Europa';
% PO_0 = PO_ICs.Jupiter_Europa.CR3BP_J2pJ4pJ6pJ2s.L2_Lyapunov;

famName = 'Saturn_Enceladus';
% PO_0 = PO_ICs.Saturn_Enceladus.CR3BP_J2pJ4pJ6pJ2s.L2_Lyapunov;
989
PO_0 = [1.0006564111824836,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0240983700639387,0.0000000000000000,5.5225980168712114];


%%%% Set primary & secondary
[primary, secondary] = assignPrimaryAndSecondary_CR3BP(famName, bodies);

% ------------------------------------------------- 
%%% Shooter options
% -------------------------------------------------
%%% Number of nodes for multiple shooter
N_nodes = 12; 

%%% Error tolerance for constraint vector in multiple shooter
error_tol = 1e-13; 

%%% Maximum number of multiple-shooter iterations 
maxIter = 500;

% -------------------------------------------------
%%% System
% -------------------------------------------------
%%% Factor for normalizing distances
rNorm = secondary.a; % n <-> km

%%% Setting parameters structure
prms.u = secondary.MR;
prms.R1 = primary.R / rNorm;
prms.R2 = secondary.R_n;
prms.J2p = primary.J2;
prms.J4p = primary.J4;
prms.J6p = primary.J6;
prms.J2s = secondary.J2;

%%% Getting normalized mean motion
tNorm = sqrt((rNorm^3)/(bodies.constants.G*(primary.mass + secondary.mass)));
vNorm = rNorm/tNorm;
prms.n = secondary.meanMot * tNorm;
% prms.n = 1 + 3*(prms.J2p*prms.R1*prms.R1 + prms.J2s*prms.R2*prms.R2)/2;

prms.J2p = primary.J2; prms.J4p = primary.J4; prms.J6p = primary.J6; prms.J2s = secondary.J2;

%%% Equillibrium Points
rLPs_n = collinearEquilibriumPoints_ZH(prms);


% ------------------------------------------------- 
%%% Integration Options
% -------------------------------------------------
tol     = 1e-13;
options = odeset('RelTol',tol,'AbsTol',tol);


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

% -------------------------------------------------
%%% Create spectrum of x0 values
% -------------------------------------------------
n_Solutions = 300;
x0s = linspace(PO_0(1), 1-prms.u+1/rNorm,n_Solutions);

% ========================================================================
%%% Continue family
% ========================================================================
% -------------------------------------------------
%%% Preallocate
% -------------------------------------------------
POs                        = NaN(n_Solutions,7);
jacobiConstants            = NaN(n_Solutions,1);
L2FlythroughVelocities_mps = NaN(n_Solutions,1);
landingVelocities_mps      = NaN(n_Solutions,1);
stabilityIndices           = NaN(n_Solutions,2);

% % -------------------------------------------------
% %%% Fill in the first slot
% % -------------------------------------------------
% JC = getJacobiConstant_ZH(PO_0(1:6),prms);
% L2FTV_mps = JC_2_L2FlyoverVelocity(JC, prms, rLPs_n(2,:), vNorm);
% landingV_mps = JC_2_approxLandingVelocity(JC,prms,vNorm);

% -------------------------------------------------
%%% Loop through x0 values and continue the family
% -------------------------------------------------
solutionIndex = 0;
t_guess = PO_0(end);
X0_guess = PO_0(1:6);

figure(555); hold all
PlotBoi3_CR3Bn(20)
for x0 = x0s
    %%% Update solution index
    solutionIndex = solutionIndex + 1;
    
    %%% Get nodes
    [~, X_nodes] = get_nodes(X0_guess', [0, t_guess], N_nodes+1, @Int_CR3Bn_ZH, options, prms);
    
    % -------------------------------------------------
    %%% Create the free-variable vector
    % -------------------------------------------------
    F_vec = zeros(N_nodes*6+1,1);
    for kk = 1:N_nodes
        F_vec((kk*6-5):(kk*6)) = X_nodes(kk,:)';
    end
    F_vec(end) = t_guess;
    F_vec = F_vec([5,7:end]);
    
    
    
    
    
    
    % -------------------------------------------------
    %%% Initialize constrain error and loop over the linear update equation
    %%% with a multiple shooter
    % -------------------------------------------------
    c_norm = 1000;
    iter = 0;

    while (c_norm > error_tol) && (iter < maxIter)
        iter = iter + 1;

        [DF, constraints_vec] = multShooter_LyapPO_stateContinuity(N_nodes, x0, F_vec, @Int_CR3BnSTM_J2pJ4pJ6pJ2s, options, prms);

        c_norm = norm(constraints_vec);

        if (c_norm > error_tol)
            warning('off','MATLAB:nearlySingularMatrix')
            F_vec = F_vec - DF'*((DF*(DF'))\constraints_vec);
            warning('on','MATLAB:nearlySingularMatrix')
        end
    end

    if c_norm > error_tol
        warning('Solution not converged on')
        fprintf('solutionIndex = %d\n',solutionIndex);
        break
    else
        POs(solutionIndex,:) = [x0; 0; 0; 0; F_vec(1); 0; F_vec(end)];
        
        
        X0_guess = POs(solutionIndex,1:6)';
        t_guess  = POs(solutionIndex,7);
    
        stm0_colVec = reshape(eye(6),36,1);
        [~, XSTM_PO] = ode113(@Int_CR3BnSTM_J2pJ4pJ6pJ2s, [0, POs(solutionIndex,7)], [POs(solutionIndex,1:6)'; stm0_colVec], options, prms);

        plot3(XSTM_PO(:,1),XSTM_PO(:,2),XSTM_PO(:,3),'r')
        drawnow
        
        
        stm_tf_t0                           = reshape(XSTM_PO(end,7:42),6,6);
        monodromy                           = stm_tf_t0;
        [eigenVectors_new, eigenValues_new] = eig(monodromy);
        [S1, S2]                            = getStabilityIndices(diag(eigenValues_new));
        
        jacobiConstants(solutionIndex)              = getJacobiConstant_ZH(POs(solutionIndex, 1:6), prms);
        L2FlythroughVelocities_mps(solutionIndex)   = JC_2_L2FlyoverVelocity(jacobiConstants(solutionIndex), prms, rLPs_n(2,:), vNorm);
        landingVelocities_mps(solutionIndex)        = JC_2_approxLandingVelocity(jacobiConstants(solutionIndex), prms, vNorm);
        stabilityIndices(solutionIndex,:)           = [S1, S2];
        
        fprintf('Solution %1d ... error %1.1e ... iterations: %1d\n', solutionIndex, c_norm, iter)
        fprintf('Elapsed time: %1.4f seconds\n',toc(ticWhole))
    end

end

% -------------------------------------------------
%%% Plot Stability Indices
% -------------------------------------------------
%%% Find stability indices
S1 = stabilityIndices(:,1);
S2 = stabilityIndices(:,2);

figure
subplot(1,2,1); hold all
p1 = plot(jacobiConstants, abs(S1),'o','markeredgecolor',colors.blue,'markerfacecolor',colors.ltblue);
p2 = plot(jacobiConstants, abs(S2),'o','markeredgecolor',colors.red,'markerfacecolor',colors.ltred);
plot([min(jacobiConstants), max(jacobiConstants)],[2 2],'k','linewidth',1)
PlotBoi2('Jacobi Constant','Stability Indices',18,'LaTex')
legend([p1 p2],'S_1','S_2')
% xlim(unique([min([jacobiConstants, jacobiConstants]) max([jacobiConstants, jacobiConstants])]))
xlim([min(jacobiConstants), max(jacobiConstants)])

subplot(1,2,2); hold all
p1 = plot(jacobiConstants, abs(S1),'o','markeredgecolor',colors.blue,'markerfacecolor',colors.ltblue);
p2 = plot(jacobiConstants, abs(S2),'o','markeredgecolor',colors.red,'markerfacecolor',colors.ltred);
% plot(unique([min([jacobiConstants, jacobiConstants]) max([jacobiConstants, jacobiConstants])]),[2 2],'k','linewidth',1)
plot([min(jacobiConstants), max(jacobiConstants)],[2 2],'k','linewidth',1)
PlotBoi2('Jacobi Constant','',18,'LaTex')
ylim([1.9 2.1])
xlim([min(jacobiConstants), max(jacobiConstants)])









% ========================================================================
%%% Saving Data
% ========================================================================
if savePOs 
    %%% Cleaning data
    POs_new = POs(~isnan(POs(:,1)),:);
    jacobiConstants_new = jacobiConstants(~isnan(jacobiConstants(:,1)),:);
    L2FlythroughVelocities_mps_new = L2FlythroughVelocities_mps(~isnan(L2FlythroughVelocities_mps(:,1)),:);
    landingVelocities_mps_new = landingVelocities_mps(~isnan(landingVelocities_mps(:,1)),:);
    stabilityIndices_new = stabilityIndices(~isnan(stabilityIndices(:,1)),:);

    % -------------------------------------------------
    %%% Preparing Save File
    % -------------------------------------------------
    fileAlreadyExists = 1;
    fileVersion       = 0;
    famName_new       = [famName,'.CR3BP_J2pJ4pJ6pJ2s.L2_Lyapunov'];
    while fileAlreadyExists == 1
        if isfile([savePath,famName_new,'.txt']) == 1
            fileVersion = fileVersion + 1;
            famName_new = sprintf('%s_%1d',famName_new,fileVersion);
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
    dataStop = size(POs_new,1);
    
    
    for kk = dataStart:dataStop
        fprintf(datafile,'%1.16f,%1.16f,%1.16f,%1.16f,%1.16f,%1.16f,%1.16f,%1.10f,%1.3f,%1.3f,%1.4f,%1.4f\n',...
            POs_new(kk,1), POs_new(kk,2), POs_new(kk,3),...
            POs_new(kk,4), POs_new(kk,5), POs_new(kk,6),...
            POs_new(kk,end), jacobiConstants_new(kk), L2FlythroughVelocities_mps_new(kk),...
            landingVelocities_mps_new(kk), stabilityIndices_new(kk,1), stabilityIndices_new(kk,2));
            
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
















