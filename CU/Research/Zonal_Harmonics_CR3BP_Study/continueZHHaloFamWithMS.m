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
step_by_zd = false;
step_by_zd_withZFixed = false;
step_by_Tp_withZFixed = true;

savePOs = false;
saveOtherHalos = false;

% ========================================================================
%%% Setup
% ========================================================================
% -------------------------------------------------
%%% Choose PO family
% -------------------------------------------------
famName = 'Jupiter_Europa';
% PO_0 = PO_ICs.Jupiter_Europa.CR3BP_J2pJ4pJ6pJ2s.L2_SHalo;
PO_0 = [1.0172272846853616,-0.0109451486318729,0.0000000000000000,-0.0066862977324815,0.0046539654987015,0.0127067739236348,3.0994358012136050];

% famName = 'Saturn_Enceladus';
% PO_0 = PO_ICs.Saturn_Enceladus.CR3BP_J2pJ4pJ6pJ2s.L2_SHalo;
% PO_0 = [1.0000233857220133,-0.0007865570464379,0.0000000000000000,-0.0028844553732887,0.0142144390927371,0.0153398310104072,1.8000000000000000];


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
n_Solutions = 600;
if step_by_zd || step_by_zd_withZFixed
    % zd0s = linspace(PO_0(6), PO_0(6)*2, n_Solutions);
    zd0s = linspace(PO_0(6), PO_0(6)*1000, n_Solutions);
    
elseif step_by_Tp_withZFixed
    Tps = linspace(PO_0(7), 0, n_Solutions);

end
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
% solutionIndex = 0;
t_guess = PO_0(end);
X0_guess = PO_0(1:6);

figure(555); hold all
PlotBoi3_CR3Bn(20)
% for zd0 = zd0s
for solutionIndex = 1:n_Solutions
    %%% Update solution index
%     solutionIndex = solutionIndex + 1;
    
    %%% Get nodes
    [~, X_nodes] = get_nodes(X0_guess', [0, t_guess], N_nodes+1, @Int_CR3Bn_ZH, options, prms);
    
    % -------------------------------------------------
    %%% Create the free-variable vector
    % -------------------------------------------------
    if step_by_zd
        zd0 = zd0s(solutionIndex);
        
        F_vec = zeros(N_nodes*6+1,1);
        for kk = 1:N_nodes
            F_vec((kk*6-5):(kk*6)) = X_nodes(kk,:)';
        end
        F_vec(end) = t_guess;
%         F_vec = F_vec([1:2,4:5, 7:end]);
        F_vec = F_vec([1:5, 7:end]);
        
    elseif step_by_zd_withZFixed
        zd0 = zd0s(solutionIndex);
        
        F_vec = zeros(N_nodes*6+1,1);
        for kk = 1:N_nodes
            F_vec((kk*6-5):(kk*6)) = X_nodes(kk,:)';
        end
        F_vec(end) = t_guess;
        F_vec = F_vec([1:2,4:5, 7:end]);
        
    elseif step_by_Tp_withZFixed
        Tp = Tps(solutionIndex);
        
        F_vec = zeros(N_nodes*6,1);
        for kk = 1:N_nodes
            F_vec((kk*6-5):(kk*6)) = X_nodes(kk,:)';
        end
        F_vec = F_vec([1,2,4:end]);
    end
    
    
    
    
    
    
    % -------------------------------------------------
    %%% Initialize constrain error and loop over the linear update equation
    %%% with a multiple shooter
    % -------------------------------------------------
    c_norm = 1000;
    iter = 0;

    while (c_norm > error_tol) && (iter < maxIter)
        iter = iter + 1;
        
        if step_by_zd
            [DF, constraints_vec] = multShooter_stateContinuity_zdFixed(N_nodes, zd0, F_vec, @Int_CR3BnSTM_J2pJ4pJ6pJ2s, options, prms);
        elseif step_by_zd_withZFixed
            [DF, constraints_vec] = multShooter_stateContinuity_zzdFixed(N_nodes, 0, zd0, F_vec, @Int_CR3BnSTM_J2pJ4pJ6pJ2s, options, prms);
        elseif step_by_Tp_withZFixed
            [DF, constraints_vec] = multShooter_stateContinuity_zTpFixed(N_nodes, F_vec, @Int_CR3BnSTM_J2pJ4pJ6pJ2s, options, prms, Tp, 0);
        end
        
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
        
        if step_by_zd
            POs(solutionIndex,:) = [F_vec(1:5); zd0; F_vec(end)];
        elseif step_by_zd_withZFixed
            POs(solutionIndex,:) = [F_vec(1:2); 0; F_vec(3:4); zd0; F_vec(end)];
        elseif step_by_Tp_withZFixed
            POs(solutionIndex,:) = [F_vec(1:2); 0; F_vec(3:5); Tp];
        end
        
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
if n_Solutions > 1
    xlim([min(jacobiConstants), max(jacobiConstants)])
end

subplot(1,2,2); hold all
p1 = plot(jacobiConstants, abs(S1),'o','markeredgecolor',colors.blue,'markerfacecolor',colors.ltblue);
p2 = plot(jacobiConstants, abs(S2),'o','markeredgecolor',colors.red,'markerfacecolor',colors.ltred);
% plot(unique([min([jacobiConstants, jacobiConstants]) max([jacobiConstants, jacobiConstants])]),[2 2],'k','linewidth',1)
plot([min(jacobiConstants), max(jacobiConstants)],[2 2],'k','linewidth',1)
PlotBoi2('Jacobi Constant','',18,'LaTex')
ylim([1.9 2.1])
if n_Solutions > 1
    xlim([min(jacobiConstants), max(jacobiConstants)])
end









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
    famName_new       = [famName,'.CR3BP_J2pJ4pJ6pJ2s.L2_SHalo'];
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


if saveOtherHalos
    %%% Load the first family (northern or southern)
    family1_path = [savePath, 'Saturn_Enceladus.CR3BP_J2pJ4pJ6pJ2s.L2_SHalo.txt'];
    family1_data = dlmread(family1_path,',',1,0);
    
    %%% Get the header line from the first family
    fid_family1 = fopen(family1_path,'r');
    headerline = fgetl(fid_family1);
    fclose(fid_family1);
 
 
    % -------------------------------------------------
    %%% Preparing Save File for new family (southern or northern)
    % -------------------------------------------------
    fileAlreadyExists   = 1;
    fileVersion         = 0;
    newHaloFileName_new = ['Saturn_Enceladus.CR3BP_J2pJ4pJ6pJ2s.L2_NHalo'];
    while fileAlreadyExists == 1
        if isfile([savePath,newHaloFileName_new,'.txt']) == 1
            fileVersion = fileVersion + 1;
            newHaloFileName_new = sprintf('%s_%1d',newHaloFileName_new,fileVersion);
        else
            fileAlreadyExists = 0;
        end
    end
    newHaloFileName = [savePath,newHaloFileName_new,'.txt'];
    
    %%% Open File
    newHalofile = fopen(newHaloFileName,'wt');
   
    % -------------------------------------------------
    %%% Writing data
    % -------------------------------------------------
    %%% Write header
    fprintf(newHalofile,[headerline,'\n']);
    
    %%% Write data
    dataStart = 1;
    dataStop = size(family1_data,1);
    
    
    for kk = dataStart:dataStop
        fprintf(newHalofile,'%1.16f,%1.16f,%1.16f,%1.16f,%1.16f,%1.16f,%1.16f,%1.10f,%1.3f,%1.3f,%1.4f,%1.4f\n',...
            family1_data(kk,1), family1_data(kk,2), -family1_data(kk,3), family1_data(kk,4), family1_data(kk,5), -family1_data(kk,6),...
            family1_data(kk,7), family1_data(kk,8), family1_data(kk,9), family1_data(kk,10), family1_data(kk,11), family1_data(kk,12));
            
    end

    %%% Close file
    fclose(newHalofile);
    
    
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



X = [1.0005309373379350,-0.0067532866161135,-0.0000000000000000,-0.0133285109266791,0.0500383855202688,-0.0611264824419726]';
stm0_colVec = reshape(eye(6),36,1);
[~, XSTM_PO] = ode113(@Int_CR3BnSTM_J2pJ4pJ6pJ2s, [0, 1.9610787456760539], [X; stm0_colVec], options, prms);

figure
plot3(XSTM_PO(:,1),XSTM_PO(:,2),XSTM_PO(:,3),'r')



stm_tf_t0                           = reshape(XSTM_PO(end,7:42),6,6);
monodromy                           = stm_tf_t0;
[eigenVectors_new, eigenValues_new] = eig(monodromy);
[S1, S2]                            = getStabilityIndices(diag(eigenValues_new));











