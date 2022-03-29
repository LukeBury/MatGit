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
% clear
% clc
% close all
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
fix_nothing          = false;
fix_JC               = false;
fix_Tp               = true;
fix_zd               = false;
fix_z_zd             = false;
fix_z_Tp             = false;
fix_ydAndTpFree_Lyap = false;

% ========================================================================
%%% Setup
% ========================================================================
% -------------------------------------------------
%%% Choose PO family
% -------------------------------------------------
famName = 'Jupiter_Europa';
PO_guess = [1.018461701526617;
 0.000000000000000;
 0.000000000000000;
 0.000000000000000;
 -0.004241100000000;
 0.000000000000000;
 7.395200000000000];


%%%% Set primary & secondary
[primary, secondary] = assignPrimaryAndSecondary_CR3BP(famName, bodies);

% ------------------------------------------------- 
%%% Shooter options
% -------------------------------------------------
%%% Number of nodes for multiple shooter
N_nodes = 36; 

%%% Error tolerance for constraint vector in multiple shooter
error_tol = 1e-13; 

%%% Maximum number of multiple-shooter iterations 
maxIter = 200;

%%% Scale the step size (<1)
stepSizeMultiplier = 0.5;
% -------------------------------------------------
%%% System
% -------------------------------------------------
%%% Factor for normalizing distances
rNorm = secondary.a; % n <-> km
tNorm = 1/secondary.meanMot;
vNorm = rNorm/tNorm;

%%% Setting parameters structure
prms.u = secondary.MR;
prms.R1 = primary.R / rNorm;
prms.R2 = secondary.R_n;


%%% Getting normalized mean motion
prms.n = 1;

%%% Equillibrium Points
rLPs_n = collinearEquilibriumPoints_ZH(prms);

% ------------------------------------------------- 
%%% Integration Options
% -------------------------------------------------
tol     = 1e-13;
options = odeset('RelTol',tol,'AbsTol',tol);

% ========================================================================
%%% Find PO
% ========================================================================
%%% Establish guess
t_guess = PO_guess(end);
X0_guess = PO_guess(1:6);

% -------------------------------------------------
%%% Plot Guess
% -------------------------------------------------
[~, X_guess] = ode113(@Int_CR3Bn, [0, t_guess], X0_guess, options, prms);

figure(555); hold all
title('Guess')
PlotBoi3_CR3Bn(23)
plot3(X_guess(:,1),X_guess(:,2),X_guess(:,3),'r')

% -------------------------------------------------
%%% Get nodes
% -------------------------------------------------
[~, X_nodes] = get_nodes(X0_guess', [0, t_guess], N_nodes+1, @Int_CR3Bn, options, prms);

% -------------------------------------------------
%%% Create the free-variable vector
% -------------------------------------------------
if fix_nothing
    F_vec = zeros(N_nodes*6+1,1);
    for kk = 1:N_nodes
        F_vec((kk*6-5):(kk*6)) = X_nodes(kk,:)';
    end
    F_vec(end) = t_guess;
    
elseif fix_JC
    JC0 = getJacobiConstant_ZH(PO_guess(1:6)',prms);
    
    F_vec = zeros(N_nodes*6+1,1);
    for kk = 1:N_nodes
        F_vec((kk*6-5):(kk*6)) = X_nodes(kk,:)';
    end
    F_vec(end) = t_guess;
    
elseif fix_Tp
    F_vec = zeros(N_nodes*6,1);
    for kk = 1:N_nodes
        F_vec((kk*6-5):(kk*6)) = X_nodes(kk,:)';
    end

elseif fix_zd
    zd0 = PO_guess(6);

    F_vec = zeros(N_nodes*6+1,1);
    for kk = 1:N_nodes
        F_vec((kk*6-5):(kk*6)) = X_nodes(kk,:)';
    end
    F_vec(end) = t_guess;
%         F_vec = F_vec([1:2,4:5, 7:end]);
    F_vec = F_vec([1:5, 7:end]);

elseif fix_z_zd
    z0  = PO_guess(3);
    zd0 = PO_guess(6);

    F_vec = zeros(N_nodes*6+1,1);
    for kk = 1:N_nodes
        F_vec((kk*6-5):(kk*6)) = X_nodes(kk,:)';
    end
    F_vec(end) = t_guess;
    F_vec = F_vec([1:2,4:5, 7:end]);

elseif fix_z_Tp
    z0  = PO_guess(3);

    F_vec = zeros(N_nodes*6,1);
    for kk = 1:N_nodes
        F_vec((kk*6-5):(kk*6)) = X_nodes(kk,:)';
    end
    F_vec = F_vec([1,2,4:end]);
    
elseif fix_ydAndTpFree_Lyap
    x0 = PO_guess(1);
    
    F_vec = zeros(N_nodes*6+1,1);
    for kk = 1:N_nodes
        F_vec((kk*6-5):(kk*6)) = X_nodes(kk,:)';
    end
    F_vec(end) = t_guess;
    F_vec = F_vec([5, 7:end]);
end
    
    
    
    
    
    
% -------------------------------------------------
%%% Initialize constrain error and loop over the linear update equation
%%% with a multiple shooter
% -------------------------------------------------
c_norm = 1000;
iter = 0;

% figure; hold all
% PlotBoi3_CR3Bn(23)
% tol = 1e-13;
% options = odeset('RelTol',tol,'AbsTol',tol);

while (c_norm > error_tol) && (iter < maxIter)
    iter = iter + 1;
    
    if ((iter > 1) && (iter <= 15)) || (mod(iter,50) == 0)
        fprintf('Error %1.1e ... iteration: %1d\n', c_norm, iter)
    end
    
    if fix_nothing
        [DF, constraints_vec] = multShooter_stateContinuity(N_nodes, F_vec, @Int_CR3BnSTM, options, prms);
    elseif fix_JC
        [DF, constraints_vec] = multShooter_stateContinuity_JCFixed(N_nodes, F_vec, @Int_CR3BnSTM, options, prms, JC0);
    elseif fix_Tp
        [DF, constraints_vec] = multShooter_stateContinuity_TpFixed(N_nodes, F_vec, @Int_CR3BnSTM, options, prms, t_guess);
    elseif fix_zd
        [DF, constraints_vec] = multShooter_stateContinuity_zdFixed(N_nodes, zd0, F_vec, @Int_CR3BnSTM, options, prms);
    elseif fix_z_zd
        [DF, constraints_vec] = multShooter_stateContinuity_zzdFixed(N_nodes, z0, zd0, F_vec, @Int_CR3BnSTM, options, prms);
    elseif fix_z_Tp
        [DF, constraints_vec] = multShooter_stateContinuity_zTpFixed(N_nodes, F_vec, @Int_CR3BnSTM, options, prms, t_guess, z0);
    elseif fix_ydAndTpFree_Lyap
        [DF, constraints_vec] = multShooter_stateContinuity_ydFree(N_nodes, x0, F_vec, @Int_CR3BnSTM, options, prms);
    end

    c_norm = norm(constraints_vec);

    if (c_norm > error_tol)
        warning('off','MATLAB:nearlySingularMatrix')
        F_vec = F_vec - (stepSizeMultiplier)*DF'*((DF*(DF'))\constraints_vec);
        warning('on','MATLAB:nearlySingularMatrix')
    end
    
    
    

%     X = F_vec([1:6]);
%     [T_test, X_test] = ode113(@Int_CR3Bn, [0 t_guess], X, options, prms);
%     plot3(X_test(:,1),X_test(:,2),X_test(:,3))
%     drawnow
    
end

fprintf('-----------------------------------\n')
fprintf('Error %1.1e ... iterations: %1d\n', c_norm, iter)
fprintf('Elapsed time: %1.4f seconds\n',toc(ticWhole))

%%% Assign PO
if fix_nothing
    PO = F_vec([1:6,end]);
elseif fix_JC
    PO = F_vec([1:6,end]);
elseif fix_Tp
    PO = [F_vec(1:6); t_guess];
elseif fix_zd
    PO = [F_vec(1:5); zd0; F_vec(end)];
elseif fix_z_zd
    PO = [F_vec(1:2); z0; F_vec(3:4); zd0; F_vec(end)];
elseif fix_z_Tp
    PO = [F_vec(1:2); z0; F_vec(3:5); t_guess];
elseif fix_ydAndTpFree_Lyap
    PO = [x0;0;0;0;F_vec(1);0;F_vec(end)];
end


prettyColVec(PO)

stm0_colVec = reshape(eye(6),36,1);
[~, XSTM_PO] = ode113(@Int_CR3BnSTM, [0, PO(7)], [PO(1:6); stm0_colVec], options, prms);


figure(556); hold all
PlotBoi3_CR3Bn(23)
plot3(XSTM_PO(:,1),XSTM_PO(:,2),XSTM_PO(:,3),'r')

stm_tf_t0                           = reshape(XSTM_PO(end,7:42),6,6);
monodromy                           = stm_tf_t0;
[eigenVectors_new, eigenValues_new] = eig(monodromy);
[S1, S2]                            = getStabilityIndices(diag(eigenValues_new));

JC                       = getJacobiConstant_ZH(PO(1:6)', prms);
L2FlythroughVelocity_mps = JC_2_L2FlyoverVelocity(JC, prms, rLPs_n(2,:), vNorm);
landingVelocity_mps      = JC_2_approxLandingVelocity(JC, prms, vNorm);
stabilityIndices         = [S1, S2];

       








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











