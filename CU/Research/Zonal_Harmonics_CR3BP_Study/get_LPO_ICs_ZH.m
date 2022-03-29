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

%%% Periodic orbit ICs
PO_ICs = get_PO_ICs();

% ========================================================================
%%% Run Switches
% ========================================================================
run_Lyap = true;
run_Vert = true;

% ========================================================================
%%% Setup
% ========================================================================
% -------------------------------------------------
%%% System
% -------------------------------------------------
%%% Choose system
systemTag = 'Saturn_Enceladus';
% systemTag = 'Saturn_Enceladus';

%%% Set primary & secondary
[primary, secondary] = assignPrimaryAndSecondary_CR3BP(systemTag, bodies);

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
prms.n = secondary.meanMot * tNorm;
% prms.n = 1 + 3*(prms.J2p*prms.R1*prms.R1 + prms.J2s*prms.R2*prms.R2)/2;

% -------------------------------------------------
%%% Integration options
% -------------------------------------------------
%%% Integration options
tol = 1e-13;
options = odeset('RelTol',tol,'AbsTol',tol);

% -------------------------------------------------
%%% Multiple shooter options
% -------------------------------------------------
N_nodes = 3;
maxIter = 300;
errorTol = 1e-13;

% ========================================================================
%%% Dynamical Analysis
% ========================================================================
% -------------------------------------------------
%%% Analyze Libration points to get eigenvectors of state-dynamics matrix
% -------------------------------------------------
%%% Location of collinear libration points
rLPs_n_ZH = collinearEquilibriumPoints_ZH(prms);

%%% Jacobian dF/dX evaulated at each equilibrium point
A_L1 = get_Amat_CR3BP_J2pJ4pJ6pJ2s(prms, rLPs_n_ZH(1,:));
A_L2 = get_Amat_CR3BP_J2pJ4pJ6pJ2s(prms, rLPs_n_ZH(2,:));

%%% Acquire eigenvalues and eigenvectors 
[vMat1, eMat1] = eig(A_L1);
[vMat2, eMat2] = eig(A_L2);

e1 = diag(eMat1);
e2 = diag(eMat2);

%%% Choose an index to pick out an eigenvalue/eigenvector pair
indexLyap = 3;
indexVert = 5;

eVal2_Lyap = e2(indexLyap);
eVec2_Lyap = vMat2(:, indexLyap);

eVal2_Vert = e2(indexVert);
eVec2_Vert = vMat2(:, indexVert);

% ========================================================================
%%% Perturb Libration Points and correct to get Lyapunov PO
% ========================================================================
if run_Lyap
    % -------------------------------------------------
    %%% Perturb the Libration point and create nodes
    % -------------------------------------------------
    scaleDown_Lyap = 1e-7;

    time_guessLyap = [0, 3.5];

    X0_guessLyap_BaCR_n = [rLPs_n_ZH(2,:),zeros(1,3)]' + real(eVec2_Lyap).*scaleDown_Lyap;

    [~, X_nodes_Lyap] = get_nodes(X0_guessLyap_BaCR_n, time_guessLyap, N_nodes+1, @Int_CR3Bn_ZH, options, prms);

    % -------------------------------------------------
    %%% Create the free-variable vector
    % -------------------------------------------------
    F_vec_Lyap = zeros(N_nodes*6+1,1);
    for kk = 1:N_nodes
        F_vec_Lyap((kk*6-5):(kk*6)) = X_nodes_Lyap(kk,:)';
    end
    F_vec_Lyap(end) = time_guessLyap(end);
    F_vec_Lyap = F_vec_Lyap([5,7:end]);
    % F_vec = zeros(N_nodes*6,1);
    % for kk = 1:N_nodes
    %     F_vec((kk*6-5):(kk*6)) = X_nodes(kk,:)';
    % end


    % -------------------------------------------------
    %%% Initialize constrain error and loop over the linear update equation
    %%% with a multiple shooter
    % -------------------------------------------------
    c_norm_Lyap = 1000;
    iter = 0;

    while (c_norm_Lyap > errorTol) && (iter < maxIter)
        iter = iter + 1;

    %     [DF, constraints_vec] = multShooter_stateContinuity(N_nodes, F_vec, @Int_CR3BnSTM_J2pJ4pJ6pJ2s, options, prms);
    %     [DF, constraints_vec] = multShooter_stateContinuity_TpFixed(N_nodes, F_vec, @Int_CR3BnSTM_J2pJ4pJ6pJ2s, options, prms, time_guess(end));
        [DF, constraints_vec] = multShooter_LyapPO_stateContinuity(N_nodes, X0_guessLyap_BaCR_n(1), F_vec_Lyap, @Int_CR3BnSTM_J2pJ4pJ6pJ2s, options, prms);

        c_norm_Lyap = norm(constraints_vec);

        if (c_norm_Lyap > errorTol)
            warning('off','MATLAB:nearlySingularMatrix')
            F_vec_Lyap = F_vec_Lyap - DF'*((DF*(DF'))\constraints_vec);
            warning('on','MATLAB:nearlySingularMatrix')
        end




    end

    if c_norm_Lyap > errorTol
        warning('Solution not converged on')
    else
        PO_IC_Lyap = [X0_guessLyap_BaCR_n(1); 0;0;0;F_vec_Lyap(1); 0; F_vec_Lyap(end)];

        [t_Lyap, X_Lyap] = ode113(@Int_CR3Bn_ZH, [0, PO_IC_Lyap(end)], PO_IC_Lyap(1:6), options, prms);

        figure; hold all
        plot3(X_Lyap(:,1),X_Lyap(:,2),X_Lyap(:,3))
        PlotBoi3_CR3Bn(20)

        fprintf('Lyapunov IC:\n')
        prettyColVec(PO_IC_Lyap)
    end


end % if run_Lyap










% ========================================================================
%%% Perturb Libration Points and correct to get Vertical PO
% ========================================================================
if run_Vert
    % -------------------------------------------------
    %%% Perturb the Libration point and create nodes
    % -------------------------------------------------
    scaleDown_Vert = 1e-4;

    time_guessVert = [0, 3];

    X0_guessVert_BaCR_n = [rLPs_n_ZH(2,:),zeros(1,3)]' + real(eVec2_Vert).*scaleDown_Vert;

    [~, X_nodes_Vert] = get_nodes(X0_guessVert_BaCR_n, time_guessVert, N_nodes+1, @Int_CR3Bn_ZH, options, prms);

    % -------------------------------------------------
    %%% Create the free-variable vector
    % -------------------------------------------------
    F_vec_Vert = zeros(N_nodes*6,1);
    F_vec_Vert(1:5) = X0_guessVert_BaCR_n(1:5);
    for kk = 2:N_nodes
        F_vec_Vert((kk*6-6):(kk*6-1)) = X_nodes_Vert(kk,:)';
    end
    F_vec_Vert(end) = time_guessVert(end);

    % -------------------------------------------------
    %%% Initialize constrain error and loop over the linear update equation
    %%% with a multiple shooter
    % -------------------------------------------------
    c_norm_Vert = 1000;
    iter = 0;

    while (c_norm_Vert > errorTol) && (iter < maxIter)
        iter = iter + 1;

%         [DF, constraints_vec] = multShooter_stateContinuity(N_nodes, F_vec_Vert, @Int_CR3BnSTM_J2pJ4pJ6pJ2s, options, prms);
        [DF, constraints_vec] = multShooter_stateContinuity_zdFixed(N_nodes, X0_guessVert_BaCR_n(6),F_vec_Vert, @Int_CR3BnSTM_J2pJ4pJ6pJ2s, options, prms);

        c_norm_Vert = norm(constraints_vec);

        if (c_norm_Vert > errorTol)
            warning('off','MATLAB:nearlySingularMatrix')
            F_vec_Vert = F_vec_Vert - DF'*((DF*(DF'))\constraints_vec);
            warning('on','MATLAB:nearlySingularMatrix')
        end




    end

    if c_norm_Vert > errorTol
        warning('Solution not converged on')
    else
%         PO_IC_Vert = [X0_guessVert_BaCR_n(1); 0;0;0;F_vec_Vert(1); 0; F_vec_Vert(end)];
%         PO_IC_Vert = [F_vec_Vert(1:6); F_vec_Vert(end)];
        PO_IC_Vert = [F_vec_Vert(1:5); X0_guessVert_BaCR_n(6); F_vec_Vert(end)];

        [t_Vert, X_Vert] = ode113(@Int_CR3Bn_ZH, [0, PO_IC_Vert(end)], PO_IC_Vert(1:6), options, prms);

        figure; hold all
        plot3(X_Vert(:,1),X_Vert(:,2),X_Vert(:,3))
        PlotBoi3_CR3Bn(20)

        fprintf('Vertical IC:\n')
        prettyColVec(PO_IC_Vert)
    end

end %if run_Vert
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
















