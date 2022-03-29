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

% ========================================================================
%%% Setup
% ========================================================================
% -------------------------------------------------
%%% System
% -------------------------------------------------
%%% Choose system
% systemTag = 'Jupiter_Europa';
% PO_0 = [1.0169952731165666,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0199746571663958,0.0000000000000000,3.1143692478140430];

systemTag = 'Saturn_Enceladus';
PO_0 = [1.0034009101733159,0.0000000000000000,0.0000000000000000,0.0000000000000000,0.0038595783739539,0.0000000000000000,3.1566436907700011];

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
N_nodes = 4;
maxIter = 300;
errorTol = 1e-13;

z_new = 1e-5;
zd_new = 0;
% ========================================================================
%%% Perturb Conditions and perturb to get halo
% ========================================================================
% -------------------------------------------------
%%% Get nodes from initial PO
% -------------------------------------------------
[~, X_nodes] = get_nodes(PO_0(1:6)', [0, PO_0(7)], N_nodes+1, @Int_CR3Bn_ZH, options, prms);

% -------------------------------------------------
%%% Create the free-variable vector
% -------------------------------------------------
F_vec = zeros(N_nodes*6+1,1);
for kk = 1:N_nodes
    F_vec((kk*6-5):(kk*6)) = X_nodes(kk,:)';
end
F_vec(end) = PO_0(7);
F_vec = F_vec([1:2,4:5, 7:end]);

% -------------------------------------------------
%%% Initialize constrain error and loop over the linear update equation
%%% with a multiple shooter
% -------------------------------------------------
c_norm = 1000;
iter = 0;

while (c_norm > errorTol) && (iter < maxIter)
    iter = iter + 1;

    [DF, constraints_vec] = multShooter_stateContinuity_zzdFixed(N_nodes, z_new, zd_new, F_vec, @Int_CR3BnSTM_J2pJ4pJ6pJ2s, options, prms);
    
    c_norm = norm(constraints_vec);

    if (c_norm > errorTol)
        warning('off','MATLAB:nearlySingularMatrix')
        F_vec = F_vec - DF'*((DF*(DF'))\constraints_vec);
        warning('on','MATLAB:nearlySingularMatrix')
    end




end

if c_norm > errorTol
    warning('Solution not converged on')
else
    PO_IC_new = [F_vec(1:2); z_new; F_vec(3:4); zd_new; F_vec(end)];

    [~, X_new] = ode113(@Int_CR3Bn_ZH, [0, PO_IC_new(end)], PO_IC_new(1:6), options, prms);
    

    figure; hold all
    plot3(X_new(:,1),X_new(:,2),X_new(:,3))
    PlotBoi3_CR3Bn(20)

    fprintf('Halo IC:\n')
    prettyColVec(PO_IC_new)
    
    
    options_zEqualsZero = odeset('Event',@event_zEqualsZero,'RelTol',tol,'AbsTol',tol);
    [~, X_newWithEvent, ~, X_event, ~] = ode113(@Int_CR3Bn_ZH, [0, PO_IC_new(end)], PO_IC_new(1:6), options_zEqualsZero, prms);
    figure; hold all
    plot3(X_newWithEvent(:,1),X_newWithEvent(:,2),X_newWithEvent(:,3))
    PlotBoi3_CR3Bn(20)

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
















