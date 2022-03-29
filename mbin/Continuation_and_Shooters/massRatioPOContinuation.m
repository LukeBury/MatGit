% ========================================================================
%%% Description
% ========================================================================
% Use a multiple shooter to find a recurrence of a PO at a new mass ratio

% Created: 01/05/21
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
shooter_allFree         = false; % 
shooter_y0Fixed         = false;
shooter_y0xd0zd0Fixed   = false;

shooter_TpFixed         = true;  %  
shooter_y0TpFixed       = false; %  
shooter_y0xd0zd0TpFixed = false;
    extraTpScale = 1+2e-5;

plot_TpScalingFormula = false;


% 
%     F_new = [F_new(1); F_new(3); F_new(5); F_new(7:end-1)];

% ========================================================================
%%% Setup
% ========================================================================
% -------------------------------------------------
%%% Inputs / ICs
% -------------------------------------------------
%%% Primary_Secondary system
% systemName = 'Jupiter_Europa';
systemName = 'Saturn_Enceladus';

% %%% PO (from europa SHalo)
% PO_0 = [1.016439777885220;
%  0.000000000000000;
%  -0.023440560089874;
%  0.000000000000000;
%  -0.030460921895392;
%  -0.000000000000001;
%  2.685055450126591];

% PO_0 = [1.016655732932955; % from L2 L T 3P2
%  0.000000000000000;
%  -0.023244074523533;
%  -0.000032010806411;
%  -0.030536812091331;
%  -0.000610023587066;
%  5.416160835728377];

% Enc try 1
% PO_0 = [0.9878707573684111;
%  0.0000000000000000;
%  -0.0166948937322335;
%  0.0000000000000000;
%  -0.0056182441223105;
%  0.0000000000000000;
%  4.4369646720470675];



% PO_0 = [0.9996165286614656;
%  0.0000000000000000;
%  -0.0009890792364341;
%  0.0000000000000000;
%  -0.0161973279788486;
%  0.0000000000000000;
%  4.8871289984803123];
PO_0 = [0.9992979397836381;
 0.0000000000000000;
 -0.0023476139142152;
 -0.0000000000000050;
 -0.0087808958860086;
 0.0000000000000033;
 5.9565678747052306];

step1 = 1;
% -------------------------------------------------
%%% Parameter setup
% -------------------------------------------------
% --------------------------
%%% Shooter options
% --------------------------
% n_Nodes    = 6;
n_Nodes    = 13;
error_tol  = 1e-11;
% error_tol  = 1e-9;989
ms_iterMax = 1200;
stepSize   = 1;
% --------------------------
%%% Set primary & secondary
% --------------------------
[primary, secondary] = assignPrimaryAndSecondary_CR3BP(systemName, bodies);

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
% prms.R2 = secondary.R_n;

%%% Mass ratios
mu_1 = bodies.enceladus.MR;
% mu_2 = bodies.ganymede.MR;
mu_2 = 1.898884589251784e-07;
% mu_2 = bodies.moon.MR;

%%% Choose how many steps of natural parameter continuation to go from mu_1
%%% to mu_2. Vectors of these steps are generated with both linear and
%%% logarithmic spacing and combined to help avoid disproportional gaps 
n_mu = 6
% n_mu = 1000;
% n_mu = 2000;
% n_mu = 5000;
% n_mu = 10000;
% mu_vec = linspace(mu_1, mu_2, n_mu);
mu_vec1 = logspace(log10(mu_1), log10(mu_2), (n_mu/2+2));
mu_vec2 = linspace(mu_1, mu_2, n_mu/2);
mu_vec = sort([mu_vec1(2:end-1), mu_vec2], 'descend');

%%% Initialize STM as column vector to be added with state
stm0_colVec = reshape(eye(6),36,1);

% ------------------------------------------------- 
%%% Integration Options
% -------------------------------------------------
tol = 1e-13;
options = odeset('RelTol',tol,'AbsTol',tol);
options_yEquals0_terminal = odeset('Event',@event_yEqualsZero,'RelTol',tol,'AbsTol',tol);

% ------------------------------------------------- 
%%% Integrate forward and backward to nearest y=0 crossings so the shortest
%%% of the two can be used - this ensures a guess closest to the initial
% ------------------------------------------------- 
if 1+1==1
    %%% Integrate forward and backward
    [Tfix_fwd, Xfix_fwd, ~, ~, ~] = ode113(@Int_CR3BnSTM, [0, PO_0(7)], [PO_0(1:6); reshape(eye(6),36,1)], options_yEquals0_terminal, prms);
    [Tfix_bkwd, Xfix_bkwd, ~, ~, ~] = ode113(@Int_CR3BnSTM, [PO_0(7), 0], [PO_0(1:6); reshape(eye(6),36,1)], options_yEquals0_terminal, prms);

    %%% Times to y=0
    dt_fwd = abs(Tfix_fwd(end) - Tfix_fwd(1));
    dt_bkwd = abs(Tfix_bkwd(end) - Tfix_bkwd(1));

    %%% Find minmum of two times and use that state as new guess
    if dt_fwd < dt_bkwd
        PO_0_guess = [Xfix_fwd(end,1:6)'; PO_0(7)];
    elseif dt_fwd >= dt_fwd
        PO_0_guess = [Xfix_bkwd(end,1:6)'; PO_0(7)];
    end
    figure(333); hold all
    p1 = plot3(Xfix_fwd(:,1),Xfix_fwd(:,2),Xfix_fwd(:,3),'b');
    plot3(Xfix_fwd(1,1),Xfix_fwd(1,2),Xfix_fwd(1,3),'bo')
    plot3(Xfix_fwd(end,1),Xfix_fwd(end,2),Xfix_fwd(end,3),'bx')
    
    p2 = plot3(Xfix_bkwd(:,1),Xfix_bkwd(:,2),Xfix_bkwd(:,3),'r');
    plot3(Xfix_bkwd(1,1),Xfix_bkwd(1,2),Xfix_bkwd(1,3),'ro')
    plot3(Xfix_bkwd(end,1),Xfix_bkwd(end,2),Xfix_bkwd(end,3),'rx')
    
    legend([p1 p2], 'fwd', 'bkwd')
    
%     PO_0 = [Xfix_fwd(end,1:6)'; PO_0(7)]

end

%%% Integrate initial PO
[T_PO_0, X_PO_0] = ode113(@Int_CR3BnSTM, linspace(0, PO_0(7), 10000), [PO_0(1:6); stm0_colVec], options, prms);


% -------------------------------------------------
%%% Shooter-specific setup
% -------------------------------------------------
if shooter_allFree
    
elseif shooter_y0Fixed
    
elseif shooter_y0xd0zd0Fixed
    if PO_0(4) ~= 0
        warning('Switching xd0 to 0')
        PO_0(4) = 0;
    end
    if PO_0(6) ~= 0
        warning('Switching zd0 to 0')
        PO_0(6) = 0;
    end
end

% ------------------------------------------------- 
%%% Loop through mass ratio values
% -------------------------------------------------

POs = zeros(n_mu, 7);
stabilityIndices = zeros(n_mu,2);

PO_new = PO_0;

if shooter_TpFixed || shooter_y0TpFixed || shooter_y0xd0zd0TpFixed
%     x = log10(prms.u);
%     p = [0.001017467385714, 0.022429766590012, 0.191742044946650, 0.769927745899565, 4.288444369877601];
%     Tp_fundamentalLyap_initial = p(1)*x^4 + p(2)*x^3 + p(3)*x^2 + p(4)*x + p(5);
%     Tp0 = PO_0(7);
    prms.u = mu_vec(step1);
    rLPs_n = EquilibriumPoints(prms.u, 1);
    A_mat_Lp = get_Amat_CR3BP(prms.u, rLPs_n(2,:), 1);
    [Tp_fundamentalLyap_last] = get_fundamentalLyapunovTp(A_mat_Lp);
end

%%% Prep for warning suppression
warning('off','MATLAB:nearlySingularMatrix')

figure; hold all
for index = step1:n_mu
% for index = 110:n_mu
    prms.u = mu_vec(index);
    
    PO_guess = PO_new;
    
    
    if shooter_TpFixed || shooter_y0TpFixed || shooter_y0xd0zd0TpFixed
        %%% figure out time scaling
%         x = log10(prms.u);
%         Tp_fundamentalLyap_new = p(1)*x^4 + p(2)*x^3 + p(3)*x^2 + p(4)*x + p(5);
%         T_scaling = Tp_fundamentalLyap_new / Tp_fundamentalLyap_initial;
%         PO_guess(7) = Tp0 * T_scaling;
        
        rLPs_n = EquilibriumPoints(prms.u, 1);
        A_mat_Lp = get_Amat_CR3BP(prms.u, rLPs_n(2,:), 1);
        [Tp_fundamentalLyap_current] = get_fundamentalLyapunovTp(A_mat_Lp);
        
        scale_Tp = Tp_fundamentalLyap_current / Tp_fundamentalLyap_last;
        
        PO_guess(7) = PO_guess(7)*scale_Tp*extraTpScale;
        

    end
    
    if shooter_allFree
        [PO_new, counter, constraint_error] = correctPO_multShooter_stateContinuity(PO_guess, n_Nodes, @Int_CR3BnSTM, options, prms, error_tol, ms_iterMax, stepSize);
    elseif shooter_y0Fixed
        [PO_new, counter, constraint_error] = correctPO_multShooter_stateContinuity_y0Fixed(PO_guess, n_Nodes, @Int_CR3BnSTM, options, prms, error_tol, ms_iterMax, stepSize, PO_guess(2));
    elseif shooter_y0xd0zd0Fixed
        [PO_new, counter, constraint_error] = correctPO_multShooter_stateContinuity_y0xd0zd0Fixed(PO_guess, n_Nodes, @Int_CR3BnSTM, options, prms, error_tol, ms_iterMax, stepSize, PO_guess(2), PO_guess(4), PO_guess(6));
    elseif shooter_TpFixed
        [PO_new, counter, constraint_error] = correctPO_multShooter_stateContinuity_TpFixed(PO_guess, n_Nodes, @Int_CR3BnSTM, options, prms, error_tol, ms_iterMax, stepSize, PO_guess(7));
    elseif shooter_y0TpFixed
        [PO_new, counter, constraint_error] = correctPO_multShooter_stateContinuity_y0TpFixed(PO_guess, n_Nodes, @Int_CR3BnSTM, options, prms, error_tol, ms_iterMax, stepSize, PO_guess(2), PO_guess(7));
    elseif shooter_y0xd0zd0TpFixed
        [PO_new, counter, constraint_error] = correctPO_mS_sC_y0xd0zd0TpFixed(PO_guess, n_Nodes, @Int_CR3BnSTM, options, prms, error_tol, ms_iterMax, stepSize, PO_guess(2), PO_guess(4), PO_guess(6), PO_guess(7));
    end
    
    POs(index,:) = PO_new';
    
    %%% Integrate and plot new solution
    [T_current, X_current] = ode113(@Int_CR3BnSTM, [0, PO_new(7)], [PO_new(1:6); stm0_colVec], options, prms);
    
    %%% Check stability of new solution
    stm_tf_t0                           = reshape(X_current(end,7:42),6,6);
    monodromy                           = stm_tf_t0;
    [eigenVectors_new, eigenValues_new] = eig(monodromy);
    [S1, S2]                            = getStabilityIndices(diag(eigenValues_new));
    stabilityIndices(index, :) = [S1, S2];
    
    plot3(X_current(:,1),X_current(:,2),X_current(:,3),'b')
    drawnow
    
    if shooter_TpFixed || shooter_y0TpFixed || shooter_y0xd0zd0TpFixed
        Tp_fundamentalLyap_last = Tp_fundamentalLyap_current;
    end
    
    prop_error = norm(X_current(end,1:6) - X_current(1,1:6)) / norm(X_current(1,1:6));
    fprintf('%3d ... constr. error: %1.1e ... prop error: %1.1e ... iter: %3d ... time: %1.1f\n', index, constraint_error, prop_error, counter, toc(ticWhole))
    989;
    
end

%%% Prep for warning suppression
warning('on','MATLAB:nearlySingularMatrix')

%%% Integrate the new guess
[T_newPO, X_newPO] = ode113(@Int_CR3BnSTM, linspace(0, POs(end,7), 10000), [POs(end,1:6)'; stm0_colVec], options, prms);

%%% Plot the results
figure; hold all
p1 = plot3(X_PO_0(:,1),X_PO_0(:,2),X_PO_0(:,3),'color', colors.blue2,'linewidth',2);
p2 = plot3(X_newPO(:,1),X_newPO(:,2),X_newPO(:,3),'color', colors.red2,'linewidth',2);
PlotBoi3_CR3Bn(26)
legend([p1, p2], 'Old PO', 'New PO','FontSize',16)

prettyColVec(POs(end,:))


stm_tf_t0                           = reshape(X_newPO(end,7:42),6,6);
monodromy                           = stm_tf_t0;
[eigenVectors_new, eigenValues_new] = eig(monodromy);
[S1, S2]                            = getStabilityIndices(diag(eigenValues_new));
[alpha, beta]                       = getBrouckeStabilityParameters(diag(eigenValues_new), monodromy);













if plot_TpScalingFormula
    % -------------------------------------------------
    %%% Looking at fit between base LyapTp and log10(MR)
    % -------------------------------------------------
    %%%[Tp Lyap1, Tp sys, MR]
%     TpLyap_Tp_MR = [bodies.enceladus.LyapTp, bodies.enceladus.Tp, bodies.enceladus.MR;...%Saturn Enceladus
%     bodies.europa.LyapTp, bodies.europa.Tp, bodies.europa.MR;...%Jupiter Europa
%     bodies.ganymede.LyapTp, bodies.ganymede.Tp, bodies.ganymede.MR;...%Jupiter Ganymede
%     bodies.triton.LyapTp, bodies.triton.Tp, bodies.triton.MR;...% Neptune Triton
%     bodies.titan.LyapTp, bodies.titan.Tp, bodies.titan.MR;...% Saturn Titan
%     bodies.moon.LyapTp, bodies.moon.Tp, bodies.moon.MR]; % Earth Moon

    MRs     = [bodies.enceladus.MR,     bodies.europa.MR,     bodies.ganymede.MR,     bodies.triton.MR,     bodies.titan.MR,     bodies.moon.MR];
    LyapTps = [bodies.enceladus.LyapTp, bodies.europa.LyapTp, bodies.ganymede.LyapTp, bodies.triton.LyapTp, bodies.titan.LyapTp, bodies.moon.LyapTp];
    [p, S] = polyfit(log10(MRs), LyapTps, 4);
%     p = [0.001017467385714, 0.022429766590012, 0.191742044946650, 0.769927745899565, 4.288444369877601];

    xs = linspace(min(log10(MRs)), max(log10(MRs)), 1000);

    yFit_Lyap = p(1).*xs.^4 + p(2).*xs.^3 + p(3).*xs.^2 + p(4).*xs + p(5);

    figure; hold all
    plot(xs, yFit_Lyap, 'r','linewidth', 2)
%     p_enceleadus = plot(log10(TpLyap_Tp_MR(1,3)), TpLyap_Tp_MR(1,1), '.', 'markersize', 35);
%     p_europa = plot(log10(TpLyap_Tp_MR(2,3)), TpLyap_Tp_MR(2,1), '.', 'markersize', 35);
%     p_ganymede = plot(log10(TpLyap_Tp_MR(3,3)), TpLyap_Tp_MR(3,1), '.', 'markersize', 35);
%     p_triton = plot(log10(TpLyap_Tp_MR(4,3)), TpLyap_Tp_MR(4,1), '.', 'markersize', 35);
%     p_titan = plot(log10(TpLyap_Tp_MR(5,3)), TpLyap_Tp_MR(5,1), '.', 'markersize', 35);
%     p_moon = plot(log10(TpLyap_Tp_MR(6,3)), TpLyap_Tp_MR(6,1), '.', 'markersize', 35);
    p_enceleadus = plot(log10(bodies.enceladus.MR), bodies.enceladus.LyapTp, '.', 'markersize', 35);
    p_europa     = plot(log10(bodies.europa.MR),    bodies.europa.LyapTp,    '.', 'markersize', 35);
    p_ganymede   = plot(log10(bodies.ganymede.MR),  bodies.ganymede.LyapTp,  '.', 'markersize', 35);
    p_triton     = plot(log10(bodies.triton.MR),    bodies.triton.LyapTp,    '.', 'markersize', 35);
    p_titan      = plot(log10(bodies.titan.MR),     bodies.titan.LyapTp,     '.', 'markersize', 35);
    p_moon       = plot(log10(bodies.moon.MR),      bodies.moon.LyapTp,      '.', 'markersize', 35);
    
    PlotBoi2('$log_{10}(\mu)$','$T_{P_{Lyap}}$',26, 'LaTex')
    legend([p_enceleadus, p_europa, p_ganymede, p_triton, p_titan, p_moon], 'Enceladus', 'Europa', 'Ganymede', 'Triton', 'Titan', 'Moon','location', 'northwest','FontSize',18)
    
    f = polyval(p,log10(MRs));
    T = table(log10(MRs),LyapTps,f,LyapTps-f,'VariableNames',{'log10(MR)','Tp_Lyap','Fit','FitError'});
    title(sprintf('Mean/Max Error: %1.1e / %1.1e', mean(abs(T.FitError)), max(abs(T.FitError))))

    % -------------------------------------------------
    %%% Look at how well the fit predicts halo/lyap bifurcation points
    % -------------------------------------------------
    bodies.enceladus.LyapHaloTp = 3.08985242;
    bodies.europa.LyapHaloTp    = 3.12421243;
    bodies.ganymede.LyapHaloTp  = 3.14356667;
    bodies.triton.LyapHaloTp    = 3.16753961;
    bodies.titan.LyapHaloTp     = 3.17115669;
    bodies.moon.LyapHaloTp      = 3.41551092;
    
    LyapHaloTp_referenceBody = bodies.titan;
    bodies.enceladus.LyapHaloTp_pred = LyapHaloTp_referenceBody.LyapHaloTp * (bodies.enceladus.LyapTp / LyapHaloTp_referenceBody.LyapTp);
    bodies.europa.LyapHaloTp_pred    = LyapHaloTp_referenceBody.LyapHaloTp * (bodies.europa.LyapTp    / LyapHaloTp_referenceBody.LyapTp);
    bodies.ganymede.LyapHaloTp_pred  = LyapHaloTp_referenceBody.LyapHaloTp * (bodies.ganymede.LyapTp  / LyapHaloTp_referenceBody.LyapTp);
    bodies.triton.LyapHaloTp_pred    = LyapHaloTp_referenceBody.LyapHaloTp * (bodies.triton.LyapTp    / LyapHaloTp_referenceBody.LyapTp);
    bodies.titan.LyapHaloTp_pred     = LyapHaloTp_referenceBody.LyapHaloTp * (bodies.titan.LyapTp     / LyapHaloTp_referenceBody.LyapTp);
    bodies.moon.LyapHaloTp_pred      = LyapHaloTp_referenceBody.LyapHaloTp * (bodies.moon.LyapTp      / LyapHaloTp_referenceBody.LyapTp);
    
    figure; hold all
    
    p_enceleadus = plot(log10(bodies.enceladus.MR), bodies.enceladus.LyapHaloTp, '.', 'markersize', 35);
    p_europa     = plot(log10(bodies.europa.MR),    bodies.europa.LyapHaloTp,    '.', 'markersize', 35);
    p_ganymede   = plot(log10(bodies.ganymede.MR),  bodies.ganymede.LyapHaloTp,  '.', 'markersize', 35);
    p_triton     = plot(log10(bodies.triton.MR),    bodies.triton.LyapHaloTp,    '.', 'markersize', 35);
    p_titan      = plot(log10(bodies.titan.MR),     bodies.titan.LyapHaloTp,     '.', 'markersize', 35);
    p_moon       = plot(log10(bodies.moon.MR),      bodies.moon.LyapHaloTp,      '.', 'markersize', 35);
    
    p_enceleadus_pred = plot(log10(bodies.enceladus.MR), bodies.enceladus.LyapHaloTp_pred, '+', 'markersize', 18, 'markeredgecolor','k');
    p_europa_pred     = plot(log10(bodies.europa.MR),    bodies.europa.LyapHaloTp_pred,    '+', 'markersize', 18, 'markeredgecolor','k');
    p_ganymede_pred   = plot(log10(bodies.ganymede.MR),  bodies.ganymede.LyapHaloTp_pred,  '+', 'markersize', 18, 'markeredgecolor','k');
    p_triton_pred     = plot(log10(bodies.triton.MR),    bodies.triton.LyapHaloTp_pred,    '+', 'markersize', 18, 'markeredgecolor','k');
    p_titan_pred      = plot(log10(bodies.titan.MR),     bodies.titan.LyapHaloTp_pred,     '+', 'markersize', 18, 'markeredgecolor','k');
    p_moon_pred       = plot(log10(bodies.moon.MR),      bodies.moon.LyapHaloTp_pred,      '+', 'markersize', 18, 'markeredgecolor','k');
    
    PlotBoi2('$log_{10}(\mu)$','$T_{P_{Lyap/Halo}}$',26, 'LaTex')
    legend([p_enceleadus, p_europa, p_ganymede, p_triton, p_titan, p_moon, p_enceleadus_pred], 'Enceladus', 'Europa', 'Ganymede', 'Triton', 'Titan', 'Moon', 'Prediction','location', 'northwest','FontSize',18)
    
    LyapHaloTps = [bodies.enceladus.LyapHaloTp, bodies.europa.LyapHaloTp, bodies.ganymede.LyapHaloTp, bodies.triton.LyapHaloTp, bodies.titan.LyapHaloTp, bodies.moon.LyapHaloTp];
    errors_LyapHaloTps = [bodies.enceladus.LyapHaloTp_pred - bodies.enceladus.LyapHaloTp;...
                          bodies.europa.LyapHaloTp_pred - bodies.europa.LyapHaloTp;...
                          bodies.ganymede.LyapHaloTp_pred - bodies.ganymede.LyapHaloTp;...
                          bodies.triton.LyapHaloTp_pred - bodies.triton.LyapHaloTp;...
                          bodies.titan.LyapHaloTp_pred - bodies.titan.LyapHaloTp;...
                          bodies.moon.LyapHaloTp_pred - bodies.moon.LyapHaloTp];
    title(sprintf('Ref. Body: %s ... Mean/Max Error: %1.1e / %1.1e', LyapHaloTp_referenceBody.title, mean(abs(errors_LyapHaloTps)), max(abs(errors_LyapHaloTps))))
%     f = polyval(p,log10(MRs));
%     T = table(log10(MRs),LyapHaloTps,f,LyapHaloTps-f,'VariableNames',{'log10(MR)','Tp_Lyap','Fit','FitError'});
%     title(sprintf('Mean/Max Error: %1.1e / %1.1e', mean(abs(T.FitError)), max(abs(T.FitError))))



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







