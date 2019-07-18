% ========================================================================
%%% Description
% ========================================================================
% Function to find and analyze Lyapunov periodic orbits in various 3-body
% systems

% Created: 06/20/19
% Author : Luke Bury, luke.bury@colorado.edu
% ========================================================================
%%% Initialization
% ========================================================================
clear
clc
close all
mbinPath    = '~/CU_Google_Drive/Documents/MatGit/mbin';
PO_SavePath = '~/CU_Google_Drive/Documents/MatGit/CU/Research/Moon_Landing/Manifold_Study/PO_X0s';
addpath(genpath(mbinPath))
ticWhole = tic;

% ========================================================================
%%% Importing Data
% ========================================================================
%%% General data on solar system bodies
bodies = getBodyData(mbinPath);

%%% Color options/schemes
colors = get_colors();

% ========================================================================
%%% Run Switches
% ========================================================================
limitEnergy    = 0; % Set upper bound on energy of PO family (chosen below)

print_PO_index = 1; % Print PO index as family is generated?
plot_POs       = 1; % Plot the generate POs?

plot_PO_energy = 1; % Plot increasing energy of computed family

plot_stability = 1; % Plot stability index of the family

save_PO_X0s    = 0; % Save the initial conditions?

% ========================================================================
%%% Setup
% ========================================================================
% -------------------------------------------------
%%% System
% -------------------------------------------------
%%% Choose bodies
% primary = bodies.earth;     secondary = bodies.moon;
% primary = bodies.jupiter;   secondary = bodies.europa;
% primary = bodies.jupiter;   secondary = bodies.ganymede;
% primary = bodies.jupiter;   secondary = bodies.callisto;
% primary = bodies.saturn;    secondary = bodies.enceladus;
primary = bodies.saturn;    secondary = bodies.titan;
% primary = bodies.neptune;   secondary = bodies.triton;

%%% Normalizing constants
rNorm = secondary.a;         % n <-> km
tNorm = 1/secondary.meanMot; % n <-> sec
vNorm = rNorm / tNorm;       % n <-> km/sec

%%% Shortcut variables
mu   = secondary.MR;
R2_n = secondary.R_n;

%%% Equillibrium Points
rLPs_n = EquilibriumPoints(mu);
% ------------------------------------------------- 
%%% Periodic Orbit options
% -------------------------------------------------
n_POs_max    = 20;
LP           = 2;    % Lagrange Point of interest
PO_plot_skip = 1;

%%% Secondary-specific
if isequal(secondary.name,'europa')
    if limitEnergy == 1
        %%% Upper bound on energy (L2 flyover velocity, m/s)
        ub_L2FlyoverVelocity_mps = 100; % 'NaN' if no upper bound is desired
    end 
    
    %%% Step size for finding PO family
    familyStepSize = .0001; % dS_PO
    
elseif isequal(secondary.name,'enceladus')
    if limitEnergy == 1
        %%% Upper bound on energy (L2 flyover velocity, m/s)
        ub_L2FlyoverVelocity_mps = 50; % 'NaN' if no upper bound is desired
    end 
    
    %%% Step size for finding PO family
    familyStepSize = .00001; % dS_PO
    
elseif isequal(secondary.name,'titan')
    if limitEnergy == 1
        %%% Upper bound on energy (L2 flyover velocity, m/s)
        ub_L2FlyoverVelocity_mps = 1000; % 'NaN' if no upper bound is desired
    end 
    
    %%% Step size for finding PO family
    familyStepSize = .0001; % dS_PO
    
elseif isequal(secondary.name,'triton')
    if limitEnergy == 1
        %%% Upper bound on energy (L2 flyover velocity, m/s)
        ub_L2FlyoverVelocity_mps = 1000; % 'NaN' if no upper bound is desired
    end 
    
    %%% Step size for finding PO family
    familyStepSize = .0001; % dS_PO

elseif isequal(secondary.name,'moon')
    %%% Step size for finding PO family
    familyStepSize = 0.0001; % dS_PO
end

% ------------------------------------------------- 
%%% Integration Options
% -------------------------------------------------
tol = 1e-13;
options = odeset('RelTol',tol,'AbsTol',tol);
prms.u    = mu;
prms.R2_n = R2_n;

% ========================================================================
%%% Evaluate the Dynamics at Equillibrium point
% ========================================================================
%%% Grab location of equillibrium point
rLP_n = rLPs_n(LP,:);

%%% Evaluate jacobian of EOM and state at equillibrium point
[A_LP] = get_Amat_CR3BP(mu, rLP_n);

% ------------------------------------------------- 
%%% Eigenvectors and Eigenvalues
% -------------------------------------------------
%%% Getting eigenvectors and eigenvalues
[eVecs, eVals] = eig(A_LP);

%%% Choose appropriate eigenvectors to perturb in direciton of
ev_i = 5;
eVal = eVals(ev_i,ev_i);
eVec = eVecs(:,ev_i);

% %%% Finding center manifold (purely imaginary eigenvalue) ... or maybe not
% for ev_i = 1:size(eVals,1)
%     if abs(real(eVals(ev_i,ev_i))) < 5e-16 && abs(imag(eVals(ev_i,ev_i))) > 5e-16 % If real() == 0 and imag() ~= 0
%         ev_i = 5;
%         eVal_CM = eVals(ev_i,ev_i);
%         eVec_CM = eVecs(:,ev_i);
% %         eVal_CM = eVals(3,3);
% %         eVec_CM = eVecs(:,3);
% 
%         break
%     end
% end

% ========================================================================
%%% Calculate Family of Periodic Orbits
% ========================================================================
% ------------------------------------------------- 
%%% Family-finding options
% -------------------------------------------------
%%% Initial perturbation from equilibrium point
initialPerturbationScale = 1e-4;
% initialPerturbationScale = 1e-6;


%%% Tangent orbit scaling parameter
tanScale = 1e-5; % Tuning parameter
% tanScale = 1e-8; % Tuning parameter

%%% Error tolerance in finding PO
error_tol = 1e-10;


% ------------------------------------------------- 
%%% Perturbing the Lagrange-point state in the direction of the center manifold
% -------------------------------------------------
%%% Perturbing the equilibrium point to find first periodic orbit
X_pert = [rLP_n'; 0; 0; 0] + real(eVec).*initialPerturbationScale;

%%% Estimating period of PO
T_pert = 2*pi/imag(eVal);

%%% First Tangent Orbit
X_tan = real(eVec)/norm(real(eVec));
T_tan = 0;

%%% Stepping in direction of first tangent orbit
X0_guess_n = X_pert + tanScale*X_tan;
T_guess_n = T_pert + tanScale*T_tan;

%%% Preallocating
POs = NaN(n_POs_max,7);
stabilityIndices = NaN(n_POs_max,2);
L2FlyoverVelocities_mps = NaN(n_POs_max,1);

%%% Storing first PO
POs(1,:) = [X0_guess_n', T_guess_n];

%%% Storing energy of first PO
JC_PO_i = JacobiConstantCalculator(mu,X0_guess_n(1:3)',X0_guess_n(4:6)');
L2_FlyoverVelocity_PO_i_mps = JC_2_L2FlyoverVelocity(JC_PO_i,mu,rLPs_n(LP,:),vNorm);
L2FlyoverVelocities_mps(1) = L2_FlyoverVelocity_PO_i_mps;

% ------------------------------------------------- 
%%% Calculating family of POs
% -------------------------------------------------
%%% Initializating
for PO_i = 1:n_POs_max
    if print_PO_index == 1
        fprintf('PO_i = %1d\n',PO_i)
    end
    
    %%% Resetting error and PO counter
    error = 1;
% %     counter = 0;

    %%% Finding new periodic orbit
    while error > error_tol
        % --------------------------
        % Preparing for integration
        % --------------------------
        %%% Initialize STM
        stm0 = eye(6);
        
        %%% Set initial conditions
        X0_n = [X0_guess_n; reshape(stm0,36,1)];
        time0_n = [0,T_guess_n];
        
        %%% Integrating
        [T_n, X_n] = ode113(@Int_CR3BnSTM, time0_n, X0_n, options, prms);
        
        %%% Retreiving updated STM
        stm_tf_t0 = reshape(X_n(end,7:42),6,6);

        %%% Correcting the Initial Guess
        F = [X_n(end,1:6)-X_n(1,1:6), (X_n(1, 1:6)' - X_pert)'*Int_CR3Bn(0,X_pert,prms), (X_n(1, 1:6)' - X_pert)'*X_tan + (T_guess_n - T_pert)*T_tan - familyStepSize]';
        D = [stm_tf_t0-eye(6), Int_CR3Bn(0,X0_guess_n,prms);...
             Int_CR3Bn(0,X_pert,prms)', 0;...
             X_tan', T_tan];
        delta = (D'*D)\D'*(-F);
        
        X0_guess_n = X0_guess_n + delta(1:6);
        T_guess_n = T_guess_n + delta(7);
        
        %%% Storing error
        error = norm(F);
        
% %         %%% Counting error tries
% %         counter = counter + 1;
    end
    
    %%% Calculate energy
    JC_PO_i = JacobiConstantCalculator(mu,X0_guess_n(1:3)',X0_guess_n(4:6)');
    L2_FlyoverVelocity_PO_i_mps = JC_2_L2FlyoverVelocity(JC_PO_i,mu,rLPs_n(LP,:),vNorm);
    
    if limitEnergy == 1    
        if L2_FlyoverVelocity_PO_i_mps > ub_L2FlyoverVelocity_mps
            break
        end
    end
    
    %%% Store energy
    L2FlyoverVelocities_mps(PO_i+1) = L2_FlyoverVelocity_PO_i_mps;
    
    %%% Check stability indices
    if plot_stability == 1
        [eVecs_PO, eVals_PO] = eig(stm_tf_t0);
        [S1, S2] = getStabilityIndices(diag(eVals_PO));
        
        stabilityIndices(PO_i+1,:) = [S1, S2];
    end
    
    %%% Storing new periodic orbit
    POs(PO_i+1,:) = [X0_guess_n', T_guess_n];

    %%% Update perturbed orbit
    X_pert = X0_guess_n;
    T_pert = T_guess_n;

    %%% Compute new tangent orbit
    tangentOrbit_hat = ((POs(PO_i+1,:) - POs(PO_i,:)) / norm(POs(PO_i+1,:)-POs(PO_i,:)))';

    %%% Update the tangent orbit
    X_tan = tangentOrbit_hat(1:6);
    T_tan = tangentOrbit_hat(7);

    %%% Update guess for first tangent orbit
    X0_guess_n = X_pert + familyStepSize*X_tan;
    T_guess_n  = T_pert + familyStepSize*T_tan;

    %%% Should this orbit be plotted?
    if plot_POs == 1
        %%% Plotting some fraction of all the orbits
        if mod(PO_i,PO_plot_skip) == 0
            figure(1); hold all
            p2 = plot3(X_n(:,1),X_n(:,2),X_n(:,3),'k');

        end
    end % plot_POs
    
end

% ------------------------------------------------- 
%%% Adding detail to plot
% -------------------------------------------------
if plot_POs == 1
    figure(1); hold all; grid on;
    plot3(rLPs_n(LP,1),rLPs_n(LP,2),rLPs_n(LP,3),'^','markeredgecolor',colors.std.black,'markerfacecolor','k')
    PlotBoi3('$X_n$','$Y_n$','$Z_n$',20,'LaTex')
    view(0,90)
    axis equal
end % if plot_POs == 1

% ------------------------------------------------- 
%%% Plotting family energy
% -------------------------------------------------
if plot_PO_energy == 1
    figure; hold all
    plot(L2FlyoverVelocities_mps,'o','markersize',6,'markeredgecolor',colors.std.blue,'markerfacecolor',colors.std.ltblue)
    PlotBoi2('PO Index','$L_2$ Flyover Velocity, $m/s$',18,'LaTex')
end % plot_PO_energy

% ------------------------------------------------- 
%%% Plotting stability indices
% -------------------------------------------------
if plot_stability == 1
    figure
    subplot(1,2,1); hold all
    p1 = plot(stabilityIndices(:,1),'o','markersize',6,'markeredgecolor',colors.std.purp,'markerfacecolor',colors.std.ltpurp);
    p2 = plot(stabilityIndices(:,2),'o','markersize',6,'markeredgecolor',colors.std.red,'markerfacecolor',colors.std.ltred);
    PlotBoi2('PO Index','Stability Index',18,'LaTex')
    xlim([0 size(stabilityIndices,1)])
    legend([p1 p2],'S1','S2')
    
    subplot(1,2,2); hold all
    p1 = plot(stabilityIndices(:,1),'o','markersize',6,'markeredgecolor',colors.std.purp,'markerfacecolor',colors.std.ltpurp);
    p2 = plot(stabilityIndices(:,2),'o','markersize',6,'markeredgecolor',colors.std.red,'markerfacecolor',colors.std.ltred);
    plot([0 size(stabilityIndices,1)],[2 2],'k','linewidth',1)
    PlotBoi2('PO Index','',18,'LaTex')
    legend([p1 p2],'S1','S2')
    ylim([1.9 2.1])
    xlim([0 size(stabilityIndices,1)])
    title('Zoomed')
end

% ========================================================================
%%% Closeout
% ========================================================================
toc(ticWhole);















