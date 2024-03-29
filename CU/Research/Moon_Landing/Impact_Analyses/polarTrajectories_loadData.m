clear
clc
close all
mbinPath = '~/CU_Google_Drive/Documents/MatGit/mbin';
moonFuncsPath = '~/CU_Google_Drive/Documents/MatGit/CU/Research/Moon_Landing/Moon_Landing_funcs';
polarTrajDataPath = '~/CU_Google_Drive/Documents/MatGit/CU/Research/Moon_Landing/Impact_Analyses/polarTrajectoryConditions';
addpath(genpath(mbinPath))
addpath(genpath(moonFuncsPath))
tic

% ========================================================================
%%% Run Switches and files
% ========================================================================
run_NominalTrajectories = 0; % if you want to integrate and plot all the cases
useEvents = 1; % Stop trajectories upon impact or escape

% ========================================================================
%%% Free variables
% ========================================================================
%%% JC scaling
% JC_scalars = linspace(1e-8, 2e-6, 10) + 1;
% JC_scalars = linspace(1e-8, 2e-5, 10) + 1;

% %%% For Europa L2 350-300-250-200-150-100-50 mps flyover energies
% JC_scalars = [3.002960593980811,3.003132670751979,3.003278274173736,3.003397404246083,...
%     3.003490060969019,3.003556244342545,3.003595954366661] ./ 3.002960593980811;

%%% For Europa L2 350-300 mps flyover energies
JC_scalars = linspace(1, 3.003132670751979/3.002960593980811,25);

%%% Load file
loadFile = '/PolarXs_eur_350mps_1000km_79v0s.txt';
polarTrajFile = [polarTrajDataPath,loadFile];

%%% Timing info
n_dt = 10000;
t_i = 0;
t_f = 4*pi;
time_fwd = linspace(t_i,t_f,n_dt);
time_bkd = linspace(t_f,t_i,n_dt);

% ========================================================================
%%% Importing Data
% ========================================================================
%%% General data on solar system bodies
bodies = getBodyData(mbinPath);

%%% Color options/schemes
colors = get_colors();

% ========================================================================
%%% Setup
% ========================================================================
% -------------------------------------------------
% Choose system
% -------------------------------------------------
%%% 3B system
primary   = bodies.jupiter;
secondary = bodies.europa;

%%% Normalizing constants
rNorm = secondary.a;         % n <-> km
tNorm = 1/secondary.meanMot; % n <-> sec
vNorm = rNorm / tNorm;       % n <-> km/sec

%%% Lagrange Points
L123 = EquilibriumPoints(secondary.MR,1:3); % [3x3] of L1, L2, L3 normalized BCR coordinates

% -------------------------------------------------
% Load and parse data
% -------------------------------------------------
%%% Load polar impact-traj data
polarImpactMat = dlmread(polarTrajFile,',',1,0);

%%% Grab indices of tX0s and tXfs ([time, state])
tX0s_ind = find(polarImpactMat(:,1) == 0);
tXfs_ind = find(polarImpactMat(:,1) ~= 0);

%%% Store tX0s and tXfs
tX0s = polarImpactMat(tX0s_ind,:);
tXfs = polarImpactMat(tXfs_ind,:);

%%% Find Jacobi Constant of nominal cases
JC_scInitial = JacobiConstantCalculator(secondary.MR,tX0s(1,2:4) ,tX0s(1,5:7));

%%% Find scaled JCs and delta-JCs
JCs_scaled = JC_scalars .* JC_scInitial;

dJCs = JCs_scaled - JC_scInitial;
% -------------------------------------------------
% Integration settings
% -------------------------------------------------
%%% Choosing ode tolerance
tol = 1e-13;

%%% Setting integrator options
options = odeset('RelTol',tol,'AbsTol',tol);
options_ImpactEscape = odeset('Events',@event_ImpactEscape_CR3Bn,'RelTol',tol,'AbsTol',tol);

%%% Parameters for integrator
prms = struct();
prms.u = secondary.MR;
prms.R2_n = secondary.R_n;
prms.L1x = L123(1,1);
prms.L2x = L123(2,1);
prms.R1_n = primary.R/rNorm;
prms.J21 = primary.J2;

% -------------------------------------------------
% Plotting colors
% -------------------------------------------------
[ colorMatrix ] = colorScale([colors.std.black; colors.std.ltgrn],length(JC_scalars) );
% ========================================================================
%%% Integrating nominal cases
% ========================================================================
if run_NominalTrajectories == 1
    %%% Setting up plot
    figure; hold all
    plotBodyTexture3(secondary.R_n,[1-secondary.MR,0,0],secondary.img)
    axis equal
    PlotBoi3('x','y','z',16,'LaTex')

    %%% Integrating X0s forward
    for kk = 1:size(tX0s,1)
        %%% Propagating X0 forward in time
        [time_X0_n, X_BCR_X0_n, time_X0_eventImpact, X_X0_eventImpact, index_X0_eventImpact] = ode113(@Int_CR3Bn,...
            time_fwd, tX0s(kk,2:7), options_ImpactEscape, prms);

        %%% Plotting trajectory
        plot3(X_BCR_X0_n(:,1),X_BCR_X0_n(:,2),X_BCR_X0_n(:,3),'b','linewidth',1.5)
    end

    %%% Integrating Xfs backward
    for kk = 1:size(tXfs,1)
        %%% Propagating Xf backward in time
        [time_Xf_n, X_BCR_Xf_n, time_Xf_eventImpact, X_Xf_eventImpact, index_Xf_eventImpact] = ode113(@Int_CR3Bn,...
            time_bkd, tXfs(kk,2:7), options_ImpactEscape, prms);

        %%% Plotting trajectory
        plot3(X_BCR_Xf_n(:,1),X_BCR_Xf_n(:,2),X_BCR_Xf_n(:,3),'r','linewidth',1.5)
    end
end % run_NominalTrajectories

% ========================================================================
%%% Integrating reduced velocity cases
% ========================================================================


for traj_id = 3
    % for traj_k = 1:size(tX0s,1)
    figure('position',[23 71 1382 719])

    %%% Integrating nominal case
    if useEvents == 0 
        [time_nom_n, X_BCR_nom_n] = ode113(@Int_CR3Bn,time_fwd, tX0s(traj_id,2:7), options, prms);
    elseif useEvents == 1
        [time_nom_n, X_BCR_nom_n, time_nom_eventImpact, X_nom_eventImpact, index_nom_eventImpact] = ode113(@Int_CR3Bn,...
            time_fwd, tX0s(traj_id,2:7), options_ImpactEscape, prms);
    end

    %%% Plotting nominal case
    subplot(1,2,1); hold all
    % plot3(X_BCR_nom_n(:,1),X_BCR_nom_n(:,2),X_BCR_nom_n(:,3),'k','linewidth',1.5)
    plotBodyTexture3(secondary.R_n,[1-secondary.MR,0,0],secondary.img)
    axis equal
    PlotBoi3('x','y','z',16,'LaTex')
    camva(8)

    subplot(1,2,2); hold all
    % plot3(X_BCR_nom_n(:,1),X_BCR_nom_n(:,2),X_BCR_nom_n(:,3),'k','linewidth',1.5)
    plotBodyTexture3(secondary.R_n,[1-secondary.MR,0,0],secondary.img)
    axis equal
    PlotBoi3('x','y','z',16,'LaTex')
    camva(8)

    %%% Integrating and plotting perturbed cases
    for dJC_i = 1:length(JC_scalars)
        % -------------------------------------------------
        % X0 forward in time
        % -------------------------------------------------
        %%% Velocity magnitude and direction of current trajectory
        vMag_X0_old = norm(tX0s(traj_id,5:7));
        vHat_X0 = tX0s(traj_id,5:7)./vMag_X0_old;

        %%% New velocity
        vMag_X0_new = sqrt(vMag_X0_old^2 - dJCs(dJC_i));
        vNew_X0 = vHat_X0 .* vMag_X0_new;
        
        %%% If it's still possible to start at this location in the neck
        if imag(vMag_X0_new) == 0

            %%% Propagating X0 forward in time
            if useEvents == 0
                [time_X0_n, X_BCR_X0_n] = ode113(@Int_CR3Bn, time_fwd, [tX0s(traj_id,2:4),vNew_X0], options, prms);
            elseif useEvents == 1
                [time_X0_n, X_BCR_X0_n, time_X0_eventImpact, X_X0_eventImpact, index_X0_eventImpact] = ode113(@Int_CR3Bn,...
                    time_fwd, [tX0s(traj_id,2:4),vNew_X0], options_ImpactEscape, prms);
            end

            %%% Plotting trajectory
            subplot(1,2,1)
    %         [p_X0] = plot3(X_BCR_X0_n(:,1),X_BCR_X0_n(:,2),X_BCR_X0_n(:,3),'b','linewidth',1.5);
            [p_X0] = plot3(X_BCR_X0_n(:,1),X_BCR_X0_n(:,2),X_BCR_X0_n(:,3),'linewidth',1.5,'color',colorMatrix(dJC_i,:));

            plotCR3BP_YZNeck( JC_scInitial+dJCs(dJC_i), secondary.MR , 2, 0, prms, colorMatrix(dJC_i,:),2)

            drawnow
        end % imag(vMag_X0_new) == 0
        % -------------------------------------------------
        % Xf backward in time
        % -------------------------------------------------
        %%% Velocity magnitude and direction of current trajectory
        vMag_Xf_old = norm(tXfs(traj_id,5:7));
        vHat_Xf = tXfs(traj_id,5:7)./vMag_Xf_old;

        %%% New velocity
        vMag_Xf_new = sqrt(vMag_Xf_old^2 - dJCs(dJC_i));
        vNew_Xf = vHat_Xf .* vMag_Xf_new;

        %%% Propagating Xf backward in time
        if useEvents == 0
            [time_Xf_n, X_BCR_Xf_n] = ode113(@Int_CR3Bn, time_bkd, [tXfs(traj_id,2:4),vNew_Xf], options, prms);
        elseif useEvents == 1
            [time_Xf_n, X_BCR_Xf_n, time_Xf_eventImpact, X_Xf_eventImpact, index_Xf_eventImpact] = ode113(@Int_CR3Bn,...
                time_bkd, [tXfs(traj_id,2:4),vNew_Xf], options_ImpactEscape, prms);
        end

        %%% Plotting trajectory
        subplot(1,2,2)
%         [p_Xf] = plot3(X_BCR_Xf_n(:,1),X_BCR_Xf_n(:,2),X_BCR_Xf_n(:,3),'r','linewidth',1.5);
        [p_Xf] = plot3(X_BCR_Xf_n(:,1),X_BCR_Xf_n(:,2),X_BCR_Xf_n(:,3),'linewidth',1.5,'color',colorMatrix(dJC_i,:));
        
        plotCR3BP_YZNeck( JC_scInitial+dJCs(dJC_i), secondary.MR , 2, 0, prms, colorMatrix(dJC_i,:),2)
        
        drawnow

        %%% Printing info
        fprintf('----------------------\n')
        fprintf('traj_id = %1.0d\t%%JC Change = %1.7f%%\n',traj_id,(JC_scalars(dJC_i)-1)*100)
        fprintf('X_fwd "final" miss distance    = %1.4f km\n',norm(X_BCR_X0_n(end,1:3)-X_BCR_Xf_n(1,1:3))*rNorm)
        fprintf('X_bkd "starting" miss distance = %1.4f km\n',norm(X_BCR_Xf_n(end,1:3)-X_BCR_X0_n(1,1:3))*rNorm)
        fprintf('Total X0 velocity reduction = %1.3f mps\n',(vMag_X0_new - vMag_X0_old)*vNorm*1000)
        fprintf('Total Xf velocity reduction = %1.3f mps\n',(vMag_Xf_new - vMag_Xf_old)*vNorm*1000)

    end % dJC_i = 1:length(JC_scalars)

    %%% Extra plot settings
    % p_end = plot3(X_BCR_Xf_n(1,1),X_BCR_Xf_n(1,2),X_BCR_Xf_n(1,3),'r.','markersize',15);
    % p_start = plot3(X_BCR_X0_n(1,1),X_BCR_X0_n(1,2),X_BCR_X0_n(1,3),'b.','markersize',15);
    % plotBodyTexture3(secondary.R_n,[1-secondary.MR,0,0],secondary.img)
    % legend([p_X0 p_Xf],'X0 Forward','Xf Backward','Nominal Start','Nominal End','location','northeast')
    % title(sprintf('JC Reduced %1.7f%%',(JC_scalars(dJC_i)-1)*100))
    % axis equal
    % PlotBoi3('x','y','z',16,'LaTex')
    subplot(1,2,1)
    [p_nom] = plot3(X_BCR_nom_n(:,1),X_BCR_nom_n(:,2),X_BCR_nom_n(:,3),'k','linewidth',3);
    title('Correct X0, Propagated forward with decreasing JCs','FontSize',20)
    h = legend([p_nom], 'Nominal Trajectory','location','northeast');
    h.FontSize = 15;
    subplot(1,2,2); hold all
    [p_nom] = plot3(X_BCR_nom_n(:,1),X_BCR_nom_n(:,2),X_BCR_nom_n(:,3),'k','linewidth',3);
    title('Correct Xf, Propagated backward with decreasing JCs','FontSize',20)
    h = legend([p_nom], 'Nominal Trajectory','location','northeast');
    h.FontSize = 15;
end % kk = 1:size(tX0s,1)







toc






