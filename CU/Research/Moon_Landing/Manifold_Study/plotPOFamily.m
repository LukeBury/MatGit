% ========================================================================
%%% Description
% ========================================================================
% Plot families of pre-computed periodic orbits

% Created: 09/30/19
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
plot_primary   = 0;
plot_secondary = 1;
plot_L1        = 0;
plot_L2        = 0;

plot_onlyOneOrbit    = 1;
    chosenOrbitIndex = 1;
plot_stabilityIndices = 1;

% ========================================================================
%%% Setup
% ========================================================================
% -------------------------------------------------
%%% Choose data
% -------------------------------------------------
% family = 'Earth_Moon.CR3BP.L1_Lyapunov.txt';
% family = 'Earth_Moon.CR3BP.L1_Vertical.txt';
% family = 'Earth_Moon.CR3BP.L1_SHalo.txt';
% family = 'Earth_Moon.CR3BP.L2_Lyapunov.txt';
% family = 'Earth_Moon.CR3BP.L2_SHalo.txt';

% family = 'Jupiter_Europa.CR3BP.L2_Lyapunov.txt';
% family = 'Jupiter_Europa.CR3BP.L2_NHalo.txt';
% family = 'Jupiter_Europa.CR3BP.L2_SHalo.txt';
% family = 'Jupiter_Europa.CR3BP.L2_Vertical.txt';
% family = 'Jupiter_Europa.CR3BP.L2_EasternAxial.txt';
% family = 'Jupiter_Europa.CR3BP.L2_WesternAxial.txt';
% family = 'Jupiter_Europa.CR3BP_J2pJ4pJ6pJ2s.L2_Lyapunov.txt';
% family = 'Jupiter_Europa.CR3BP_J2pJ4pJ6pJ2s.L2_SHalo.txt';
% family = 'Jupiter_Europa.CR3BP_J2pJ4pJ6pJ2s.L2_Vertical.txt';

% family = 'Jupiter_Ganymede.CR3BP.L2_Lyapunov.txt';
% family = 'Jupiter_Ganymede.CR3BP.L2_Vertical.txt';
% family = 'Jupiter_Ganymede.CR3BP.L2_SHalo.txt';

% family = 'Saturn_Enceladus.CR3BP.L1_Lyapunov.txt';
% family = 'Saturn_Enceladus.CR3BP.L1_NHalo.txt';
% family = 'Saturn_Enceladus.CR3BP.L1_SHalo.txt';
% family = 'Saturn_Enceladus.CR3BP.L2_EasternAxial.txt';
% family = 'Saturn_Enceladus.CR3BP.L2_Lyapunov.txt';
% family = 'Saturn_Enceladus.CR3BP.L2_NHalo.txt';
% family = 'Saturn_Enceladus.CR3BP.L2_SHalo.txt';
% family = 'Saturn_Enceladus.CR3BP.L2_Vertical.txt';
% family = 'Saturn_Enceladus.CR3BP.L2_WesternAxial.txt';
% family = 'Saturn_Enceladus.CR3BP.SLeadingSneakerHeel.txt';
% family = 'Saturn_Enceladus.CR3BP.SLeadingSneakerToe.txt';
% family = 'Saturn_Enceladus.CR3BP.unlabeled1LeadingS.txt';
% family = 'Saturn_Enceladus.CR3BP_J2pJ4pJ6pJ2s.L2_Lyapunov.txt';
% family = 'Saturn_Enceladus.CR3BP_J2pJ4pJ6pJ2s.L2_SHalo.txt';
% family = 'Saturn_Enceladus.CR3BP_J2pJ4pJ6pJ2s.L2_Vertical.txt';

% family = 'Saturn_Titan.CR3BP.L2_Lyapunov.txt';
% family = 'Saturn_Titan.CR3BP.L2_Vertical.txt';

% family = 'Neptune_Triton.CR3BP.L2_Lyapunov.txt';
% family = 'Neptune_Triton.CR3BP.L2_Vertical.txt';
family = 'Neptune_Triton.CR3BP.L2_SHalo.txt';

%%% Path from mbin to data
dataPathFromMBin = '/Data/InitialConditions/PO_Families/';

%%% PO data file
PO_datafile = [mbinPath, dataPathFromMBin, family];

% -------------------------------------------------
%%% Plot Options
% -------------------------------------------------
lw = 2; % linewidth

% -------------------------------------------------
%%% Set up the system
% -------------------------------------------------
% --------------------------
% Set primary & secondary
% --------------------------
[primary, secondary] = assignPrimaryAndSecondary_CR3BP(family, bodies);

%%% Normalizing constants
rNorm = secondary.a;         % n <-> km
tNorm = 1/secondary.meanMot; % n <-> sec
vNorm = rNorm / tNorm;       % n <-> km/sec

%%% prms for integration
prms.u    = secondary.MR;
prms.R2_n = secondary.R_n;
prms.n    = 1;



%%% Equillibrium Points
if contains(family,'.CR3BP.')
    rLPs_n = EquilibriumPoints(prms.u, 1);
elseif contains(family,'.CR3BP_J2pJ4pJ6pJ2s.')
    rLPs_n = collinearEquilibriumPoints_ZH(prms);
    prms.J2p = primary.J2; 
    prms.J4p = primary.J4; 
    prms.J6p = primary.J6; 
    prms.J2s = secondary.J2;
    prms.R2 = secondary.R_n;
    prms.R1 = primary.R / rNorm;
    
    tNorm = sqrt((rNorm^3)/(bodies.constants.G*(primary.mass + secondary.mass)));
    vNorm = rNorm/tNorm;
    prms.n = secondary.meanMot * tNorm;

end


% -------------------------------------------------
%%% Integration options
% -------------------------------------------------
tol = 1e-13;
options = odeset('RelTol',tol,'AbsTol',tol);

% -------------------------------------------------
%%% Load data
% -------------------------------------------------
%%% Load the data file
PO_data = dlmread(PO_datafile,',',1,0);

%%% Number of ICs
n_POs = size(PO_data,1);

%%% Set column specifiers
c_X0   = 1:6;
c_Tp   = 7;
c_JC   = 8;
c_L2EV = 9;

% ========================================================================
%%% Integrate
% ========================================================================
% -------------------------------------------------
%%% Preallocate
% -------------------------------------------------
%%% Preallocate cell array to hold trajectories
trajs = cell(n_POs,1);

%%% Create figure
figure(837920); hold all
% 989
% figure(123456); hold all

% -------------------------------------------------
%%% Loop through conditions and integrate
% -------------------------------------------------
if plot_onlyOneOrbit == 1
    orbitRange = chosenOrbitIndex;
else
    orbitRange = 1:n_POs;
%     989
%     orbitRange = [6, 35, 50, 80, 95, 105, 112, 120, 138];
end

%%% Preallocate
if plot_stabilityIndices
    stabilityIndices = NaN(length(orbitRange),2);
    stm0_vec = reshape(eye(6), 36, 1);
end

for kk = orbitRange
    %%% Get current X0 and T
    X0n_i = PO_data(kk,c_X0)';
    Tpn_i = PO_data(kk,c_Tp);
    
    %%% Integrate
    if contains(family,'.CR3BP.')
        if plot_stabilityIndices == 1
            [~, Xn_i] = ode113(@Int_CR3BnSTM, [0 Tpn_i], [X0n_i; stm0_vec], options, prms);
        else
            [~, Xn_i] = ode113(@Int_CR3Bn, [0 Tpn_i], X0n_i, options, prms);
        end
    elseif contains(family,'.CR3BP_J2pJ4pJ6pJ2s.')
        if plot_stabilityIndices  == 1
            [~, Xn_i] = ode113(@Int_CR3BnSTM_J2pJ4pJ6pJ2s, [0 Tpn_i], [X0n_i; stm0_vec], options, prms);
        else
            [~, Xn_i] = ode113(@Int_CR3Bn_J2pJ4pJ6pJ2s, [0 Tpn_i], X0n_i, options, prms);
        end
    end
    
    if plot_stabilityIndices
        stm_tf_t0                           = reshape(Xn_i(end,7:42),6,6);
        monodromy                           = stm_tf_t0;
        [eigenVectors_new, eigenValues_new] = eig(monodromy);
        [S1, S2]                            = getStabilityIndices(diag(eigenValues_new));
        stabilityIndices(kk,:)              = [S1, S2];
    end
    
    
    %%% Plot
    if contains(family,'.CR3BP.')
        plot3(Xn_i(:,1),Xn_i(:,2),Xn_i(:,3),'b','linewidth',lw)
    elseif contains(family,'.CR3BP_J2pJ4pJ6pJ2s.')
        plot3(Xn_i(:,1),Xn_i(:,2),Xn_i(:,3),'r','linewidth',lw)
    end
    
    
    %%% Store trajectory
    trajs{kk}.Xn = Xn_i(:,1:6);
    
end

% -------------------------------------------------
%%% Format the plot and add detail
% -------------------------------------------------
% PlotBoi3('$x_n$','$y_n$','$z_n$',20,'LaTex')
PlotBoi3_CR3Bn(20)
if plot_primary  == 1
    plotPrimary(primary,secondary);
end
if plot_secondary == 1
    plotSecondary(secondary);
end
axis equal

if plot_L1 == 1
    plot3(rLPs_n(1,1), rLPs_n(1,2), rLPs_n(1,3),'^','markerfacecolor',colors.grn,'markeredgecolor',colors.black,'markersize',8)
end
if plot_L2 == 1
    plot3(rLPs_n(2,1), rLPs_n(2,2), rLPs_n(2,3),'^','markerfacecolor',colors.ltgrn,'markeredgecolor',colors.black,'markersize',8)
end

% -------------------------------------------------
%%% Plot stability indices
% -------------------------------------------------
if plot_stabilityIndices == 1
    s1_ind = ~isnan(stabilityIndices(:,1));
    s2_ind = ~isnan(stabilityIndices(:,2));
    S1 = stabilityIndices(s1_ind,1);
    S2 = stabilityIndices(s2_ind,2);
    JC_temp1 = PO_data(s1_ind,c_JC);
    JC_temp2 = PO_data(s2_ind,c_JC);
    
    figure('position',[209 322 948 302])
    subplot(1,2,1); hold all
%     p1 = plot(abs(S1),'o','markeredgecolor',colors.blue,'markerfacecolor',colors.ltblue);
%     p2 = plot(abs(S2),'o','markeredgecolor',colors.red,'markerfacecolor',colors.ltred);
%     plot([0 length(S1)],[2 2],'k','linewidth',1)
    p1 = plot(JC_temp1, abs(S1),'o','markeredgecolor',colors.blue,'markerfacecolor',colors.ltblue);
    p2 = plot(JC_temp2, abs(S2),'o','markeredgecolor',colors.red,'markerfacecolor',colors.ltred);
    plot(unique([min([JC_temp2, JC_temp1]) max([JC_temp2, JC_temp1])]),[2 2],'k','linewidth',1)
    PlotBoi2('Jacobi Constant','Stability Indices',18,'LaTex')
    legend([p1 p2],'S_1','S_2')
    xlim(unique([min([JC_temp2, JC_temp1]) max([JC_temp2, JC_temp1])]))
    
    subplot(1,2,2); hold all
%     p1 = plot(abs(S1),'o','markeredgecolor',colors.blue,'markerfacecolor',colors.ltblue);
%     p2 = plot(abs(S2),'o','markeredgecolor',colors.red,'markerfacecolor',colors.ltred);
%     plot([0 length(S1)],[2 2],'k','linewidth',1)
%     PlotBoi2('PO index','',18,'LaTex')
%     ylim([1.9 2.1])
%     xlim([0 length(S1)])
    p1 = plot(JC_temp1, abs(S1),'o','markeredgecolor',colors.blue,'markerfacecolor',colors.ltblue);
    p2 = plot(JC_temp2, abs(S2),'o','markeredgecolor',colors.red,'markerfacecolor',colors.ltred);
    plot(unique([min([JC_temp2, JC_temp1]) max([JC_temp2, JC_temp1])]),[2 2],'k','linewidth',1)
    PlotBoi2('Jacobi Constant','',18,'LaTex')
    ylim([1.9 2.1])
    xlim(unique([min([JC_temp2, JC_temp1]) max([JC_temp2, JC_temp1])]))

%     JC_temp1 = jacobiConstants(s1_ind);
%     JC_temp2 = jacobiConstants(s2_ind);
%     
%     figure('position',[209 322 948 302])
%     subplot(1,2,1); hold all
%     p1 = plot(JC_temp1, abs(S1),'o','markeredgecolor',colors.blue,'markerfacecolor',colors.ltblue);
%     p2 = plot(JC_temp2, abs(S2),'o','markeredgecolor',colors.red,'markerfacecolor',colors.ltred);
%     plot(unique([min([JC_temp2, JC_temp1]) max([JC_temp2, JC_temp1])]),[2 2],'k','linewidth',1)
%     PlotBoi2('Jacobi Constant','Stability Indices',18,'LaTex')
%     legend([p1 p2],'S_1','S_2')
%     xlim(unique([min([JC_temp2, JC_temp1]) max([JC_temp2, JC_temp1])]))
%     
%     subplot(1,2,2); hold all
%     p1 = plot(JC_temp1, abs(S1),'o','markeredgecolor',colors.blue,'markerfacecolor',colors.ltblue);
%     p2 = plot(JC_temp2, abs(S2),'o','markeredgecolor',colors.red,'markerfacecolor',colors.ltred);
%     plot(unique([min([JC_temp2, JC_temp1]) max([JC_temp2, JC_temp1])]),[2 2],'k','linewidth',1)
%     PlotBoi2('Jacobi Constant','',18,'LaTex')
%     ylim([1.9 2.1])
%     xlim(unique([min([JC_temp2, JC_temp1]) max([JC_temp2, JC_temp1])]))

end

% ========================================================================
%%% Closeout
% ========================================================================
tocWhole = toc(ticWhole);
fprintf('Elapsed time: %1.4f seconds\n',tocWhole)




%===== testing something
% As it's coded in the nCR3BP, norm(Xn_i(1,1:3) - Xn_i(end,1:3)) = 1.019935959809836e-11

































