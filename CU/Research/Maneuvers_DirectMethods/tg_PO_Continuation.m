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
print_PO_index  = 1;
plot_current_PO = 1;

plot_stability  = 1;
plot_energy     = 1;

save_PO_X0_database = 0;

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
primary = bodies.saturn;    secondary = bodies.enceladus;
% primary = bodies.saturn;    secondary = bodies.titan;
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
%%% Integration Options
% -------------------------------------------------
tol = 1e-13;
options = odeset('RelTol',tol,'AbsTol',tol);
options_XZStop = odeset('Event',@event_yEqualsZeroPastL2,'RelTol',tol,'AbsTol',tol);
prms.u    = mu;
prms.R2_n = R2_n;
t0 = 0;


% ------------------------------------------------- 
%%% PO Options
% -------------------------------------------------
n_POs_max = 1000;

% myPO_ICs = PO_ICs.JupiterEuropa.CR3BP.L2_Lyapunov; famName = 'JupiterEuropa.CR3BP.L2_Lyapunov';
% myPO_ICs = PO_ICs.JupiterEuropa.CR3BP.L2_Vertical; famName = 'JupiterEuropa.CR3BP.L2_Vertical';
% myPO_ICs = PO_ICs.JupiterEuropa.CR3BP.L2_SHalo; famName = 'JupiterEuropa.CR3BP.L2_SHalo';
% myPO_ICs = PO_ICs.JupiterEuropa.CR3BP.L2_EasternAxial; famName = 'JupiterEuropa.CR3BP.L2_EasternAxial';
% myPO_ICs = PO_ICs.JupiterEuropa.CR3BP.L2_WesternAxial; famName = 'JupiterEuropa.CR3BP.L2_WesternAxial';

% myPO_ICs = PO_ICs.SaturnEnceladus.CR3BP.L2_Lyapunov; famName = 'SaturnEnceladus.CR3BP.L2_Lyapunov';
% myPO_ICs = PO_ICs.SaturnEnceladus.CR3BP.L1_Lyapunov; famName = 'SaturnEnceladus.CR3BP.L1_Lyapunov';
% myPO_ICs = PO_ICs.SaturnEnceladus.CR3BP.L2_SHalo; famName = 'SaturnEnceladus.CR3BP.L2_SHalo';
% myPO_ICs = PO_ICs.SaturnEnceladus.CR3BP.L2_NHalo; famName = 'SaturnEnceladus.CR3BP.L2_NHalo';
% myPO_ICs = PO_ICs.SaturnEnceladus.CR3BP.L1_SHalo; famName = 'SaturnEnceladus.CR3BP.L1_SHalo';
% myPO_ICs = PO_ICs.SaturnEnceladus.CR3BP.L1_NHalo; famName = 'SaturnEnceladus.CR3BP.L1_NHalo';
% myPO_ICs = PO_ICs.SaturnEnceladus.CR3BP.L2_Vertical; famName = 'SaturnEnceladus.CR3BP.L2_Vertical';
% myPO_ICs = PO_ICs.SaturnEnceladus.CR3BP.L2_WesternAxial; famName = 'SaturnEnceladus.CR3BP.L2_WesternAxial';
% myPO_ICs = PO_ICs.SaturnEnceladus.CR3BP.L2_EasternAxial; famName = 'SaturnEnceladus.CR3BP.L2_EasternAxial';
% myPO_ICs = PO_ICs.SaturnEnceladus.CR3BP.DRO; famName = 'SaturnEnceladus.CR3BP.DRO';
% myPO_ICs = PO_ICs.SaturnEnceladus.CR3BP.L1_NHalo_DoublePeriod; famName = 'SaturnEnceladus.CR3BP.L1_NHalo_DoublePeriod';

% myPO_ICs = PO_ICs.SaturnEnceladus.CR3BP.unlabeled1LeadingS; famName = 'SaturnEnceladus.CR3BP.unlabeled1LeadingS';
% myPO_ICs = PO_ICs.SaturnEnceladus.CR3BP.unlabeled1TrailingS; famName = 'SaturnEnceladus.CR3BP.unlabeled1TrailingS';
% myPO_ICs = PO_ICs.SaturnEnceladus.CR3BP.unlabeled1LeadingN; famName = 'SaturnEnceladus.CR3BP.unlabeled1LeadingN';
% myPO_ICs = PO_ICs.SaturnEnceladus.CR3BP.unlabeled1TrailingN; famName = 'SaturnEnceladus.CR3BP.unlabeled1TrailingN';
% myPO_ICs = PO_ICs.SaturnEnceladus.CR3BP.SLeadingSneakerToe; famName = 'SaturnEnceladus.CR3BP.SLeadingSneakerToe';
myPO_ICs = PO_ICs.SaturnEnceladus.CR3BP.SLeadingSneakerHeel; famName = 'SaturnEnceladus.CR3BP.SLeadingSneakerHeel';
% myPO_ICs = PO_ICs.SaturnEnceladus.CR3BP.unknown1; famName = 'SaturnEnceladus.CR3BP.unknown1';


% myPO_ICs = PO_ICs.SaturnTitan.CR3BP.L2_Lyapunov; famName = 'SaturnTitan.CR3BP.L2_Lyapunov';
% myPO_ICs = PO_ICs.SaturnTitan.CR3BP.L2_Vertical; famName = 'SaturnTitan.CR3BP.L2_Vertical';

LP = 2;

% ds_PO = 1e-2;
% ds_PO = 1e-3; % Europa 
ds_PO = 1e-4;
% ds_PO = 1e-5;
% ds_PO = 1e-6; % Enceladus 
% ds_PO = 1e-7;
% ds_PO = 1e-8;
% ds_PO = 1e-9;

% error_tol = 1e-11; 
% error_tol = 1e-10;
error_tol = 1e-9;

familyDirection = -1; % 1 or -1

n_Nodes = 6; % 3

stepSize = 2; % 3 often seems to work well

iterMax = 100;

counter = 0;
lostFamilyCounter = 0;

plotSkip = 1; % Only plot every X POs
% ------------------------------------------------- 
%%% Preparing to find family
% -------------------------------------------------
X0_guess_n = myPO_ICs(1:6);
T_guess_n_new  = myPO_ICs(7);

%%% Preallocating
POs                     = NaN(n_POs_max,7);
stabilityIndices        = NaN(n_POs_max,2);
jacobiConstants         = NaN(n_POs_max,1);
L2ExcessVelocities_mps  = NaN(n_POs_max,1);

%%% Storing first PO
POs(1,:) = [X0_guess_n', T_guess_n_new];

%%% Storing first indices
JC_PO_i                     = JacobiConstantCalculator(mu,X0_guess_n(1:3)',X0_guess_n(4:6)');
L2_FlyoverVelocity_PO_i_mps = JC_2_L2FlyoverVelocity(JC_PO_i,mu,rLPs_n(LP,:),vNorm);
jacobiConstants(1)          = JC_PO_i;
L2ExcessVelocities_mps(1)   = L2_FlyoverVelocity_PO_i_mps;

%%% Initialize STM
stm0        = eye(6);
stm0_colVec = reshape(stm0,36,1);

% ------------------------------------------------- 
%%% Integrate reference and divide into nodes
% -------------------------------------------------
%%% Integrate reference PO
[Tref_n, Xref_n] = ode113(@Int_CR3BnSTM, [0, myPO_ICs(7)], [myPO_ICs(1:6); stm0_colVec], options, prms);

%%% Integrate again to discretize into nodes evenly spaced in time
[T_nodes, X_nodes] = ode113(@Int_CR3BnSTM, linspace(0, myPO_ICs(7),n_Nodes+1), [myPO_ICs(1:6); stm0_colVec], options, prms);

%%% Setting initial nodes and node-time
X_nodes = X_nodes(1:n_Nodes,1:6);
T_PO_i  = myPO_ICs(7);

%%% Plotting reference PO
if plot_current_PO == 1
     figure(1); hold all
     plot3(Xref_n(:,1),Xref_n(:,2),Xref_n(:,3),'b','linewidth',2)
     PlotBoi3('$X_n$','$Y_n$','$Z_n$',20,'LaTex')
     view(-10,50)
end % plot_current_PO

% ------------------------------------------------- 
%%% Preallocating before loop
% -------------------------------------------------
%%% Preallocate space for PO free variables
Fs      = NaN(n_POs_max, 6*n_Nodes + 1);

%%% Create first free-variable column vector
F_new   = X_nodes';
F_new   = F_new(:);
F_new   = [F_new; T_PO_i];
Fs(1,:) = F_new';

% ------------------------------------------------- 
%%% Enter loop to populate family of POs
% -------------------------------------------------
% for PO_i = 1:n_POs_max-1
PO_i = 1;
while PO_i < n_POs_max
    if print_PO_index == 1
        fprintf('PO_i = %1d ... iterations: %1d ... ds_PO = %1.2e\n',PO_i, counter, ds_PO)
    end % print_PO_index
    
    % --------------------------
    % Setup for new PO
    % --------------------------
    %%% Set free variable vector
    F_new = Fs(PO_i,:)';
    
    %%% Reset while-loop/convergence variables
    counter          = 0;
    constraint_error = 100;
    
    % --------------------------
    % Enter while loop - attempt to converge on next PO in family
    % --------------------------
    while (constraint_error > error_tol) && (counter < iterMax)
        
        %%% Count iteration
        counter = counter + 1;
        
        % --------------------------
        % Loop through nodes, integrate, and populate constraint vector and
        % create DF matrix
        % --------------------------
        %%% Preallocate constraint vector and DF matrix
        constraints = NaN(6*n_Nodes,1);
        DF          = zeros(6*n_Nodes,6*n_Nodes+1);
        
        %%% Loop through nodes
        for node_i = 1:n_Nodes
            %%% Integrate node            
            [~, X_node] = ode113(@Int_CR3BnSTM, [0, F_new(end)/n_Nodes], [F_new((6*(node_i-1)+1):(6*(node_i-1)+6)); stm0_colVec], options, prms);
            
            %%% Populate constraint vector
            if node_i < n_Nodes % If not the last node
                %%% Differencing end-state of current node with beginning
                %%% state of next node
                constraints((6*(node_i-1)+1):(6*(node_i-1)+6)) = X_node(end,1:6)' - F_new((node_i*6+1):(node_i*6+6));
                
                %%% Add Identity to the DF matrix
                DF((6*(node_i-1)+1):(6*(node_i-1)+6),(6*(node_i)+1):(6*(node_i)+6)) = -eye(6);

            elseif node_i == n_Nodes % If the last node
                %%% Differencing end-state of final node with beginning
                %%% of first node
                constraints((6*(node_i-1)+1):(6*(node_i-1)+6)) = X_node(end,1:6)' - F_new(1:6);

                %%% Add Identity to the DF matrix
                DF((6*(node_i-1)+1):(6*(node_i-1)+6),1:6) = -eye(6);

            end
            
            %%% Add STM to the DF matrix
            stm_tf_t0 = reshape(X_node(end,7:42),6,6);
            DF((6*(node_i-1)+1):(6*(node_i-1)+6),(6*(node_i-1)+1):(6*(node_i-1)+6)) = stm_tf_t0;
            
            %%% Add end-state dynamics to the DF matrix
            DF((6*(node_i-1)+1):(6*(node_i-1)+6),(6*n_Nodes+1)) = Int_CR3Bn(0,X_node(end,1:6)',prms);
            
        end % for node_i = 1:n_Nodes

        % --------------------------
        % If PO_i == 1, do a simplified version of things to find next F
        % --------------------------
        if PO_i == 1
            %%% Compute error
            constraint_error = norm(constraints);
            
            %%% Compute new free-variable vector if error is not converged
            if (constraint_error > error_tol)
                F_new = F_new - stepSize * DF'*((DF*(DF'))\constraints);
            end

        % --------------------------
        % If PO_i > 1, do full process for finding next F
        % --------------------------
        elseif PO_i > 1
            %%% Compute z
            z = (F_new - Fs(PO_i,:)')'*nullVecDF - ds_PO;
            
            %%% Compute augemented constraint vector
            constraints_aug = [constraints; z'];
            
            %%% Compute DH matrix
            DH = [DF; nullVecDF'];
            
            %%% Compute error
            constraint_error = norm(constraints_aug);
            
%             if (PO_i + 1) == 110
%                 989
%             end
                
            %%% Compute new free-variable vector if error is not converged
            if (constraint_error >= error_tol)
                F_new = F_new - (DH\constraints_aug);
            end

        end
        
    end % while constraint_error > error_tol
    
    % --------------------------
    % Checks
    % --------------------------


    % --------------------------
    % If checks are cleared and solution meets error tolerance, then
    % compute null space of DF, adjust step-size if necessary, and
    % store converged F
    % --------------------------
    if (constraint_error < error_tol)
        %%% Compute null space of DF
        if PO_i == 1
            nullVecDF = null(DF).*familyDirection;

        %%% If PO_i > 1, make sure null space hasn't flipped sign
        elseif PO_i > 1
            nullVecDFPrevious = nullVecDF;
            nullVecDF = null(DF).*familyDirection;
            
            %%% Check continuity of null space dimension and quit if it has
            %%% changed
            if isequal(size(nullVecDFPrevious),size(nullVecDF)) == 0
                warning('Null space changed dimensions')
                break
            end

            %%% Check sign
            nullDot = nullVecDFPrevious' * nullVecDF;

            %%% Modify null vector so there's no sign change
            for kk = 1:size(nullVecDF,2)
                if isequal(sign(nullDot(kk,kk)), -1)
                    nullVecDF(:,kk) = -1*nullVecDF(:,kk);
                end
            end
        end
   
        %%% Check for family continuity (no jumps to other equilibria)
        percentChangePosition = norm(F_new(1:3) - Fs(PO_i,1:3)')/norm(Fs(PO_i,1:3)');
        if (percentChangePosition > 0.1) || (norm(F_new(1:3) - rLPs_n(LP,:)') < 1e-14)
            if lostFamilyCounter < 5
                %%% Count the lost family iteration
                lostFamilyCounter = lostFamilyCounter + 1;

                %%% Decrease step size and try to find family again
                ds_PO = ds_PO * 0.9;
                continue
            
            %%% If algorithm repeated`ly fails to continue family, then end
            %%% the search
            elseif lostFamilyCounter >= 5
                warning('Lost track of family ... It possibly ended')
                break
            end
        end
        
        %%% Reset counter for when the algorithm loses the family
        lostFamilyCounter = 0;
        
        %%% Store converged F
        Fs(PO_i+1,:) = F_new';
        POs(PO_i+1,:) = [F_new(1:6)', F_new(end)];
        
        
        %%% Integrate new PO
        [~, X_new] = ode113(@Int_CR3BnSTM, [0, F_new(end)], [F_new(1:6); stm0_colVec], options, prms);
        
        %%% Plotting PO
        if plot_current_PO
            if mod(PO_i,plotSkip) == 0
                plot3(X_new(:,1),X_new(:,2),X_new(:,3),'r')
                drawnow
            end
        end
        
        %%% Stability indices of new PO
        stm_tf_t0 = reshape(X_new(end,7:42),6,6);
        monodromy = stm_tf_t0;
        [eigenVectors_new, eigenValues_new] = eig(monodromy);
        [S1, S2] = getStabilityIndices(diag(eigenValues_new));
        stabilityIndices(PO_i+1,:) = [S1, S2];
        
        %%% Energy of new PO
        JC_PO_new = JacobiConstantCalculator(mu,X_new(1,1:3),X_new(1,4:6));
        L2_FlyoverVelocity_mps_new = JC_2_L2FlyoverVelocity(JC_PO_new,secondary.MR,rLPs_n(LP,:),vNorm);
        
        jacobiConstants(PO_i+1)         = JC_PO_new;
        L2ExcessVelocities_mps(PO_i+1) = L2_FlyoverVelocity_mps_new;

        %%% Make corrections to dS
%         if counter <= 5
        if counter <= 30
            if isequal(secondary.name,'enceladus')
                if PO_i == 66
                    ds_PO = 9e-5;
%                     if PO_i > 137
%                         ds_PO = 0.00001;
%                     end
                else
                    ds_PO = 1.05 * ds_PO;
                end
                if PO_i == 192
                    break
                end
            elseif isequal(secondary.name,'titan')
                ds_PO = 1.2*ds_PO;
            else % if none of these bodies
                ds_PO = 1.1*ds_PO;
            end
        end
        
        %%% Update PO_i
        PO_i = PO_i + 1;
    else
        if PO_i == 1
            warning('Failed to converge on a 2nd PO')
        end
        
        %%% Lower the dS
        if isequal(secondary.name,'titan')
            ds_PO = 0.7*ds_PO;
        else % if none of these bodies
            ds_PO = 0.5*ds_PO;
        end
    end
    
end % while PO_i < n_POs_max

%% =======================================================================
%%% Clean up the data
% ========================================================================
% -------------------------------------------------
%%% Get rid of extra 'NaN' states
% -------------------------------------------------
%%% Grab all non-NaN data
nonNaN_ind = ~isnan(POs(:,1));
POs_trimmed = POs(nonNaN_ind,:);

% -------------------------------------------------
%%% Integrating solution to x-y plane crossing for state simplicity
% -------------------------------------------------
POs_final = NaN(size(POs_trimmed));

for kk = 1:size(POs_trimmed,1)
    [~, X_BCR_n_XY, time_event, X_event, index_event] = ...
        ode113(@Int_CR3Bn, [0, POs_trimmed(kk,end)], POs_trimmed(kk,1:6)', options_XZStop, prms);
    
    %%% Attempt, but if trajectory doesn't cross x-y plane, then forget it
    if isempty(X_event) == 0
        POs_final(kk,:) = [X_event(end,1:6),POs_trimmed(kk,end)];
    else
        POs_final(kk,:) = [POs_trimmed(kk,1:6),POs_trimmed(kk,end)];
    end
end
    

% ========================================================================
%%% Plots/Studies
% ========================================================================
% -------------------------------------------------
%%% Stability Indices
% -------------------------------------------------
if plot_stability == 1
    s1_ind = ~isnan(stabilityIndices(:,1));
    s2_ind = ~isnan(stabilityIndices(:,2));
    S1 = stabilityIndices(s1_ind,1);
    S2 = stabilityIndices(s2_ind,2);
    
    figure('position',[209 322 948 302])
    subplot(1,2,1); hold all
    p1 = plot(abs(S1),'o','markeredgecolor',colors.std.blue,'markerfacecolor',colors.std.ltblue);
    p2 = plot(abs(S2),'o','markeredgecolor',colors.std.red,'markerfacecolor',colors.std.ltred);
    plot([0 length(S1)],[2 2],'m','linewidth',0.5)
    PlotBoi2('PO index','Stability Indices',18,'LaTex')
    legend([p1 p2],'S_1','S_2')
    xlim([0 length(S1)])
    
    subplot(1,2,2); hold all
    p1 = plot(abs(S1),'o','markeredgecolor',colors.std.blue,'markerfacecolor',colors.std.ltblue);
    p2 = plot(abs(S2),'o','markeredgecolor',colors.std.red,'markerfacecolor',colors.std.ltred);
    plot([0 length(S1)],[2 2],'m','linewidth',0.5)
    PlotBoi2('PO index','',18,'LaTex')
    ylim([1.9 2.1])
    xlim([0 length(S1)])
end

% -------------------------------------------------
%%% Energy
% -------------------------------------------------
if plot_energy
    figure; hold all
    yyaxis left
    plot(jacobiConstants,'o','markeredgecolor',colors.std.blue,'markerfacecolor',colors.std.ltblue)
    PlotBoi2('','Jacobi Constant',18,'LaTex')
    yyaxis right
    plot(L2ExcessVelocities_mps,'o','markeredgecolor',colors.std.red,'markerfacecolor',colors.std.ltred)
    PlotBoi2('PO Index','$L_2$ Excess Velocity, $m/s$',18,'LaTex')
end






% ========================================================================
%%% Saving Data
% ========================================================================
if save_PO_X0_database == 1
    % -------------------------------------------------
    %%% Preparing Save File
    % -------------------------------------------------
    fileAlreadyExists = 1;
    fileVersion       = 0;
    famName_new       = famName;
    while fileAlreadyExists == 1
        if isfile([savePath,famName_new,'.txt']) == 1
            fileVersion = fileVersion + 1;
            famName_new = sprintf('%s_%1d',famName,fileVersion);
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
    headerString = ['x0_n,y0_n,z0_n,xd0_n,yd0_n,zd0_n,Tp,JC,L2ExcessVelocity_mps,...,',...
                    sprintf('error_tol=%1.0e,n_Nodes=%1d,stepsize=%1.1f\n',error_tol, n_Nodes, stepSize)];
    fprintf(datafile,headerString);
    
    %%% Write data
    dataStart = 1;
    dataStop = size(POs_final,1);
    
%     warning('Not full data set')
%     dataStart = 1;
%     dataStop = 123;
    
    for kk = dataStart:dataStop
        fprintf(datafile,'%1.16f,%1.16f,%1.16f,%1.16f,%1.16f,%1.16f,%1.16f,%1.10f,%1.3f\n',...
            POs_final(kk,1), POs_final(kk,2), POs_final(kk,3),...
            POs_final(kk,4), POs_final(kk,5), POs_final(kk,6),...
            POs_final(kk,end), jacobiConstants(kk), L2ExcessVelocities_mps(kk));
            
    end
%     for kk = 1:size(POs_trimmed,1)
%         fprintf(datafile,'%1.16f,%1.16f,%1.16f,%1.16f,%1.16f,%1.16f,%1.16f,%1.10f,%1.3f\n',...
%             POs_trimmed(kk,1), POs_trimmed(kk,2), POs_trimmed(kk,3),...
%             POs_trimmed(kk,4), POs_trimmed(kk,5), POs_trimmed(kk,6),...
%             POs_trimmed(kk,end), jacobiConstants(kk), L2ExcessVelocities_mps(kk));
%             
%     end
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




% % 
% % figure; hold all
% % 
% % PlotBoi3('$X_n$','$Y_n$','$Z_n$',20,'LaTex')
% % axis equal
% % plotSecondary(secondary)
% % for idx = 20
% %     X0 = POs(idx,1:6)';
% %      T0 = POs(idx,7);
% % 
% %      X0(3) = - X0(3);
% %      X0(6) = -X0(6);
% % 
% %     [T_n, X_n] = ode113(@Int_CR3Bn, linspace(0, T0,500), X0, options, prms);
% %     plot3(X_n(:,1),X_n(:,2),X_n(:,3),'b')
% %     
% % end
% % view(0,-90);



% 
% 
% X0m = X0 + [1;0;0;0;0;0].*1e-5;
% T0m = T0*10;
% [T_n, X_n] = ode113(@Int_CR3Bn, [T0m, 0], X0m, options, prms);
% plot3(X_n(:,1),X_n(:,2),X_n(:,3),'r')
% figure
% subplot(1,2,1); hold all
% plot3(X_n(:,1),X_n(:,2),X_n(:,3),'b')
% PlotBoi3('$X_n$','$Y_n$','$Z_n$',20,'LaTex')
% axis equal
% plotSecondary(secondary)
% 
% 
% subplot(1,2,2); hold all
% plot3(X_n(:,1),X_n(:,2),X_n(:,3),'r')
% PlotBoi3('$X_n$','$Y_n$','$Z_n$',20,'LaTex')
% axis equal
% plotSecondary(secondary)
% xlim([.993 1.0015])
% ylim([-1 1].*5e-3)
% zlim([-0.006 0.003])







% F = [1.003278945892397;
%  -0.000000000000000;
%  -0.001160178622759;
%  -0.008518652394546;
%  -0.009390206737848;
%  -0.004730702228170;
%  8.797668183419933];
% 
% [T_n, X_n] = ode113(@Int_CR3Bn, linspace(0, F(end),500), F(1:6), options, prms);
% 
%     
% 
% figure; hold all
% PlotBoi3('$X_n$','$Y_n$','$Z_n$',20,'LaTex')
% plot3(X_n(:,1),X_n(:,2),X_n(:,3),'b')
% axis equal
% plotSecondary(secondary)





