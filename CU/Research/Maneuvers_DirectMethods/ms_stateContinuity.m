% ========================================================================
%%% Description
% ========================================================================
% Finding POs by also correcting time vector

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

% ========================================================================
%%% Run Switches
% ========================================================================


% ======================================================================== 
%%% Choose Initial conditions
% ========================================================================
% --------------------------
% X0 and Tp corresponding to L2 Lyapunov biffurcation at Europa
% --------------------------
% primary = bodies.jupiter;   secondary = bodies.europa;
% X0_n = [1.017026611070739; 0.000000014441540; 0;...
%         0.000000002451836; 0.020022207280852; 0];
% T0_PO_n = 3.124263922352081;


% % --------------------------
% % X0 and Tp corresponding to L2 Lyapunov biffurcation at Moon
% % --------------------------
% primary = bodies.earth;     secondary = bodies.moon;
% X0_n = [1.120360542366735; 0.000036610925873; 0;...
%         0.000002187136222; 0.176125048984917; 0];
% T0_PO_n = 3.415558008018916;


% % --------------------------
% % X0 and Tp corresponding to L2 Vertical biffurcation at Europa
% % --------------------------
% primary = bodies.jupiter;   secondary = bodies.europa;
% X0_n = [1.015816379745207; 0.000000000411681; 0.000000001297984;...
%         0.000000002061456; -0.017939765578591; -0.056562126699769];
% T0_PO_n = 4.266711475227815;
% 
% % --------------------------
% % X0 and Tp corresponding to L2 Lyapunov biffurcation at Enceladus
% % --------------------------
% primary = bodies.saturn;    secondary = bodies.enceladus;
% X0_n = [1.003390444557048; -0.000805481191125; 0;...
%         -0.000224307203511; 0.003235968095636; 0];
% T0_PO_n = 3.086292283602011;
% 

% % --------------------------
% % X0 and Tp corresponding to L1 Lyapunov biffurcation at Enceladus
% % --------------------------
% primary = bodies.saturn;    secondary = bodies.enceladus;
% X0_n = [0.996640782195733; -0.000000000000000; 0.000000000000000; 
%         0.000000000000020; -0.003745177417673; 0.000000000000000];
% T0_PO_n = 3.070641955511902;


% % --------------------------
% % X0 and Tp corresponding to L2 Vertical biffurcation at Enceladus
% % --------------------------
% primary = bodies.saturn;    secondary = bodies.enceladus;
% X0_n = [1.003104911775108; 0; -0.000000000003900;...
%         0.000000000011816; -0.003397299105638; 0.011089587166677];
% T0_PO_n = 4.248617774497027;

% % --------------------------
% % X0 and Tp corresponding to L1NHalo bifurcation at Enceladus
% % --------------------------
% primary = bodies.saturn;    secondary = bodies.enceladus;
% X0_n = [1.000004091375796;
%  -0.000000000000000;
%  0.000957917115369;
%  -0.000000000000769;
%  -0.018968378843832;
%  -0.000000000001634];
% T0_PO_n = 2.129860133440525;




% % --------------------------
% % X0 and Tp corresponding to L4-untitled1 biffurcation at Enceladus
% % --------------------------
primary = bodies.saturn;    secondary = bodies.enceladus;

F0 = [1.003220319415447;
 -0.000000000000000;
 -0.001100356595935;
 -0.007831042158205;
 -0.009467123528172;
 -0.004790308571176;
 3.976603869623910];

X0_n = F0(1:6);
T0_PO_n = F0(7);
% %^^^ may want a perturbation on the order of 1e-2


% [0.999995639364235;
%  -0.000000000000000;
%  -0.000997088724235;
%  -0.000000000000463;
%  -0.018557749287345;
%  0.000000000000901;
%  2.156603173730311];





% % --------------------------
% % X0 and Tp corresponding to unlabeled1TrailingS biffurcation 1 at Enceladus
% % --------------------------
% primary = bodies.saturn;    secondary = bodies.enceladus;
% X0_n = [1.003867374501306; 0; -0.001615356503498;...
%         0.007420489615297; -0.010140084240107; 0.004290008196305];
% T0_PO_n = 4.449570709517759;

% % % --------------------------
% % % X0 and Tp corresponding to unlabeled1TrailingS biffurcation 2 at Enceladus
% % % --------------------------
% % primary = bodies.saturn;    secondary = bodies.enceladus;
% % X0_n = [1.003328461055074; 0; -0.001179557781950;...
% %         0.007734754444788; -0.009573295236007; 0.004683361178424];
% % T0_PO_n = 4.055805113650098;

% % --------------------------
% % X0 and Tp corresponding to unlabeled1LeadingS biffurcation 1 at Enceladus
% % --------------------------
% primary = bodies.saturn;    secondary = bodies.enceladus;
% X0_n = [; 0; ;...
%         ; ; ];
% T0_PO_n = ;

% % --------------------------
% % X0 and Tp corresponding to unlabeled1LeadingS biffurcation 2 at Enceladus
% % --------------------------
% primary = bodies.saturn;    secondary = bodies.enceladus;
% X0_n = [; 0; ;...
%         ; ; ];
% T0_PO_n = ;


% From unlabeled1TrailingS
% 1.003814095574193
% -0.000000000000000
% -0.001569451247257
% 0.007439897701095
% -0.010081037303422
% 0.004319399519643
% 4.411213942244500
% 
% 
% 1.003225933063576
% -0.000000000000000
% -0.001104395924098
% 0.007825754765614
% -0.009472577747200
% 0.004784503685013
% 3.980705501546650



% From unlabeled1LeadingS
% [1.003820317696885;
% -0.000000000000000;
% -0.001574783007593;
% -0.007437504977750;
% -0.010087896968718;
% -0.004315867086289;
%  4.415704819325041];
% 
% [1.003245763478651;
% -0.000000000000000;
% -0.001118728443973;
% -0.007807332249385;
% -0.009491894560994;
% -0.004764222611185;
%  3.995204356308705];

% % --------------------------
% % X0 and Tp corresponding to L2 Vertical biffurcation at Titan
% % --------------------------
% primary = bodies.saturn;    secondary = bodies.titan;
% X0_n = [1.033297105540779; 0; 0.000000000000013;...
%         -0.000000000000012; -0.037907926186969; 0.117177360844365];
% T0_PO_n = 4.226170603679089;



% ========================================================================
%%% Setup
% ========================================================================
% -------------------------------------------------
%%% System
% -------------------------------------------------
%%% Normalizing constants
rNorm = secondary.a;         % n <-> km
tNorm = 1/secondary.meanMot; % n <-> sec
vNorm = rNorm / tNorm;       % n <-> km/sec

%%% Shortcut variables
prms.u   = secondary.MR;
prms.R2_n = secondary.R_n;

%%% Equillibrium Points
rLPs_n = EquilibriumPoints(prms.u);

% ------------------------------------------------- 
%%% Integration Options
% -------------------------------------------------
tol = 1e-12;
options = odeset('RelTol',tol,'AbsTol',tol);
t0 = 0;

% ========================================================================
%%% Find new family
% ========================================================================
% ------------------------------------------------- 
%%% Perturb X0 in z direction to find biffurcated family
% -------------------------------------------------
% X0_guess_n = X0_n - [0; 0; 5; 0; 0; 0].*1e-7;
% X0_guess_n = X0_n - [0; 0; 0; 0; 3; 0].*1e-4;
X0_guess_n = X0_n - [0; 0; 0; 1; 1; 0].*(6).*1e-5;
% X0_guess_n = X0_n - [0; 0; 4; 0; 0; 0].*1e-3;
% X0_guess_n = X0_n + [1; 0; 6; 0; 0; 0].*(1).*1e-6;
% X0_guess_n = X0_n;

T0_guess_n = T0_PO_n*2;
% T0_guess_n = T0_PO_n;
% ------------------------------------------------- 
%%% Integrate first guess
% -------------------------------------------------
%%% Initial stm is identity (t0 to t0)
stm0_colVec = reshape(eye(6),36,1);

% ------------------------------------------------- 
%%% Shooter Options
% -------------------------------------------------
n_Nodes = 6;
iterMax = 100;
error_tol = 1e-13;
stepSize = 1;
maxLostSolutionAttempts = 10;

%%% Turning off nearly-singular-matrix warning
warning_id = 'MATLAB:nearlySingularMatrix';
warning('off',warning_id)
% -------------------------------------------------
%%% Initial integration
% -------------------------------------------------
%%% Integrate Reference
[time_n_ref, X_BCR_n_ref] = ode113(@Int_CR3Bn, [0 T0_PO_n], X0_n, options, prms);

% ------------------------------------------------- 
%%% Integrate guess and divide into nodes
% -------------------------------------------------
%%% Integrate again to discretize into evenly spaced nodes
[~, X_nodes0] = ode113(@Int_CR3BnSTM, linspace(0, T0_guess_n,n_Nodes+1), [X0_guess_n; stm0_colVec], options, prms);

%%% Setting initial nodes and node-time
X_nodes0 = X_nodes0(1:n_Nodes,1:6);

% ------------------------------------------------- 
%%% Preallocating before loop
% -------------------------------------------------
%%% Create first free-variable column vector
F_new = X_nodes0';
F_new = F_new(:);
F_new = [F_new; T0_guess_n];
F0 = F_new;

% ------------------------------------------------- 
%%% Enter loop to find new PO
% -------------------------------------------------
counter           = 0;
lostSolutionCounter = 0;
constraint_error  = 100;
while (constraint_error > error_tol) && (counter < iterMax)
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
        
%         989
%         figure(1); hold all
%         plot3(X_node(:,1),X_node(:,2),X_node(:,3))

    end % for node_i = 1:n_Nodes
    
    %%% Compute error
    constraint_error = norm(constraints);
    
    %%% Compute new free-variable vector if error is not converged
    if (constraint_error > error_tol)
        F_new = F_new - stepSize * DF'*((DF*(DF'))\constraints);
    end
    
    %%% Check for solution continuity (no jumps to other equilibria)
    percentChangePosition = norm(F_new(1:3) - F0(1:3))/norm(F0(1:3));
    if (percentChangePosition > 0.4) 
        %%% Count the lost family iteration
        lostSolutionCounter = lostSolutionCounter + 1;
        
        if (lostSolutionCounter == 1)
            %%% Decrease perturbation size and recompute the guess
            perturbation0    = X0_guess_n - X0_n;
            perturbation_new = perturbation0 * 0.9;
            
        elseif (lostSolutionCounter > 1) && (lostSolutionCounter < maxLostSolutionAttempts)
            perturbation_new = perturbation_new * 0.9;
            
        elseif lostSolutionCounter >=maxLostSolutionAttempts
            warning('Not finding a nearby solution')
            break
        end
        
        %%% Print a status update
        fprintf('\n*** Lost track of the solution - trying again with smaller perturbation\n\n')
        
        %%% Use new perturbation to obtain a new guess and restart the
        %%% search
        X0_guess_n = X0_n + perturbation_new;
        [~, X_nodes0] = ode113(@Int_CR3BnSTM, linspace(0, T0_guess_n,n_Nodes+1), [X0_guess_n; stm0_colVec], options, prms);
        X_nodes0 = X_nodes0(1:n_Nodes,1:6);
        F_new = X_nodes0';
        F_new = F_new(:);
        F_new = [F_new; T0_guess_n];
        counter = 0;
        continue
    end
% %     %%% Check for family continuity (no jumps to other equilibria)
% %         percentChangePosition = norm(F_new(1:3) - Fs(PO_i,1:3)')/norm(Fs(PO_i,1:3)');
% %         if (percentChangePosition > 0.1) || (norm(F_new(1:3) - rLPs_n(LP,:)') < 1e-14)
% %             if lostFamilyCounter < 5
% %                 %%% Count the lost family iteration
% %                 lostFamilyCounter = lostFamilyCounter + 1;
% % 
% %                 %%% Decrease step size and try to find family again
% %                 ds_PO = ds_PO * 0.9;
% %                 continue
% %             
% %             %%% If algorithm repeatedly fails to continue family, then end
% %             %%% the search
% %             elseif lostFamilyCounter > 5
% %                 warning('Lost track of family ... It possibly ended')
% %                 break
% %             end
% %         end
    
    fprintf('iteration: %1d ... constraint error: %1.0e\n',counter, constraint_error)
%     fprintf('Norm of constraint vector:  \t%1.4e\n',)
    
end

%%% Turning nearly-singular-matrix warning back on
warning('on',warning_id);



% -------------------------------------------------
%%% Integrating solution to see what it looks like
% -------------------------------------------------
[time_n_full, X_BCR_n_full] = ode113(@Int_CR3Bn, [0, F_new(end)], F_new(1:6), options, prms);
figure; hold all
p1 = plot3(X_BCR_n_full(:,1),X_BCR_n_full(:,2),X_BCR_n_full(:,3),'m');
p2 = plot3(X_BCR_n_ref(:,1),X_BCR_n_ref(:,2),X_BCR_n_ref(:,3),'k');
legend([p2, p1],'ref','new')
PlotBoi3('$X_n$','$Y_n$','$Z_n$',20,'LaTex')

% -------------------------------------------------
%%% Integrating solution to x-y plane crossing for state simplicity
% -------------------------------------------------
options_XZStop = odeset('Event',@event_yEqualsZero,'RelTol',tol,'AbsTol',tol);
[time_n_XY, X_BCR_n_XY] = ode113(@Int_CR3Bn, [0, F_new(end)], F_new(1:6), options_XZStop, prms);


fprintf('\nNew Initial Conditions:\n')
for kk = 1:6
    if kk == 1
        fprintf('[')
        fprintf('%1.15f;\n',X_BCR_n_XY(end,kk))
    else
        fprintf(' %1.15f;\n',X_BCR_n_XY(end,kk))
    end
end
fprintf(' %1.15f];\n\n',F_new(end))

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
toc(ticWhole)

















