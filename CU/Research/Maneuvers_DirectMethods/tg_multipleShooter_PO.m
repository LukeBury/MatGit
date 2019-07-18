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
%%% Setup
% ========================================================================
% ------------------------------------------------- 
%%% Initial conditions
% -------------------------------------------------
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
% % X0 and Tp corresponding to L2 Vertical biffurcation at Enceladus
% % --------------------------
% primary = bodies.saturn;    secondary = bodies.enceladus;
% X0_n = [1.003104911775108; 0; -0.000000000003900;...
%         0.000000000011816; -0.003397299105638; 0.011089587166677];
% T0_PO_n = 4.248617774497027;




% % --------------------------
% % X0 and Tp corresponding to L4-untitled1 biffurcation at Enceladus
% % --------------------------
% primary = bodies.saturn;    secondary = bodies.enceladus;
% X0_n = [0.564947073975989; 0.302902430320629; -0.877900340143571;...
%         -0.454995222410374; -0.108438133722420; -0.238901417152178];
% T0_PO_n = 6.283185489835216;
% %^^^ may want a perturbation on the order of 1e-2









% % --------------------------
% % X0 and Tp corresponding to unlabeled1TrailingS biffurcation 1 at Enceladus
% % --------------------------
% primary = bodies.saturn;    secondary = bodies.enceladus;
% X0_n = [1.003867374501306; 0; -0.001615356503498;...
%         0.007420489615297; -0.010140084240107; 0.004290008196305];
% T0_PO_n = 4.449570709517759;

% --------------------------
% X0 and Tp corresponding to unlabeled1TrailingS biffurcation 2 at Enceladus
% --------------------------
primary = bodies.saturn;    secondary = bodies.enceladus;
X0_n = [1.003328461055074; 0; -0.001179557781950;...
        0.007734754444788; -0.009573295236007; 0.004683361178424];
T0_PO_n = 4.055805113650098;

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








% % --------------------------
% % X0 and Tp corresponding to L2 Vertical biffurcation at Titan
% % --------------------------
% primary = bodies.saturn;    secondary = bodies.titan;
% X0_n = [1.033297105540779; 0; 0.000000000000013;...
%         -0.000000000000012; -0.037907926186969; 0.117177360844365];
% T0_PO_n = 4.226170603679089;


% -------------------------------------------------
%%% System
% -------------------------------------------------
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
tol = 1e-12;
options = odeset('RelTol',tol,'AbsTol',tol);
prms.u    = mu;
prms.R2_n = R2_n;
t0 = 0;

% ========================================================================
%%% Find new family
% ========================================================================
% ------------------------------------------------- 
%%% Perturb X0 in z direction to find biffurcated family
% -------------------------------------------------
% X0_guess_n = X0_n - [0; 0; 2.1; 0; 0; 0].*1e-5;
% X0_guess_n = X0_n - [0; 0; 3; 0; 0; 0].*1e-3;
X0_guess_n = X0_n - [0; 0; 1; 0; 0; 0].*1e-3;

% ------------------------------------------------- 
%%% Integrate first guess
% -------------------------------------------------
%%% Initial stm is identity (t0 to t0)
stm0_vec = reshape(eye(6),36,1);

% X0_n = [X0_n; stm0];

% ------------------------------------------------- 
%%% Shooter Options
% -------------------------------------------------
n_iterations = 100;
error_tol = 1e-13;
dX0_scale = 1;

% -------------------------------------------------
%%% Initial integration
% -------------------------------------------------
%%% Integrate Reference
1
[time_n_ref, X_BCR_n_ref] = ode113(@Int_CR3Bn, [0 T0_PO_n], X0_n, options, prms);
2
%%% Integrate trajectory
time0_n = [0, T0_PO_n];
[time_n_guess, X_BCR_n_guess] = ode113(@Int_CR3BnSTM, time0_n, [X0_guess_n; stm0_vec], options, prms);
3






n_Nodes = 3;


% ------------------------------------------------- 
%%% Integrate reference and divide into nodes
% -------------------------------------------------
%%% Integrate again to discretize into evenly spaced nodes
[~, X_nodes] = ode113(@Int_CR3BnSTM, linspace(0, T0_PO_n,n_Nodes+1), [X0_guess_n; stm0_vec], options, prms);

%%% Setting initial nodes and node-time
X_nodes = X_nodes(1:n_Nodes,1:6);
T_PO_i  = time0_n(end);





% 
% % -------------------------------------------------
% %%% Discretize the guess into nodes
% % -------------------------------------------------
% n_Nodes = 3;
% [discretizedTrajectory] = discretize_by_path_distance(X_BCR_n_guess(:,1:6), n_Nodes + 1);
% X_nodes = discretizedTrajectory(1:n_Nodes,:);
% T_nodes = NaN(n_Nodes,1);
% for kk = 1:n_Nodes-1
% %     X_BCR_n_guess == X_nodes(kk,:)
%     [tf, index]=ismember(X_nodes(kk+1,:),X_BCR_n_guess(:,1:6),'rows');
%     
%     if kk == 1
%         T_nodes(kk) = time_n(index);
%     elseif kk > 1
%         T_nodes(kk) = time_n(index) - T_nodes(kk-1);
%     end
%     
%     if kk == (n_Nodes-1)
%         T_nodes(kk+1) = time_n(end) - time_n(index);
%     end
% end


% -------------------------------------------------
%%% Integrating Nodes
% -------------------------------------------------
%%% Integrate new state
[~, X_N1] = ode113(@Int_CR3BnSTM, [0, T_PO_i/3], [X_nodes(1,:)'; stm0_vec], options, prms);
[~, X_N2] = ode113(@Int_CR3BnSTM, [0, T_PO_i/3], [X_nodes(2,:)'; stm0_vec], options, prms);
[~, X_N3] = ode113(@Int_CR3BnSTM, [0, T_PO_i/3], [X_nodes(3,:)'; stm0_vec], options, prms);

% -------------------------------------------------
%%% Setting up for shooting
% -------------------------------------------------
%%% Initialize new variables to iterate over
warning('Need this to account for nodes')
% X0_new = X0_guess_n;
% X_BCR_n_new = X_BCR_n_guess;

F_new = X_nodes';
F_new = F_new(:);
F_new = [F_new; T_PO_i];

z66 = zeros(6,6);
z6 = zeros(6,1);

%%% Iterate
% for kk = 1:n_iterations
error = 100;
counter = 0;
while error > error_tol
    counter = counter + 1;
    % -------------------------------------------------
    %%% Correct state
    % -------------------------------------------------
    %%% Set constraint vector
%     c = X_BCR_n_new(end,1:6)' - X_BCR_n_new(1,1:6)';


    c = [X_N1(end,1:6)' - F_new(7:12);...
         X_N2(end,1:6)' - F_new(13:18);...
         X_N3(end,1:6)' - F_new(1:6)];
%     c = [X_N1(end,1:6)' - X_N2(1,1:6)';...
%          X_N2(end,1:6)' - X_N3(1,1:6)';...
%          X_N3(end,1:6)' - X_N1(1,1:6)'];
    
     
    %%% Grab STM(tf, t0)
%     stm_tf_t0 = reshape(X_BCR_n_new(end,7:42),6,6);
    stm1_tf_t0 = reshape(X_N1(end,7:42),6,6);
    stm2_tf_t0 = reshape(X_N2(end,7:42),6,6);
    stm3_tf_t0 = reshape(X_N3(end,7:42),6,6);
    
    %%% Grab EOMs at X0s and Xfs
    Xdot_N1_0 = Int_CR3Bn(0,X_N1(1,1:6)',prms);
    Xdot_N1_f = Int_CR3Bn(0,X_N1(end,1:6)',prms);
    Xdot_N2_0 = Int_CR3Bn(0,X_N2(1,1:6)',prms);
    Xdot_N2_f = Int_CR3Bn(0,X_N2(end,1:6)',prms);
    Xdot_N3_0 = Int_CR3Bn(0,X_N3(1,1:6)',prms);
    Xdot_N3_f = Int_CR3Bn(0,X_N3(end,1:6)',prms);
    
    %%% Create DF matrix
    DF = [stm1_tf_t0, -eye(6), z66, Xdot_N1_f;...
          z66, stm2_tf_t0, -eye(6), Xdot_N2_f;...
          -eye(6), z66, stm3_tf_t0, Xdot_N3_f];
%     DF = [stm1_tf_t0, -eye(6), z66, Xdot_N1_f, -Xdot_N2_0, z6;...
%           z66, stm2_tf_t0, -eye(6), z6, Xdot_N2_f, -Xdot_N3_0;...
%           -eye(6), z66, stm3_tf_t0, -Xdot_N1_0, z6, Xdot_N3_f];
    
    %%% Correct Free Variables
    F_old = F_new;
    F_new = F_new - (DF'/(DF*(DF')))*c;
%     F1_new = [F_new(1:6); F_new(19)];
%     F2_new = [F_new(7:12); F_new(20)];
%     F3_new = [F_new(13:18); F_new(21)];

    %%% Check for family continuity (no jumps to other equilibria)
        percentChangePosition = norm(F_new(1:3) - F_old(1:3))/norm(F_old(1:3));
        if percentChangePosition > 0.1
            warning('Lost track of family ... It possibly ended')
            break
        end
    

    % -------------------------------------------------
    %%% Integrate new correction
    % -------------------------------------------------
    %%% Integrate new state
%     [time_n, X_BCR_n_new] = ode113(@Int_CR3BnSTM, [0, T_PO_new], X0_new, options, prms);
    [~, X_N1] = ode113(@Int_CR3BnSTM, [0, F_new(end)/3], [F_new(1:6); stm0_vec], options, prms);
    [~, X_N2] = ode113(@Int_CR3BnSTM, [0, F_new(end)/3], [F_new(7:12); stm0_vec], options, prms);
    [~, X_N3] = ode113(@Int_CR3BnSTM, [0, F_new(end)/3], [F_new(13:18); stm0_vec], options, prms);
    
    error = norm(c);
    %%% Print results from current iteration
    fprintf('-------------------------------------------- %1d\n',counter)
    fprintf('Norm of constraint vector:  \t%1.4e\n',error)

end



% % % % -------------------------------------------------
% % % %%% Discretize the guess into nodes
% % % % -------------------------------------------------
% % % n_Nodes = 3;
% % % [discretizedTrajectory] = discretize_by_path_distance(X_BCR_n_guess(:,1:6), n_Nodes + 1);
% % % X_nodes = discretizedTrajectory(1:n_Nodes,:);
% % % T_nodes = NaN(n_Nodes,1);
% % % for kk = 1:n_Nodes-1
% % % %     X_BCR_n_guess == X_nodes(kk,:)
% % %     [tf, index]=ismember(X_nodes(kk+1,:),X_BCR_n_guess(:,1:6),'rows');
% % %     
% % %     if kk == 1
% % %         T_nodes(kk) = time_n(index);
% % %     elseif kk > 1
% % %         T_nodes(kk) = time_n(index) - T_nodes(kk-1);
% % %     end
% % %     
% % %     if kk == (n_Nodes-1)
% % %         T_nodes(kk+1) = time_n(end) - time_n(index);
% % %     end
% % % end
% % % 
% % % 
% % % % -------------------------------------------------
% % % %%% Integrating Nodes
% % % % -------------------------------------------------
% % % %%% Integrate new state
% % % [~, X_N1] = ode113(@Int_CR3BnSTM, [0, T_nodes(1)], [X_nodes(1,:)'; stm0_vec], options, prms);
% % % [~, X_N2] = ode113(@Int_CR3BnSTM, [0, T_nodes(2)], [X_nodes(2,:)'; stm0_vec], options, prms);
% % % [~, X_N3] = ode113(@Int_CR3BnSTM, [0, T_nodes(3)], [X_nodes(3,:)'; stm0_vec], options, prms);
% % % 
% % % % -------------------------------------------------
% % % %%% Setting up for shooting
% % % % -------------------------------------------------
% % % %%% Initialize new variables to iterate over
% % % warning('Need this to account for nodes')
% % % % X0_new = X0_guess_n;
% % % % X_BCR_n_new = X_BCR_n_guess;
% % % 
% % % F_new = X_nodes';
% % % F_new = F_new(:);
% % % F_new = [F_new; T_nodes];
% % % 
% % % z66 = zeros(6,6);
% % % z6 = zeros(6,1);
% % % 
% % % %%% Iterate
% % % % for kk = 1:n_iterations
% % % error = 100;
% % % counter = 0;
% % % while error > error_tol
% % %     counter = counter + 1;
% % %     % -------------------------------------------------
% % %     %%% Correct state
% % %     % -------------------------------------------------
% % %     %%% Set constraint vector
% % % %     c = X_BCR_n_new(end,1:6)' - X_BCR_n_new(1,1:6)';
% % %     c = [X_N1(end,1:6)' - X_N2(1,1:6)';...
% % %          X_N2(end,1:6)' - X_N3(1,1:6)';...
% % %          X_N3(end,1:6)' - X_N1(1,1:6)'];
% % % 
% % %     %%% Grab STM(tf, t0)
% % % %     stm_tf_t0 = reshape(X_BCR_n_new(end,7:42),6,6);
% % %     stm1_tf_t0 = reshape(X_N1(end,7:42),6,6);
% % %     stm2_tf_t0 = reshape(X_N2(end,7:42),6,6);
% % %     stm3_tf_t0 = reshape(X_N3(end,7:42),6,6);
% % %     
% % %     %%% Grab EOMs at X0s and Xfs
% % %     Xdot_N1_0 = Int_CR3Bn(0,X_N1(1,1:6)',prms);
% % %     Xdot_N1_f = Int_CR3Bn(0,X_N1(end,1:6)',prms);
% % %     Xdot_N2_0 = Int_CR3Bn(0,X_N2(1,1:6)',prms);
% % %     Xdot_N2_f = Int_CR3Bn(0,X_N2(end,1:6)',prms);
% % %     Xdot_N3_0 = Int_CR3Bn(0,X_N3(1,1:6)',prms);
% % %     Xdot_N3_f = Int_CR3Bn(0,X_N3(end,1:6)',prms);
% % %     
% % %     %%% Create DF matrix
% % % %     DF = [stm_tf_t0 - eye(6), Xdot_f];
% % %     DF = [stm1_tf_t0, -eye(6), z66, Xdot_N1_f, -Xdot_N2_0, z6;...
% % %           z66, stm2_tf_t0, -eye(6), z6, Xdot_N2_f, -Xdot_N3_0;...
% % %           -eye(6), z66, stm3_tf_t0, -Xdot_N1_0, z6, Xdot_N3_f];
% % %     
% % %     %%% Correct Free Variables
% % %     F_new = F_new - (DF'/(DF*(DF')))*c;
% % % %     X0_new = [F_new(1:6); stm0_vec];
% % %     F1_new = [F_new(1:6); F_new(19)];
% % %     F2_new = [F_new(7:12); F_new(20)];
% % %     F3_new = [F_new(13:18); F_new(21)];
% % %     
% % % %     T_PO_new = F_new(7);
% % % 
% % %     % -------------------------------------------------
% % %     %%% Integrate new correction
% % %     % -------------------------------------------------
% % %     %%% Integrate new state
% % % %     [time_n, X_BCR_n_new] = ode113(@Int_CR3BnSTM, [0, T_PO_new], X0_new, options, prms);
% % %     [~, X_N1] = ode113(@Int_CR3BnSTM, [0, F1_new(end)], [F1_new(1:6); stm0_vec], options, prms);
% % %     [~, X_N2] = ode113(@Int_CR3BnSTM, [0, F2_new(end)], [F2_new(1:6); stm0_vec], options, prms);
% % %     [~, X_N3] = ode113(@Int_CR3BnSTM, [0, F3_new(end)], [F3_new(1:6); stm0_vec], options, prms);
% % %     
% % %     error = norm(c);
% % %     %%% Print results from current iteration
% % %     fprintf('-------------------------------------------- %1d\n',counter)
% % %     fprintf('Norm of constraint vector:  \t%1.4e\n',error)
% % %     F1_new(3)
% % % 
% % % end


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
% [time_n_XY, X_BCR_n_XY] = ode113(@Int_CR3Bn, [0, T_PO_new], X0_new(1:6), options_XZStop, prms);
[time_n_XY, X_BCR_n_XY] = ode113(@Int_CR3Bn, [0, F_new(end)], F_new(1:6), options_XZStop, prms);



X_BCR_n_XY(end,1:6)'
F_new(end)

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

















