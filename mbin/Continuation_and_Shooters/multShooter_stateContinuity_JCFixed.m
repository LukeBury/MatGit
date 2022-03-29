function [DF_mat, constraints_vec] = multShooter_stateContinuity_JCFixed(n_Nodes, F, dynamicsWithSTM, options, prms, JC_des)
%%% Description
%       Multiple shooter for state continuity of a 6-d state (ex, [x; y; z;
%       xd; yd; zd]) with a fixed jacobi constant
%       
% ------------------------------------------------------------------------
%%% Inputs
%       n_Nodes         - [1x1] Number of nodes in multiple shooter 
%       F               - [nx1] Free-variable vector (stacked states of
%                           each node with total time at end (so 
%                           n = n_nodes*6+1, generally))
%       dynamicsWithSTM - [@FunctionHandle] Function handle to dynamics 
%                           for ode integrator to use. Needs to include the
%                           dynamics for the state transition matrix
%       options         - [struct] options struct for ode integration 
%       prms            - [struct] parameters needed for integration (here,
%                           prms.u for the mass ratio)
%       JC_des          - [1x1] Desired jacobi constant for solution
% ------------------------------------------------------------------------
%%% Outputs
%       DF_mat          - [(n_nodes*6)xn] Matrix of the partials of the
%                           constraints with respect to the free variables
%       constraints_vec - [(n_nodes*6)x1] Vector of constraints related to
%                           state continuity
% ------------------------------------------------------------------------
% Created: 09/06/19
% Author : Luke Bury, luke.bury@colorado.edu
% ========================================================================
% -------------------------------------------------
%%% Creating necessary variables
% -------------------------------------------------
stm0_colVec = reshape(eye(6),36,1);

%%% Shortcuts
u = prms.u;

% -------------------------------------------------
%%% Preallocation
% -------------------------------------------------
constraints_vec = zeros(6*n_Nodes+1, 1);
DF_mat          = zeros(6*n_Nodes+1, 6*n_Nodes+1);

% -------------------------------------------------
%%% Loop through the nodes to populate DF and constraint vector
% -------------------------------------------------
%%% Loop through nodes
% for node_i = 1:n_Nodes
%     %%% Integrate node            
%     [~, X_node] = ode113(dynamicsWithSTM, [0, F(end)/n_Nodes], [F((6*(node_i-1)+1):(6*(node_i-1)+6)); stm0_colVec], options, prms);
% 
%     %%% Populate constraint vector
%     if node_i < n_Nodes % If not the last node
%         if node_i == 1 % If first node, store x10 for later use
%             x10  = X_node(1,1);
%             y10  = X_node(1,2);
%             z10  = X_node(1,3);
%             xd10 = X_node(1,4);
%             yd10 = X_node(1,5);
%             zd10 = X_node(1,6);
%             r1   = ((x10 + u)^2 + y10^2 + z10^2)^(1/2);
%             r2   = ((x10 - 1 + u)^2 + y10^2 + z10^2)^(1/2);
%             
%             JC_X10 = getJacobiConstant_ZH(X_node(1,1:6),prms);
%         end
%         
%         %%% Differencing end-state of current node with beginning
%         %%% state of next node
%         constraints_vec((6*(node_i-1)+1):(6*(node_i-1)+6)) = X_node(end,1:6)' - F((node_i*6+1):(node_i*6+6));
% 
%         %%% Add Identity to the DF matrix
%         DF_mat((6*(node_i-1)+1):(6*(node_i-1)+6),(6*(node_i)+1):(6*(node_i)+6)) = -eye(6);
%     
%     elseif node_i == n_Nodes % If the last node
%         %%% Differencing end-state of final node with beginning
%         %%% of first node
%         constraints_vec((6*(node_i-1)+1):(6*(node_i-1)+6)) = X_node(end,1:6)' - F(1:6);
% 
%         %%% Add Identity to the DF matrix
%         DF_mat((6*(node_i-1)+1):(6*(node_i-1)+6),1:6) = -eye(6);
%         
%     end
%     
%     %%% Add STM to the DF matrix
%     stm_tf_t0 = reshape(X_node(end,7:42),6,6);
%     DF_mat((6*(node_i-1)+1):(6*(node_i-1)+6),(6*(node_i-1)+1):(6*(node_i-1)+6)) = stm_tf_t0;
% 
%     %%% Add end-state dynamics to the DF matrix
%     endstateDynamics = dynamicsWithSTM(0,[X_node(end,1:6)'; zeros(36,1)], prms);
%     DF_mat((6*(node_i-1)+1):(6*(node_i-1)+6),(6*n_Nodes+1)) = endstateDynamics(1:6);
%     
%     %%% Filling in last row of constraint vector and DF (based on JC
%     %%% constraint)
%     constraints_vec(6*n_Nodes + 1) = JC_des - JC_X10;
%     
%     %%% Filling in last row of DF matrix with -JC partials
%     if isfield(prms, 'J2p')
%         dJCJ2p_dx = (prms.J2p*prms.R1^2*(2*u + 2*x10)*(u - 1))/((u + x10)^2 + y10^2 + z10^2)^(5/2) - (5*prms.J2p*prms.R1^2*(2*u + 2*x10)*(u - 1)*((u + x10)^2 + y10^2 - 2*z10^2))/(2*((u + x10)^2 + y10^2 + z10^2)^(7/2));
%         dJCJ2p_dy = (2*prms.J2p*prms.R1^2*y10*(u - 1))/((u + x10)^2 + y10^2 + z10^2)^(5/2) - (5*prms.J2p*prms.R1^2*y10*(u - 1)*((u + x10)^2 + y10^2 - 2*z10^2))/((u + x10)^2 + y10^2 + z10^2)^(7/2);
%         dJCJ2p_dz = -(4*prms.J2p*prms.R1^2*z10*(u - 1))/((u + x10)^2 + y10^2 + z10^2)^(5/2) - (5*prms.J2p*prms.R1^2*z10*(u - 1)*((u + x10)^2 + y10^2 - 2*z10^2))/((u + x10)^2 + y10^2 + z10^2)^(7/2);
%     else
%         dJCJ2p_dx = 0;
%         dJCJ2p_dy = 0;
%         dJCJ2p_dz = 0;
%     end
% 
%     if isfield(prms, 'J2s')
%         dJCJ2s_dx = (5*prms.J2s*prms.R2^2*u*(2*u + 2*x10 - 2)*((u + x10 - 1)^2 + y10^2 - 2*z10^2))/(2*((u + x10 - 1)^2 + y10^2 + z10^2)^(7/2)) - (prms.J2s*prms.R2^2*u*(2*u + 2*x10 - 2))/((u + x10 - 1)^2 + y10^2 + z10^2)^(5/2);
%         dJCJ2s_dy = (5*prms.J2s*prms.R2^2*u*y10*((u + x10 - 1)^2 + y10^2 - 2*z10^2))/((u + x10 - 1)^2 + y10^2 + z10^2)^(7/2) - (2*prms.J2s*prms.R2^2*u*y10)/((u + x10 - 1)^2 + y10^2 + z10^2)^(5/2);
%         dJCJ2s_dz = (4*prms.J2s*prms.R2^2*u*z10)/((u + x10 - 1)^2 + y10^2 + z10^2)^(5/2) + (5*prms.J2s*prms.R2^2*u*z10*((u + x10 - 1)^2 + y10^2 - 2*z10^2))/((u + x10 - 1)^2 + y10^2 + z10^2)^(7/2);
%     else
%         dJCJ2s_dx = 0;
%         dJCJ2s_dy = 0;
%         dJCJ2s_dz = 0;
%     end
%     
%     if isfield(prms, 'J3p')
%         warning('Haven''t added this equation y10et')
%     end
% 
%     if isfield(prms, 'J3s')
%         warning('Haven''t added this equation y10et')
%     end
% 
%     if isfield(prms, 'J4p')
%         dJCJ4p_dx = (9*prms.J4p*prms.R1^4*(2*u + 2*x10)*(u - 1)*(3*((u + x10)^2 + y10^2 + z10^2)^2 + 35*z10^4 - 30*z10^2*((u + x10)^2 + y10^2 + z10^2)))/(8*((u + x10)^2 + y10^2 + z10^2)^(11/2)) - (prms.J4p*prms.R1^4*(6*(2*u + 2*x10)*((u + x10)^2 + y10^2 + z10^2) - 30*z10^2*(2*u + 2*x10))*(u - 1))/(4*((u + x10)^2 + y10^2 + z10^2)^(9/2));
%         dJCJ4p_dy = (9*prms.J4p*prms.R1^4*y10*(u - 1)*(3*((u + x10)^2 + y10^2 + z10^2)^2 + 35*z10^4 - 30*z10^2*((u + x10)^2 + y10^2 + z10^2)))/(4*((u + x10)^2 + y10^2 + z10^2)^(11/2)) - (prms.J4p*prms.R1^4*(12*y10*((u + x10)^2 + y10^2 + z10^2) - 60*y10*z10^2)*(u - 1))/(4*((u + x10)^2 + y10^2 + z10^2)^(9/2));
%         dJCJ4p_dz = (prms.J4p*prms.R1^4*(48*z10*((u + x10)^2 + y10^2 + z10^2) - 80*z10^3)*(u - 1))/(4*((u + x10)^2 + y10^2 + z10^2)^(9/2)) + (9*prms.J4p*prms.R1^4*z10*(u - 1)*(3*((u + x10)^2 + y10^2 + z10^2)^2 + 35*z10^4 - 30*z10^2*((u + x10)^2 + y10^2 + z10^2)))/(4*((u + x10)^2 + y10^2 + z10^2)^(11/2));
%     else
%         dJCJ4p_dx = 0;
%         dJCJ4p_dy = 0;
%         dJCJ4p_dz = 0;
%     end
% 
%     if isfield(prms, 'J4s')
%         warning('Haven''t added this equation y10et')
%     end
% 
%     if isfield(prms, 'J5p')
%         warning('Haven''t added this equation y10et')
%     end
% 
%     if isfield(prms, 'J5s')
%         warning('Haven''t added this equation y10et')
%     end
%     
%     if isfield(prms, 'J6p')
%         dJCJ6p_dx = (prms.J6p*prms.R1^6*(u - 1)*(315*z10^4*(2*u + 2*x10) + 15*(2*u + 2*x10)*((u + x10)^2 + y10^2 + z10^2)^2 - 210*z10^2*(2*u + 2*x10)*((u + x10)^2 + y10^2 + z10^2)))/(8*((u + x10)^2 + y10^2 + z10^2)^(13/2)) - (13*prms.J6p*prms.R1^6*(2*u + 2*x10)*(u - 1)*(5*((u + x10)^2 + y10^2 + z10^2)^3 - 105*z10^2*((u + x10)^2 + y10^2 + z10^2)^2 - 231*z10^6 + 315*z10^4*((u + x10)^2 + y10^2 + z10^2)))/(16*((u + x10)^2 + y10^2 + z10^2)^(15/2));
%         dJCJ6p_dy = (prms.J6p*prms.R1^6*(u - 1)*(630*y10*z10^4 + 30*y10*((u + x10)^2 + y10^2 + z10^2)^2 - 420*y10*z10^2*((u + x10)^2 + y10^2 + z10^2)))/(8*((u + x10)^2 + y10^2 + z10^2)^(13/2)) - (13*prms.J6p*prms.R1^6*y10*(u - 1)*(5*((u + x10)^2 + y10^2 + z10^2)^3 - 105*z10^2*((u + x10)^2 + y10^2 + z10^2)^2 - 231*z10^6 + 315*z10^4*((u + x10)^2 + y10^2 + z10^2)))/(8*((u + x10)^2 + y10^2 + z10^2)^(15/2));
%         dJCJ6p_dz = -(prms.J6p*prms.R1^6*(u - 1)*(756*z10^5 + 180*z10*((u + x10)^2 + y10^2 + z10^2)^2 - 840*z10^3*((u + x10)^2 + y10^2 + z10^2)))/(8*((u + x10)^2 + y10^2 + z10^2)^(13/2)) - (13*prms.J6p*prms.R1^6*z10*(u - 1)*(5*((u + x10)^2 + y10^2 + z10^2)^3 - 105*z10^2*((u + x10)^2 + y10^2 + z10^2)^2 - 231*z10^6 + 315*z10^4*((u + x10)^2 + y10^2 + z10^2)))/(8*((u + x10)^2 + y10^2 + z10^2)^(15/2));
%     else
%         dJCJ6p_dx = 0;
%         dJCJ6p_dy = 0;
%         dJCJ6p_dz = 0;
%     end
% 
%     if isfield(prms, 'J6s')
%         warning('Haven''t added this equation yet')
%     end
%     
%     
%     DF_mat(6*n_Nodes + 1, 1) = -2*x10 + (1-u)*2*(x10 + u)/(r1^3) + u*2*(x10 - 1 + u)/(r2^3) + dJCJ2p_dx + dJCJ2s_dx + dJCJ4p_dx + dJCJ6p_dx;
%     DF_mat(6*n_Nodes + 1, 2) = -2*y10 + (1-u)*2*y10/(r1^3) + u*2*y10/(r2^3)                 + dJCJ2p_dy + dJCJ2s_dy + dJCJ4p_dy + dJCJ6p_dy;
%     DF_mat(6*n_Nodes + 1, 3) = (1-u)*2*z10/(r1^3) + u*2*z10/(r2^3)                          + dJCJ2p_dz + dJCJ2s_dz + dJCJ4p_dz + dJCJ6p_dz;
%     DF_mat(6*n_Nodes + 1, 4) = 2*xd10;
%     DF_mat(6*n_Nodes + 1, 5) = 2*yd10;
%     DF_mat(6*n_Nodes + 1, 6) = 2*zd10;
% 
% end % for node_i = 1:n_Nodes










% x10  = F(1);
% y10  = F(2);
% z10  = F(3);
% xd10 = F(4);
% yd10 = F(5);
% zd10 = F(6);
r1   = ((F(1) + u)^2 + F(2)^2 + F(3)^2)^(1/2);
r2   = ((F(1) - 1 + u)^2 + F(2)^2 + F(3)^2)^(1/2);
JC_X10 = getJacobiConstant_ZH(F(1:6)',prms);

nodeData = cell(n_Nodes);

parfor node_i = 1:n_Nodes
    
    constraints_vec_temp = zeros(6*n_Nodes+1, 1);
    DF_mat_temp          = zeros(6*n_Nodes+1, 6*n_Nodes+1);
    
    %%% Integrate node            
    [~, X_node] = ode113(dynamicsWithSTM, [0, F(end)/n_Nodes], [F((6*(node_i-1)+1):(6*(node_i-1)+6)); stm0_colVec], options, prms);

    %%% Populate constraint vector
    if node_i < n_Nodes % If not the last node

        %%% Differencing end-state of current node with beginning
        %%% state of next node
        constraints_vec_temp((6*(node_i-1)+1):(6*(node_i-1)+6)) = X_node(end,1:6)' - F((node_i*6+1):(node_i*6+6));

        %%% Add Identity to the DF matrix
        DF_mat_temp((6*(node_i-1)+1):(6*(node_i-1)+6),(6*(node_i)+1):(6*(node_i)+6)) = -eye(6);
    
    elseif node_i == n_Nodes % If the last node
        %%% Differencing end-state of final node with beginning
        %%% of first node
        constraints_vec_temp((6*(node_i-1)+1):(6*(node_i-1)+6)) = X_node(end,1:6)' - F(1:6);

        %%% Add Identity to the DF matrix
        DF_mat_temp((6*(node_i-1)+1):(6*(node_i-1)+6),1:6) = -eye(6);
        
    end
    
    %%% Add STM to the DF matrix
    stm_tf_t0 = reshape(X_node(end,7:42),6,6);
    DF_mat_temp((6*(node_i-1)+1):(6*(node_i-1)+6),(6*(node_i-1)+1):(6*(node_i-1)+6)) = stm_tf_t0;

    %%% Add end-state dynamics to the DF matrix
    endstateDynamics = dynamicsWithSTM(0,[X_node(end,1:6)'; zeros(36,1)], prms);
    DF_mat_temp((6*(node_i-1)+1):(6*(node_i-1)+6),(6*n_Nodes+1)) = endstateDynamics(1:6);
    
    %%% Filling in last row of constraint vector and DF (based on JC
    %%% constraint)
    constraints_vec_temp(6*n_Nodes + 1) = JC_des - JC_X10;
    
%     %%% Filling in last row of DF matrix with -JC partials
%     if isfield(prms, 'J2p')
%         dJCJ2p_dx = (prms.J2p*prms.R1^2*(2*u + 2*x10)*(u - 1))/((u + x10)^2 + y10^2 + z10^2)^(5/2) - (5*prms.J2p*prms.R1^2*(2*u + 2*x10)*(u - 1)*((u + x10)^2 + y10^2 - 2*z10^2))/(2*((u + x10)^2 + y10^2 + z10^2)^(7/2));
%         dJCJ2p_dy = (2*prms.J2p*prms.R1^2*y10*(u - 1))/((u + x10)^2 + y10^2 + z10^2)^(5/2) - (5*prms.J2p*prms.R1^2*y10*(u - 1)*((u + x10)^2 + y10^2 - 2*z10^2))/((u + x10)^2 + y10^2 + z10^2)^(7/2);
%         dJCJ2p_dz = -(4*prms.J2p*prms.R1^2*z10*(u - 1))/((u + x10)^2 + y10^2 + z10^2)^(5/2) - (5*prms.J2p*prms.R1^2*z10*(u - 1)*((u + x10)^2 + y10^2 - 2*z10^2))/((u + x10)^2 + y10^2 + z10^2)^(7/2);
%     else
%         dJCJ2p_dx = 0;
%         dJCJ2p_dy = 0;
%         dJCJ2p_dz = 0;
%     end
% 
%     if isfield(prms, 'J2s')
%         dJCJ2s_dx = (5*prms.J2s*prms.R2^2*u*(2*u + 2*x10 - 2)*((u + x10 - 1)^2 + y10^2 - 2*z10^2))/(2*((u + x10 - 1)^2 + y10^2 + z10^2)^(7/2)) - (prms.J2s*prms.R2^2*u*(2*u + 2*x10 - 2))/((u + x10 - 1)^2 + y10^2 + z10^2)^(5/2);
%         dJCJ2s_dy = (5*prms.J2s*prms.R2^2*u*y10*((u + x10 - 1)^2 + y10^2 - 2*z10^2))/((u + x10 - 1)^2 + y10^2 + z10^2)^(7/2) - (2*prms.J2s*prms.R2^2*u*y10)/((u + x10 - 1)^2 + y10^2 + z10^2)^(5/2);
%         dJCJ2s_dz = (4*prms.J2s*prms.R2^2*u*z10)/((u + x10 - 1)^2 + y10^2 + z10^2)^(5/2) + (5*prms.J2s*prms.R2^2*u*z10*((u + x10 - 1)^2 + y10^2 - 2*z10^2))/((u + x10 - 1)^2 + y10^2 + z10^2)^(7/2);
%     else
%         dJCJ2s_dx = 0;
%         dJCJ2s_dy = 0;
%         dJCJ2s_dz = 0;
%     end
%     
%     if isfield(prms, 'J3p')
%         warning('Haven''t added this equation y10et')
%     end
% 
%     if isfield(prms, 'J3s')
%         warning('Haven''t added this equation y10et')
%     end
% 
%     if isfield(prms, 'J4p')
%         dJCJ4p_dx = (9*prms.J4p*prms.R1^4*(2*u + 2*x10)*(u - 1)*(3*((u + x10)^2 + y10^2 + z10^2)^2 + 35*z10^4 - 30*z10^2*((u + x10)^2 + y10^2 + z10^2)))/(8*((u + x10)^2 + y10^2 + z10^2)^(11/2)) - (prms.J4p*prms.R1^4*(6*(2*u + 2*x10)*((u + x10)^2 + y10^2 + z10^2) - 30*z10^2*(2*u + 2*x10))*(u - 1))/(4*((u + x10)^2 + y10^2 + z10^2)^(9/2));
%         dJCJ4p_dy = (9*prms.J4p*prms.R1^4*y10*(u - 1)*(3*((u + x10)^2 + y10^2 + z10^2)^2 + 35*z10^4 - 30*z10^2*((u + x10)^2 + y10^2 + z10^2)))/(4*((u + x10)^2 + y10^2 + z10^2)^(11/2)) - (prms.J4p*prms.R1^4*(12*y10*((u + x10)^2 + y10^2 + z10^2) - 60*y10*z10^2)*(u - 1))/(4*((u + x10)^2 + y10^2 + z10^2)^(9/2));
%         dJCJ4p_dz = (prms.J4p*prms.R1^4*(48*z10*((u + x10)^2 + y10^2 + z10^2) - 80*z10^3)*(u - 1))/(4*((u + x10)^2 + y10^2 + z10^2)^(9/2)) + (9*prms.J4p*prms.R1^4*z10*(u - 1)*(3*((u + x10)^2 + y10^2 + z10^2)^2 + 35*z10^4 - 30*z10^2*((u + x10)^2 + y10^2 + z10^2)))/(4*((u + x10)^2 + y10^2 + z10^2)^(11/2));
%     else
%         dJCJ4p_dx = 0;
%         dJCJ4p_dy = 0;
%         dJCJ4p_dz = 0;
%     end
% 
%     if isfield(prms, 'J4s')
%         warning('Haven''t added this equation y10et')
%     end
% 
%     if isfield(prms, 'J5p')
%         warning('Haven''t added this equation y10et')
%     end
% 
%     if isfield(prms, 'J5s')
%         warning('Haven''t added this equation y10et')
%     end
%     
%     if isfield(prms, 'J6p')
%         dJCJ6p_dx = (prms.J6p*prms.R1^6*(u - 1)*(315*z10^4*(2*u + 2*x10) + 15*(2*u + 2*x10)*((u + x10)^2 + y10^2 + z10^2)^2 - 210*z10^2*(2*u + 2*x10)*((u + x10)^2 + y10^2 + z10^2)))/(8*((u + x10)^2 + y10^2 + z10^2)^(13/2)) - (13*prms.J6p*prms.R1^6*(2*u + 2*x10)*(u - 1)*(5*((u + x10)^2 + y10^2 + z10^2)^3 - 105*z10^2*((u + x10)^2 + y10^2 + z10^2)^2 - 231*z10^6 + 315*z10^4*((u + x10)^2 + y10^2 + z10^2)))/(16*((u + x10)^2 + y10^2 + z10^2)^(15/2));
%         dJCJ6p_dy = (prms.J6p*prms.R1^6*(u - 1)*(630*y10*z10^4 + 30*y10*((u + x10)^2 + y10^2 + z10^2)^2 - 420*y10*z10^2*((u + x10)^2 + y10^2 + z10^2)))/(8*((u + x10)^2 + y10^2 + z10^2)^(13/2)) - (13*prms.J6p*prms.R1^6*y10*(u - 1)*(5*((u + x10)^2 + y10^2 + z10^2)^3 - 105*z10^2*((u + x10)^2 + y10^2 + z10^2)^2 - 231*z10^6 + 315*z10^4*((u + x10)^2 + y10^2 + z10^2)))/(8*((u + x10)^2 + y10^2 + z10^2)^(15/2));
%         dJCJ6p_dz = -(prms.J6p*prms.R1^6*(u - 1)*(756*z10^5 + 180*z10*((u + x10)^2 + y10^2 + z10^2)^2 - 840*z10^3*((u + x10)^2 + y10^2 + z10^2)))/(8*((u + x10)^2 + y10^2 + z10^2)^(13/2)) - (13*prms.J6p*prms.R1^6*z10*(u - 1)*(5*((u + x10)^2 + y10^2 + z10^2)^3 - 105*z10^2*((u + x10)^2 + y10^2 + z10^2)^2 - 231*z10^6 + 315*z10^4*((u + x10)^2 + y10^2 + z10^2)))/(8*((u + x10)^2 + y10^2 + z10^2)^(15/2));
%     else
%         dJCJ6p_dx = 0;
%         dJCJ6p_dy = 0;
%         dJCJ6p_dz = 0;
%     end
% 
%     if isfield(prms, 'J6s')
%         warning('Haven''t added this equation yet')
%     end
    
    
%     DF_mat_temp(6*n_Nodes + 1, 1) = -2*x10 + (1-u)*2*(x10 + u)/(r1^3) + u*2*(x10 - 1 + u)/(r2^3) + dJCJ2p_dx + dJCJ2s_dx + dJCJ4p_dx + dJCJ6p_dx;
%     DF_mat_temp(6*n_Nodes + 1, 2) = -2*y10 + (1-u)*2*y10/(r1^3) + u*2*y10/(r2^3)                 + dJCJ2p_dy + dJCJ2s_dy + dJCJ4p_dy + dJCJ6p_dy;
%     DF_mat_temp(6*n_Nodes + 1, 3) = (1-u)*2*z10/(r1^3) + u*2*z10/(r2^3)                          + dJCJ2p_dz + dJCJ2s_dz + dJCJ4p_dz + dJCJ6p_dz;
%     DF_mat_temp(6*n_Nodes + 1, 4) = 2*xd10;
%     DF_mat_temp(6*n_Nodes + 1, 5) = 2*yd10;
%     DF_mat_temp(6*n_Nodes + 1, 6) = 2*zd10;
    DF_mat_temp(6*n_Nodes + 1, 1) = -2*F(1) + (1-u)*2*(F(1) + u)/(r1^3) + u*2*(F(1) - 1 + u)/(r2^3);
    DF_mat_temp(6*n_Nodes + 1, 2) = -2*F(2) + (1-u)*2*F(2)/(r1^3) + u*2*F(2)/(r2^3);
    DF_mat_temp(6*n_Nodes + 1, 3) = (1-u)*2*F(3)/(r1^3) + u*2*F(3)/(r2^3);
    DF_mat_temp(6*n_Nodes + 1, 4) = 2*F(4);
    DF_mat_temp(6*n_Nodes + 1, 5) = 2*F(5);
    DF_mat_temp(6*n_Nodes + 1, 6) = 2*F(6);    
    nodeData{node_i}.constraints = constraints_vec_temp;
    nodeData{node_i}.DF          = DF_mat_temp;

end % for node_i = 1:n_Nodes

for node_i = 1:n_Nodes
    constraints_vec = constraints_vec + nodeData{node_i}.constraints;
    DF_mat          = DF_mat + nodeData{node_i}.DF;
end



end % function