function [DF_mat, constraints_vec] = multShooter_stateContinuity_yzxdzdTpFixed(n_Nodes, F, dynamicsWithSTM, options, prms, Tp_des)
%%% Description
%       Multiple shooter for state continuity of a 6-d state (ex, [x; y; z;
%       xd; yd; zd]) with a fixed time period, and y=z=xd=zd=0. (Good for
%       Lyapunov POs)
%       
% ------------------------------------------------------------------------
%%% Inputs
%       n_Nodes         - [1x1] Number of nodes in multiple shooter 
%       F               - [(6N-4)x1] Free-variable vector (stacked states of
%                           each node
%       dynamicsWithSTM - [@FunctionHandle] Function handle to dynamics 
%                           for ode integrator to use. Needs to include the
%                           dynamics for the state transition matrix
%       options         - [struct] options struct for ode integration 
%       prms            - [struct] parameters needed for integration (here,
%                           prms.u for the mass ratio)
%       Tp_des          - [1x1] Desired time period for solution
% ------------------------------------------------------------------------
%%% Outputs
%       DF_mat          - [(6N)x(6N-4)] Matrix of the partials of the
%                           constraints with respect to the free variables
%       constraints_vec - [(6N)x1] Vector of constraints related to
%                           state continuity
% ------------------------------------------------------------------------
% Created: 7/10/20
% Author : Luke Bury, luke.bury@colorado.edu
% ========================================================================
% -------------------------------------------------
%%% Creating necessary variables
% -------------------------------------------------
stm0_colVec = reshape(eye(6),36,1);

% %%% Shortcuts
% u = prms.u;

% -------------------------------------------------
%%% Preallocation
% -------------------------------------------------
constraints_vec = zeros(6*n_Nodes, 1);
DF_mat          = zeros(6*n_Nodes, 6*n_Nodes-4);

% -------------------------------------------------
%%% Loop through the nodes to populate DF and constraint vector
% -------------------------------------------------

nodeData = cell(n_Nodes);

[~, X1_node] = ode113(dynamicsWithSTM, [0, Tp_des/n_Nodes], [F(1); 0; 0; 0; F(2); 0; stm0_colVec], options, prms);

%%% Loop through nodes
parfor node_i = 1:n_Nodes
    constraints_vec_temp = zeros(6*n_Nodes, 1);
    DF_mat_temp          = zeros(6*n_Nodes, 6*n_Nodes-4);

    %%% Integrate node            
    if node_i == 1
        X_node = X1_node;
    else
        [~, X_node] = ode113(dynamicsWithSTM, [0, Tp_des/n_Nodes], [F((6*(node_i-1)-3):(6*(node_i-1)+2)); stm0_colVec], options, prms);
    end
    
    %%% Populate constraint vector
    if node_i < n_Nodes % If not the last node
        
        %%% Differencing end-state of current node with beginning
        %%% state of next node
        constraints_vec_temp((6*(node_i-1)+1):(6*(node_i-1)+6)) = X_node(end,1:6)' - F((node_i*6-3):(node_i*6+2));

        %%% Add Identity to the DF matrix
        DF_mat_temp((6*(node_i-1)+1):(6*(node_i-1)+6),(6*(node_i)-3):(6*(node_i)+2)) = -eye(6);
    
    elseif node_i == n_Nodes % If the last node
        %%% Differencing end-state of final node with beginning
        %%% of first node
        constraints_vec_temp((6*(node_i-1)+1):(6*(node_i-1)+6)) = X_node(end,1:6)' - X1_node(1,1:6)';

        %%% Add Identity to the DF matrix
        neye6 = -eye(6);
        DF_mat_temp((6*(node_i-1)+1):(6*(node_i-1)+6),1:2) = [-1, 0; 0, 0; 0, 0; 0, 0; 0, -1; 0, 0];
        
    end
    
    %%% Add STM to the DF matrix
    stm_tf_t0 = reshape(X_node(end,7:42),6,6);
    if node_i == 1
        DF_mat_temp(1:6, 1) = stm_tf_t0(:,1);
        DF_mat_temp(1:6, 2) = stm_tf_t0(:,5);
    else
        DF_mat_temp((6*(node_i-1)+1):(6*(node_i-1)+6),(6*(node_i-1)-3):(6*(node_i-1)+2)) = stm_tf_t0;
    end
    
    nodeData{node_i}.constraints = constraints_vec_temp;
    nodeData{node_i}.DF          = DF_mat_temp;
end % for node_i = 1:n_Nodes

for node_i = 1:n_Nodes
    constraints_vec = constraints_vec + nodeData{node_i}.constraints;
    DF_mat          = DF_mat + nodeData{node_i}.DF;
end

end % function