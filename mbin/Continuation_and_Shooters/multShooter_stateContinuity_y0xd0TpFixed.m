function [DF_mat, constraints_vec] = multShooter_stateContinuity_y0xd0TpFixed(n_Nodes, F, dynamicsWithSTM, options, prms, y0, xd0, Tp_des)
%%% Description
%       Multiple shooter for state continuity of a 6-d state (ex, [x; y; z;
%       xd; yd; zd]) with a fixed time period and y0 of first node
%       
% ------------------------------------------------------------------------
%%% Inputs
%       n_Nodes         - [1x1] Number of nodes in multiple shooter 
%       F               - [nx1] Free-variable vector (stacked states of
%                           each node
%       dynamicsWithSTM - [@FunctionHandle] Function handle to dynamics 
%                           for ode integrator to use. Needs to include the
%                           dynamics for the state transition matrix
%       options         - [struct] options struct for ode integration 
%       prms            - [struct] parameters needed for integration (here,
%                           prms.u for the mass ratio)
%       y0              - [1x1] y0 of first node
%       xd0             - [1x1] xd0 of first node
%       Tp_des          - [1x1] Desired time period for solution
% ------------------------------------------------------------------------
%%% Outputs
%       DF_mat          - [(n_nodes*6)xn] Matrix of the partials of the
%                           constraints with respect to the free variables
%       constraints_vec - [(n_nodes*6)x1] Vector of constraints related to
%                           state continuity
% ------------------------------------------------------------------------
% Created: 10/02/19
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
DF_mat          = zeros(6*n_Nodes, 6*n_Nodes-2);

% -------------------------------------------------
%%% Loop through the nodes to populate DF and constraint vector
% -------------------------------------------------
nodeData = cell(n_Nodes);

parfor node_i = 1:n_Nodes
    constraints_vec_temp = zeros(6*n_Nodes, 1);
    DF_mat_temp          = zeros(6*n_Nodes, 6*n_Nodes-2);

    %%% Integrate node            
    if node_i == 1
        [~, X_node] = ode113(dynamicsWithSTM, [0, Tp_des/n_Nodes], [F(1); y0; F(2); xd0; F(3); F(4); stm0_colVec], options, prms);
    else
        [~, X_node] = ode113(dynamicsWithSTM, [0, Tp_des/n_Nodes], [F((6*(node_i-1)-1):(6*(node_i-1)+4)); stm0_colVec], options, prms);
    end

    %%% Populate constraint vector
    if node_i < n_Nodes % If not the last node
        
        %%% Differencing end-state of current node with beginning
        %%% state of next node
        constraints_vec_temp((6*(node_i-1)+1):(6*(node_i-1)+6)) = X_node(end,1:6)' - F((node_i*6-1):(node_i*6+4));

        %%% Add Identity to the DF matrix
        DF_mat_temp((6*(node_i-1)+1):(6*(node_i-1)+6),(6*(node_i-1)+5):(6*(node_i-1)+10)) = -eye(6);
    
    elseif node_i == n_Nodes % If the last node
        %%% Differencing end-state of final node with beginning
        %%% of first node
        constraints_vec_temp((6*(node_i-1)+1):(6*(node_i-1)+6)) = X_node(end,1:6)' - [F(1); y0; F(2); xd0; F(3); F(4)];

        %%% Add Identity to the DF matrix
        neye6 = -eye(6);
        DF_mat_temp((6*(node_i-1)+1):(6*(node_i-1)+6),1:4) = neye6(:,[1,3,5,6]);
        
    end
    
    %%% Add STM to the DF matrix
    stm_tf_t0 = reshape(X_node(end,7:42),6,6);
    if node_i == 1
        DF_mat_temp(1:6, 1:4) = stm_tf_t0(:,[1,3,5,6]);
    else
        DF_mat_temp((6*(node_i-1)+1):(6*(node_i-1)+6),(6*(node_i-1)-1):(6*(node_i-1)+4)) = stm_tf_t0;
    end
    
    nodeData{node_i}.constraints = constraints_vec_temp;
    nodeData{node_i}.DF          = DF_mat_temp;
end % for node_i = 1:n_Nodes

for node_i = 1:n_Nodes
    constraints_vec = constraints_vec + nodeData{node_i}.constraints;
    DF_mat          = DF_mat + nodeData{node_i}.DF;
end

end % function