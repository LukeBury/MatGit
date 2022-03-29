function [DF_mat, constraints_vec] = multShooter_stateContinuity_x0y0Fixed(n_Nodes, F, dynamicsWithSTM, options, prms, x0, y0)
%%% Description
%       Multiple shooter for state continuity of a 6-d state (ex, [x; y; z;
%       xd; yd; zd]) with x0 and y0 of the first node fixed
%       
% ------------------------------------------------------------------------
%%% Inputs
%       n_Nodes         - [scalar] Number of nodes in multiple shooter 
%       F               - [nx1] Free-variable vector (stacked states of
%                           each node with total time at end)
%       dynamicsWithSTM - [@FunctionHandle] Function handle to dynamics 
%                           for ode integrator to use. Needs to include the
%                           dynamics for the state transition matrix
%       options         - [struct] options struct for ode integration 
%       prms            - [struct] parameters needed for integration (here,
%                           prms.u for the mass ratio)
%       x0              - [scalar] Value for x0 of first node
%       y0              - [scalar] Value for y0 of first node
% ------------------------------------------------------------------------
%%% Outputs
%       DF_mat          - [(n_nodes*6)xn] Matrix of the partials of the
%                           constraints with respect to the free variables
%       constraints_vec - [(n_nodes*6)x1] Vector of constraints related to
%                           state continuity
% ------------------------------------------------------------------------
% Created: 02/03/21
% Author : Luke Bury, luke.bury@colorado.edu
% ========================================================================
% -------------------------------------------------
%%% Creating necessary variables
% -------------------------------------------------
stm0_colVec = reshape(eye(6),36,1);

% -------------------------------------------------
%%% Preallocation
% -------------------------------------------------
constraints_vec = zeros(6*n_Nodes,1);
DF_mat          = zeros(6*n_Nodes,6*n_Nodes-1);

%%% Initialize cell structure to temporarily hold data from each parallel 
%%% node
nodeData = cell(n_Nodes);

% -------------------------------------------------
%%% Loop through the nodes to populate DF and constraint vector
% -------------------------------------------------
%%% Loop through nodes
parfor node_i = 1:n_Nodes
    %%% Temporary constraints and DF structures that will only gather the
    %%% data for the current node
    constraints_vec_temp = zeros(6*n_Nodes, 1);
    DF_mat_temp          = zeros(6*n_Nodes, 6*n_Nodes-1);
    
    %%% Integrate node            
    if node_i == 1
        [~, X_node] = ode113(dynamicsWithSTM, [0, F(end)/n_Nodes], [x0; y0; F(1:4); stm0_colVec], options, prms);
    else
        [~, X_node] = ode113(dynamicsWithSTM, [0, F(end)/n_Nodes], [F((6*(node_i-1)-1):(6*(node_i-1)+4)); stm0_colVec], options, prms);
    end
    
    %%% Populate temporary constraint vector
    if node_i < n_Nodes % If not the last node
        %%% Differencing end-state of current node with beginning
        %%% state of next node
        constraints_vec_temp((6*(node_i-1)+1):(6*(node_i-1)+6)) = X_node(end,1:6)' - F((node_i*6-1):(node_i*6+4));

        %%% Add Identity to the temporary DF matrix
        DF_mat_temp((6*(node_i-1)+1):(6*(node_i-1)+6),(6*(node_i-1)+5):(6*(node_i-1)+10)) = -eye(6);

    elseif node_i == n_Nodes % If the last node
        %%% Differencing end-state of final node with beginning
        %%% of first node
        constraints_vec_temp((6*(node_i-1)+1):(6*(node_i-1)+6)) = X_node(end,1:6)' - [x0; y0; F(1:4)];

        %%% Add Identity to the temporary DF matrix
        neye6 = -eye(6);
        DF_mat_temp((6*(node_i-1)+1):(6*(node_i-1)+6),1:4) = neye6(1:6,3:6);
    end

    %%% Add STM to the temporary DF matrix
    stm_tf_t0 = reshape(X_node(end,7:42),6,6);
    
    if node_i == 1
        DF_mat_temp(1:6, 1:4) = stm_tf_t0(:,3:6);
    else
        DF_mat_temp((6*(node_i-1)+1):(6*(node_i-1)+6),(6*(node_i-1)-1):(6*(node_i-1)+4)) = stm_tf_t0;
    end

    %%% Add end-state dynamics to the temporary DF matrix
    endstateDynamics = dynamicsWithSTM(0,[X_node(end,1:6)'; zeros(36,1)], prms);
    DF_mat_temp((6*(node_i-1)+1):(6*(node_i-1)+6),(6*n_Nodes-1)) = endstateDynamics(1:6);
    
    %%% Collect the temporary data which is specific to this node in the
    %%% cell structure that holds data for all nodes
    nodeData{node_i}.constraints = constraints_vec_temp;
    nodeData{node_i}.DF          = DF_mat_temp;
end % for node_i = 1:n_Nodes

%%% Combine data from cell structure into single constraints vector and DF
%%% matrix
for node_i = 1:n_Nodes
    constraints_vec = constraints_vec + nodeData{node_i}.constraints;
    DF_mat          = DF_mat + nodeData{node_i}.DF;
end

end % function








