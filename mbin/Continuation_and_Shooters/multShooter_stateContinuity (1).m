function [DF_mat, constraints_vec] = multShooter_stateContinuity(n_Nodes, F, dynamicsWithSTM, options, prms)
%%% Description
%       Multiple shooter for state continuity of a 6-d state (ex, [x; y; z;
%       xd; yd; zd]). For use with pseudo arclength continuation
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

% -------------------------------------------------
%%% Preallocation
% -------------------------------------------------
constraints_vec = NaN(6*n_Nodes,1);
DF_mat          = zeros(6*n_Nodes,6*n_Nodes+1);

% -------------------------------------------------
%%% Loop through the nodes to populate DF and constraint vector
% -------------------------------------------------
%%% Loop through nodes
for node_i = 1:n_Nodes
    %%% Integrate node            
    [~, X_node] = ode113(dynamicsWithSTM, [0, F(end)/n_Nodes], [F((6*(node_i-1)+1):(6*(node_i-1)+6)); stm0_colVec], options, prms);

    %%% Populate constraint vector
    if node_i < n_Nodes % If not the last node
        %%% Differencing end-state of current node with beginning
        %%% state of next node
        constraints_vec((6*(node_i-1)+1):(6*(node_i-1)+6)) = X_node(end,1:6)' - F((node_i*6+1):(node_i*6+6));

        %%% Add Identity to the DF matrix
        DF_mat((6*(node_i-1)+1):(6*(node_i-1)+6),(6*(node_i)+1):(6*(node_i)+6)) = -eye(6);

    elseif node_i == n_Nodes % If the last node
        %%% Differencing end-state of final node with beginning
        %%% of first node
        constraints_vec((6*(node_i-1)+1):(6*(node_i-1)+6)) = X_node(end,1:6)' - F(1:6);

        %%% Add Identity to the DF matrix
        DF_mat((6*(node_i-1)+1):(6*(node_i-1)+6),1:6) = -eye(6);

    end

    %%% Add STM to the DF matrix
    stm_tf_t0 = reshape(X_node(end,7:42),6,6);
    DF_mat((6*(node_i-1)+1):(6*(node_i-1)+6),(6*(node_i-1)+1):(6*(node_i-1)+6)) = stm_tf_t0;

    %%% Add end-state dynamics to the DF matrix
    endstateDynamics = dynamicsWithSTM(0,[X_node(end,1:6)'; zeros(36,1)], prms);
    DF_mat((6*(node_i-1)+1):(6*(node_i-1)+6),(6*n_Nodes+1)) = endstateDynamics(1:6);

end % for node_i = 1:n_Nodes

end % function