function [DF_mat, constraints_vec] = multShooter_stateContinuity_JCy0Fixed(n_Nodes, F, dynamicsWithSTM, options, prms, JC_des, y0)
%%% Description
%       Multiple shooter for state continuity of a 6-d state (ex, [x; y; z;
%       xd; yd; zd]) with a fixed jacobi constant and y0 of the first node 
%       
% ------------------------------------------------------------------------
%%% Inputs
%       n_Nodes         - [scalar] Number of nodes in multiple shooter 
%       F               - [6Nx1] Free-variable vector (stacked states of
%                           each node with total time at end)
%       dynamicsWithSTM - [@FunctionHandle] Function handle to dynamics 
%                           for ode integrator to use. Needs to include the
%                           dynamics for the state transition matrix
%       options         - [struct] options struct for ode integration 
%       prms            - [struct] parameters needed for integration (here,
%                           prms.u for the mass ratio)
%       JC_des          - [1x1] Desired jacobi constant for solution
%       y0              - [scalar] Value for y0 of first node
% ------------------------------------------------------------------------
%%% Outputs
%       DF_mat          - [(6N+1)x6N] Matrix of the partials of the
%                           constraints with respect to the free variables
%       constraints_vec - [(6N+1)x1] Vector of constraints related to
%                           state continuity
% ------------------------------------------------------------------------
% Created: 04/06/21
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
DF_mat          = zeros(6*n_Nodes+1, 6*n_Nodes);

%%% Initialize cell structure to temporarily hold data from each parallel 
%%% node
nodeData = cell(n_Nodes);

%%% Determine r1, r2, and JC of the first node
% x10  = F(1);
% y10  = y0;
% z10  = F(2);
% xd10 = F(3);
% yd10 = F(4);
% zd10 = F(5);
r1   = ((F(1) + u)^2 + y0^2 + F(2)^2)^(1/2);
r2   = ((F(1) - 1 + u)^2 + y0^2 + F(2)^2)^(1/2);
JC_X10 = getJacobiConstant_ZH([F(1); y0; F(2:5)]',prms);
% -------------------------------------------------
%%% Loop through the nodes to populate DF and constraint vector
% -------------------------------------------------
parfor node_i = 1:n_Nodes
    
    constraints_vec_temp = zeros(6*n_Nodes+1, 1);
    DF_mat_temp          = zeros(6*n_Nodes+1, 6*n_Nodes);
    
    %%% Integrate node            
    if node_i == 1
        [~, X_node] = ode113(dynamicsWithSTM, [0, F(end)/n_Nodes], [F(1); y0; F(2:5); stm0_colVec], options, prms);
    else
        [~, X_node] = ode113(dynamicsWithSTM, [0, F(end)/n_Nodes], [F((6*(node_i-1)):(6*(node_i-1)+5)); stm0_colVec], options, prms);
    end

    %%% Populate constraint vector
    if node_i < n_Nodes % If not the last node

        %%% Differencing end-state of current node with beginning
        %%% state of next node
        constraints_vec_temp((6*(node_i-1)+1):(6*(node_i-1)+6)) = X_node(end,1:6)' - F((node_i*6):(node_i*6+5));

        %%% Add Identity to the DF matrix
        DF_mat_temp((6*(node_i-1)+1):(6*(node_i-1)+6),(6*(node_i-1)+6):(6*(node_i-1)+11)) = -eye(6);
    
    elseif node_i == n_Nodes % If the last node
        %%% Differencing end-state of final node with beginning
        %%% of first node
        constraints_vec_temp((6*(node_i-1)+1):(6*(node_i-1)+6)) = X_node(end,1:6)' - [F(1); y0; F(2:5)];

        %%% Add Identity to the DF matrix
        neye6 = -eye(6);
        DF_mat_temp((6*(node_i-1)+1):(6*(node_i-1)+6),1:5) = neye6(1:6,[1, 3:6]);

    end
    
    %%% Add STM to the DF matrix
    stm_tf_t0 = reshape(X_node(end,7:42),6,6);
    if node_i == 1
        DF_mat_temp(1:6, 1:5) = stm_tf_t0(:,[1,3:6]);
    else
        DF_mat_temp((6*(node_i-1)+1):(6*(node_i-1)+6),(6*(node_i-1)):(6*(node_i-1)+5)) = stm_tf_t0;
    end

    %%% Add end-state dynamics to the DF matrix
    endstateDynamics = dynamicsWithSTM(0,[X_node(end,1:6)'; zeros(36,1)], prms);
    DF_mat_temp((6*(node_i-1)+1):(6*(node_i-1)+6),(6*n_Nodes)) = endstateDynamics(1:6);
    
    %%% Filling in last row of constraint vector and DF (based on JC
    %%% constraint)
    constraints_vec_temp(6*n_Nodes + 1) = JC_des - JC_X10;
    
    %%% Add the JC partials to the DF matrix
    DF_mat_temp(6*n_Nodes + 1, 1) = -2*F(1) + (1-u)*2*(F(1) + u)/(r1^3) + u*2*(F(1) - 1 + u)/(r2^3);
    DF_mat_temp(6*n_Nodes + 1, 2) = (1-u)*2*F(2)/(r1^3) + u*2*F(2)/(r2^3);
    DF_mat_temp(6*n_Nodes + 1, 3) = 2*F(3);
    DF_mat_temp(6*n_Nodes + 1, 4) = 2*F(4);
    DF_mat_temp(6*n_Nodes + 1, 5) = 2*F(5);    
    
    %%% Collect the temporary data which is specific to this node in the
    %%% cell structure that holds data for all nodes
    nodeData{node_i}.constraints = constraints_vec_temp;
    nodeData{node_i}.DF          = DF_mat_temp;

end % for node_i = 1:n_Nodes

for node_i = 1:n_Nodes
    constraints_vec = constraints_vec + nodeData{node_i}.constraints;
    DF_mat          = DF_mat + nodeData{node_i}.DF;
end



end % function