function [PO_new, counter, constraint_error] = correctPO_mS_sC_JCy0Fixed(PO_guess, n_Nodes, integrator, integrator_options, prms, constraint_error_tol, ms_iterMax, stepSize, JC_des,y0)
%%% Description
% Multiple shooting algorithm to correct a periodic orbit for state
% continuity at a target JC with a fixed y0
%       
% ------------------------------------------------------------------------
%%% Inputs
% PO_guess             - [7x1] Vector of guess state and time period
% n_Nodes              - [scalar] number of nodes for multiple shooter
% integrator           - [func] Handle to desired integrator, 
%                         ex: '@Int_CR3BnSTM'
% integrator_options   - 'options' for ode
% prms                 - [struct] fields (u: mass ratio, n: normalized 
%                         mean motio)
% constraint_error_tol - [scalar] Error tolerance for norm of constraint
%                         vector
% ms_iterMax           - [scalar] Max number of iterations for multiple 
%                         shooter to converge in
% stepSize             - [scalar] Artificial scaling down of update that
%                         may be useful in particularly sensitive systems
% JC_des               - [scalar] target jacobi constant
% y0                   - [scalar] fixed y0 value of first node
% ------------------------------------------------------------------------
%%% Outputs
% PO_new           - [7x1] Corrected PO
% counter          - [scalar] Iterations required to converge
% constraint_error - [scalar] Norm of final constraint vector
% ------------------------------------------------------------------------
% Created: 12/11/20
% Author : Luke Bury, luke.bury@colorado.edu
% ========================================================================
%%% Acquire nodes
if n_Nodes ==1
    X_nodes = PO_guess(1:6)';
else
    stm0_colvec = reshape(eye(6),36,1);
    [~, X_nodes] = get_nodes([PO_guess(1:6); stm0_colvec], [0, PO_guess(7)], n_Nodes+1, integrator, integrator_options, prms);
    X_nodes = X_nodes(1:n_Nodes,1:6);
end

%%% Create free-variable vector
F_new   = reshape(X_nodes',6*n_Nodes,1);
F_new   = [F_new; PO_guess(7)];
F_new   = [F_new(1); F_new(3:end)];

%%% Initialize iteration counter and error value
counter = 0;
constraint_error = 100;

%%% Enter multiple-shooting loop to correct guess
while (constraint_error > constraint_error_tol) && (counter < ms_iterMax)
    %%% Count iteration
    counter = counter + 1;

    % --------------------------
    % Using multiple shooting - loop through nodes, integrate, and 
    % populate constraint vector and create DF matrix
    % --------------------------
    [DF_mat, constraints] = multShooter_stateContinuity_JCy0Fixed(n_Nodes, F_new, integrator, integrator_options, prms, JC_des, y0);
    
    %%% Compute error
    constraint_error = norm(constraints);

    %%% Compute new free-variable vector if error is not converged
    if (constraint_error > constraint_error_tol)
        F_new = F_new - stepSize * DF_mat'*((DF_mat*(DF_mat'))\constraints);
    end

end

%%% Make sure the result is good, then store as new PO_1
if counter == ms_iterMax
    warning('Initial guess failed to converge')
    return
elseif constraint_error <= constraint_error_tol
    PO_new = [F_new(1), y0, F_new(2:5)', F_new(end)]';
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


end % function