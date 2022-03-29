function [PO_new, counter, constraint_error] = correctPO_mS_sC_y0xd0TpFixed(PO_guess, n_Nodes, integrator, integrator_options, prms, constraint_error_tol, ms_iterMax, stepSize, y0, xd0, Tp)
%%% Description
% Multiple shooting algorithm to correct a periodic orbit for state
% continuity with a fixed time period and y0, and xd0 of first node
%       
% ------------------------------------------------------------------------
%%% Inputs
% PO_guess             - [6x1] Vector of guess state and time period
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
% y0                   - [scalar] Fixed y0 of first node
% xd0                  - [scalar] Fixed xd0 of first node
% Tp                   - [scalar] Fixed time period
% ------------------------------------------------------------------------
%%% Outputs
% PO_new           - [7x1] Corrected PO
% counter          - [scalar] Iterations required to converge
% constraint_error - [scalar] Norm of final constraint vector
% ------------------------------------------------------------------------
% Created: 1/15/21
% Author : Luke Bury, luke.bury@colorado.edu
% ========================================================================


%%% Acquire nodes
if n_Nodes == 1
    X_nodes = PO_guess(1:6)';
else
    stm0_colvec = reshape(eye(6),36,1);
    [~, X_nodes] = get_nodes([PO_guess(1:6); stm0_colvec], [0, Tp], n_Nodes+1, integrator, integrator_options, prms);
    X_nodes = X_nodes(1:n_Nodes,1:6);
end

%%% Create free-variable vector
F_new   = reshape(X_nodes',6*n_Nodes,1);
F_new   = [F_new(1); F_new(3); F_new(5:end)];

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
    [DF_mat, constraints] = multShooter_stateContinuity_y0xd0TpFixed(n_Nodes, F_new, integrator, integrator_options, prms, y0, xd0, Tp);
    
    
    %%% Compute error
    constraint_error = norm(constraints);

    %%% Compute new free-variable vector if error is not converged
    if (constraint_error > constraint_error_tol)
        warning('off','MATLAB:nearlySingularMatrix')
        F_new = F_new - stepSize * DF_mat'*((DF_mat*(DF_mat'))\constraints);
        warning('on','MATLAB:nearlySingularMatrix')
    end

end

%%% Make sure the result is good, then store as new PO_1
if counter == ms_iterMax
    warning('Initial guess failed to converge')
    return
elseif constraint_error <= constraint_error_tol
    PO_new = [F_new(1), y0, F_new(2), xd0, F_new(3), F_new(4), Tp]';
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