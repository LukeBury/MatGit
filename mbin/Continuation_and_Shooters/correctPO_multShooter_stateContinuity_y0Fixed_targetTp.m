function [PO_new, counter, constraint_error] = correctPO_multShooter_stateContinuity_y0Fixed_targetTp(PO_guess, N_Nodes, integrator, integrator_options, prms, constraint_error_tol, ms_iterMax, stepSize, Tdes)
%%% Description
% Multiple shooting algorithm to correct a periodic orbit for state
% continuity with a fixed y0 and a target, but not fixed, Tp
%       
% ------------------------------------------------------------------------
%%% Inputs
% PO_guess             - [7x1] Vector of guess state and time period
% N_Nodes              - [scalar] number of nodes for multiple shooter
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
% Tdes                 - [1x1] Desired time period

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
if N_Nodes == 1
    X_nodes = PO_guess(1:6)';
else
    stm0_colvec = reshape(eye(6),36,1);
    [~, X_nodes] = get_nodes([PO_guess(1:6); stm0_colvec], [0, PO_guess(7)], N_Nodes+1, integrator, integrator_options, prms);
    X_nodes = X_nodes(1:N_Nodes,1:6);
end

%%% Create free-variable vector
F_new   = reshape(X_nodes',6*N_Nodes,1);
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
    [DF_mat, constraints] = multShooter_stateContinuity_y0Fixed_targetTp(N_Nodes, F_new, integrator, integrator_options, prms, PO_guess(2), Tdes);

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
    PO_new = [F_new(1), 0, F_new(2:5)', F_new(end)]';
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