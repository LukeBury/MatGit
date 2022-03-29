function [F_new, constraint_error, DF_mat, iteration_counter] = pseudoArclengthContinuation_multShoot_stateContinuity(error_tol, iterMax, n_Nodes, F_new, dynamicsWithSTM, options, prms, familyIndex, stepSize, nullVecDF, ds_PO)
%%% Description
%       This script runs a while-loop to converge on the next member of a
%       family of a free-variable state. The algorithm in
%       place is a pseudo-arclength continuation algorithm. The
%       continuation algorithms generates a guess, and a multiple-shooting
%       algorithm converges the guess to a state continuity within a
%       provided error tolerance
%       
% ------------------------------------------------------------------------
%%% Inputs
%       error_tol       - [1x1] Tolerance for magnitude of the constraint 
%                           error vector
%       iterMax         - [1x1] Maximum number of iterations  for while
%                           loop when trying to converge on new family
%                           member
%       n_Nodes         - [1x1] Number of nodes in multiple shooter 
%       F_new           - [nx1] Free-variable vector (stacked states of
%                           each node with total time at end (so 
%                           n = n_nodes*6+1, generally))
%       dynamicsWithSTM - [@FunctionHandle] Function handle to dynamics 
%                           for ode integrator to use. Needs to include the
%                           dynamics for the state transition matrix
%       options         - [struct] options struct for ode integration 
%       prms            - [struct] parameters needed for integration (here,
%                           prms.u for the mass ratio)
%       familyIndex     - [1x1] Index of current family member from the
%                           surrounding loop
%       stepSize        - [1x1] Artificial scalar for increasing or
%                           decreasing the step size for guessing at the
%                           state of the 2nd family member. 1 has no
%                           effect, 2 and 3 are somewhat common
%       nullVecDF       - [nx1] Null space of the previous DF matrix (this
%                           vector can just be initialized as zeros if the
%                           familyIndex is 1 so there's no previous DF
%       ds_PO           - [1x1] Essentially, a step size parameter for how
%                           finding the next family member (probably
%                           somewhere between 1e-9 and 1e-2)
% ------------------------------------------------------------------------
%%% Outputs
%       F_new             - [nx1] New free-variable vector
%       constraint_error  - [1x1] Norm of the augmented constraint vector.
%                             Based on the value of this error, the
%                             parameters will either be adjusted and 
%                             re-run, or the F_new will be stored and the
%                             next family member will be found
%       DF_mat            - [(n_nodes*6)xn] Matrix of the partials of the
%                             constraints with respect to the free 
%                             variables
%       iteration_counter - [1x1] Total number of iterations it took for
%                             contraint_error to converge to the error_tol
% ------------------------------------------------------------------------
% Created: 09/06/19
% Author : Luke Bury, luke.bury@colorado.edu
% ========================================================================
% -------------------------------------------------
%%% Pre-While-Loop work
% -------------------------------------------------
%%% Set while-loop/convergence variables
iteration_counter = 0;
constraint_error  = 100;

%%% Set F_old to compare against updated F_new
F_old = F_new;
   
% --------------------------
% Enter while loop - attempt to converge on next PO in family
% --------------------------
while (constraint_error > error_tol) && (iteration_counter < iterMax)

    %%% Count iteration
    iteration_counter = iteration_counter + 1;

    % --------------------------
    % Using multiple shooting - loop through nodes, integrate, and 
    % populate constraint vector and create DF matrix
    % --------------------------
    [DF_mat, constraints] = multShooter_stateContinuity(n_Nodes, F_new, dynamicsWithSTM, options, prms);

    % --------------------------
    % If familyIndex == 1, do a simplified version of things to find next
    % F... This part is really only here to make sure that your first
    % solution is within the error bound (in case your PO_IC wasn't)
    % --------------------------
    if familyIndex == 1
        %%% Compute error
        constraint_error = norm(constraints);

        %%% Compute new free-variable vector if error is not converged
        if (constraint_error > error_tol)
            F_new = F_new - stepSize * DF_mat'*((DF_mat*(DF_mat'))\constraints);
        end

    % --------------------------
    % If familyIndex > 1, do full process for finding next F
    % --------------------------
    elseif familyIndex > 1
        %%% Compute z, which perturbs the constraints to continue the
        %%% family rather than converge on the same answer
        z = (F_new - F_old)'*nullVecDF - ds_PO;

        %%% Compute augemented constraint vector
        constraints_aug = [constraints; z];

        %%% Compute DH matrix
        DH = [DF_mat; nullVecDF'];

        %%% Compute error
        constraint_error = norm(constraints_aug);

        %%% Compute new free-variable vector if error is not converged
        if (constraint_error > error_tol)
            F_new = F_new - (DH\constraints_aug);
        end

    end

end % while constraint_error > error_tol







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