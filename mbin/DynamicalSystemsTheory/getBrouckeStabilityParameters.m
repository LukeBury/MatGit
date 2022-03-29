function [alpha, beta] = getBrouckeStabilityParameters(eigenValues_row, monodromy)
%%% Description
% Given eigenvalues of monodromy matrix for a periodic orbit, returns the
% two parameters, alpha and beta, necessary for use of the Broucke
% stabiltiy diagram
%       
% ------------------------------------------------------------------------
%%% Inputs
% eigenValues - [1x6] Diagonal matrix of eigenvalues from monodromy 
%               matrices 
% ------------------------------------------------------------------------
%%% Outputs
% alpha - [scalar] Broucke stability parameter "alpha"
% beta - [scalar] Broucke stability parameter "beta"
% ------------------------------------------------------------------------
% Created: 12/08/20
% Author : Luke Bury, luke.bury@colorado.edu
% ========================================================================
% -------------------------------------------------
%%% Grab eigenvalue pairs ([1,2] and [3,4]) and compute alpha and beta
%%% parameters
% -------------------------------------------------
% [EV_compConj_1, EV_compConj_2, EV_inverse_1, EV_inverse_2, ~, ~] = getMonodromyEigenvalueIdentification(eigenValues_row);
% 
% alpha = -(EV_compConj_1 + EV_compConj_2 + EV_inverse_1 + EV_inverse_2);
% beta  = 0.5 * (alpha^2 - (EV_compConj_1^2 + 1/(EV_compConj_1^2) + EV_compConj_2^2 + 1/(EV_compConj_2^2)));
% % alpha = -(real(EV_compConj_1) + real(EV_compConj_2) + real(EV_inverse_1) + real(EV_inverse_2));
% % beta  = 0.5 * (alpha^2 - (real(EV_compConj_1)^2 + 1/(real(EV_compConj_1)^2) + real(EV_compConj_2)^2 + 1/(real(EV_compConj_2)^2)));


alpha = 2 - trace(monodromy);
beta = 0.5 * (alpha^2 + 2 - trace(monodromy^2));


end % function