function [S1_vec, S2_vec] = getStabilityIndices(eigenValues)
%%% Description
%       Calculate the stability indicies for 6-dimensional sets of
%       eigenvalues of the monodromy matrix (stm_tf_t0)
%       
% ------------------------------------------------------------------------
%%% Inputs
%       eigenValues - [nx6] Rows of eigenvalues from monodromy matrices
% ------------------------------------------------------------------------
%%% Outputs
%       S1_vec - [nx1] 1st stability index (s1 = real(lamda1)+real(lamda2)) 
%       S2_vec - [nx1] 2nd stability index (s2 = real(lamda3)+real(lamda4)) 
% ------------------------------------------------------------------------
% Created: 06/20/19
% Author : Luke Bury, luke.bury@colorado.edu
% ========================================================================
% -------------------------------------------------
%%% Check input dimension
% -------------------------------------------------
if ~isequal(size(eigenValues,2),6)
    %%% Try correcting
    eigenValues = eigenValues';
    
    if ~isequal(size(eigenValues,2),6)
        warning('Error with input dimension')
        return
    end
end

% -------------------------------------------------
%%% Loop through input matrix of eigenvalues and find stability indices
% -------------------------------------------------
%%% Number of 1x6 eigenvalue vectors
n = size(eigenValues,1);

%%% Preallocate output
S1_vec = NaN(n,1);
S2_vec = NaN(n,1);

%%% Loop through eigenvalue sets
for kk = 1:n
    %%% Acquire the six eigenvalues identified by pairs
    [EV_compConj_1, EV_compConj_2, EV_inverse_1, EV_inverse_2, ~, ~] = getMonodromyEigenvalueIdentification(eigenValues(kk,:));
            
    %%% Calculate stability indices and store output
    S1_vec(kk) = real(EV_compConj_1) + real(EV_compConj_2);
    S2_vec(kk) = real(EV_inverse_1) + real(EV_inverse_2);
end % kk = 1:n

end % function



