function [S1_vec, S2_vec] = getStabilityIndices(eigenValues)
%%% Description
%       Calculate the stability indicies for 6-dimensional sets of
%       eigenvalues of the monodromy matrix (stm_tf_t0)
%       
% ------------------------------------------------------------------------
%%% Inputs
%       in1 - [nx6] n sets of eigenvectors from monodromy matrices 
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
if isequal(size(eigenValues,2),6) == 0
    %%% Try correcting
    eigenValues = eigenValues';
    
    if isequal(size(eigenValues,2),6) == 0
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
    %%% Set current eigenvalues
    remainingEigenvalues = eigenValues(kk,:);
    
    % -------------------------------------------------
    %%% Find and remove the two cases of lambda = 1 (note, it's never equal
    %%% to 1 exactly(
    % -------------------------------------------------
    for jj = 1:2
        %%% Difference of all EVs with 1
        tempVec = abs(remainingEigenvalues - 1);
        
        %%% Only keep those that aren't the minimum
        remainingEigenvalues = remainingEigenvalues(~(tempVec == min(tempVec)));
        
        %%% Set this in case the first two minima are equal and both
        %%% get removed
        if isequal(size(remainingEigenvalues),[1,4])
            break
        end
    end
    
    % -------------------------------------------------
    %%% Find EV1-EV4 for current eigenvalue set
    % -------------------------------------------------
    %%% Find lambda1 (we'll define it as the EV with the largest real
    %%% component)
    maxRealEV = remainingEigenvalues(real(remainingEigenvalues) == max(real(remainingEigenvalues)));
    
    %%% Check if there were two maximums - it's possible that the largest 
    %%% real component belonged to a complex conjugate pair who have the 
    %%% same real component
    if isequal(size(maxRealEV),[1,1]) == 1 %%% Single Value - EV2 must follow real(EV1) = 1/real(EV2)
        %%% Set EV1
        EV1 = maxRealEV;
        
        %%% Remove EV1 from remainingEigenvalues
        remainingEigenvalues = remainingEigenvalues(~(remainingEigenvalues == EV1));
        
        %%% Finding EV2 from real(EV1) = 1/real(EV2)
        tempVec = abs(real(EV1) - 1./real(remainingEigenvalues));
        EV2 = remainingEigenvalues(tempVec == min(tempVec));
        
        %%% Remove EV2 from remainingEigenvalues if it's [1x1]
        remainingEigenvalues = remainingEigenvalues(~(remainingEigenvalues == EV2));
        
    elseif isequal(size(maxRealEV),[1,2]) == 1 %%% EV1 has complex conjugate 
        %%% Pick out the complex conjugate pair
        EV1 = maxRealEV(1);
        EV2 = maxRealEV(2);
        
        %%% Remove EV2 from remainingEigenvalues
        remainingEigenvalues = remainingEigenvalues(~(real(remainingEigenvalues) == max(real(remainingEigenvalues))));
    end
    
    %%% Only two EVs left - they must be a pair
    EV3 = remainingEigenvalues(1);
    EV4 = remainingEigenvalues(2);
        
    %%% Calculate stability indices and store output
    S1_vec(kk) = real(EV1) + real(EV2);
    S2_vec(kk) = real(EV3) + real(EV4);
end % kk = 1:n

end % function



