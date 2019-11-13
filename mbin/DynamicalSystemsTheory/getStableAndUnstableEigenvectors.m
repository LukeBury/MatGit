function [eVec_stable, eVec_unstable] = getStableAndUnstableEigenvectors(monodromy_mat)
%%% Description
%       Return the stable and unstable eigenvectors associated with the
%       monodromy matrix of a periodic orbit
%       
% ------------------------------------------------------------------------
%%% Inputs
%       monodromy_mat - [6x6] Monodromy matrix for a 6-state periodic orbit
% ------------------------------------------------------------------------
%%% Outputs
%       eVec_stable   - [6x1] Stable eigenvector
%       eVec_unstable - [6x1] Unstable eigenvector
% ------------------------------------------------------------------------
% Created: 09/13/19
% Author : Luke Bury, luke.bury@colorado.edu
% ========================================================================
% -------------------------------------------------
%%% Check input dimension
% -------------------------------------------------
if isequal(size(monodromy_mat),[6,6]) == 0
    warning('Error with input dimension')
    return
end

% -------------------------------------------------
%%% Determine stable and unstable eigenvectors
% -------------------------------------------------
%%% Get the eigenvalues and eigenvectors of the monodromy matrix
[eigenVectors, eigenValues] = eig(monodromy_mat);

%%% Take diagonal of eigenvalues to make it 1D rather than 2D
diag_eigenValues_new = diag(eigenValues);

%%% Find index of unstable eigenvalue (largest real part)
unstableEigenvalue_index = find(abs(real(diag_eigenValues_new)) == max(abs(real(diag_eigenValues_new))));

%%% If the largest real part comes in a complex conjugate pair, such that
%%% the maximum real part shows up more than once, break
if length(unstableEigenvalue_index) == 2
    %%% ---------------------
    %%% Try removing the two eigenvalues closest to one since those are
    %%% periodic, and see if there is now an unstable pair
    %%% ---------------------
    %%% Grab real parts of current eigenvalues
    realParts = real(diag_eigenValues_new);
    
    %%% Remove the complex pair closest to 1 by differencing the real
    %%% values by one and removing the smallest remainders
    realPartsDifferencedWithOne = 1-abs(realParts);
    logicalIndicesOfValuesClosestToOne = (realPartsDifferencedWithOne == min(realPartsDifferencedWithOne));
    logicalIndicesOfRemainingValues = ~logicalIndicesOfValuesClosestToOne;
    diag_eigenValues_new = diag_eigenValues_new(logicalIndicesOfRemainingValues);
    
    %%% Try to find an unstable eigenvalue from the 4 remaining values
    unstableEigenvalue_index = find(abs(real(diag_eigenValues_new)) == max(abs(real(diag_eigenValues_new))));
    
    %%% If the largest real part still belons to a complex conjugate pair,
    %%% then there are no manifolds and the orbit is probably just stable
    if length(unstableEigenvalue_index) == 2
%         warning('This orbit doesn''t have stable/unstable manifolds')
        eVec_stable   = [];
        eVec_unstable = [];
        return
    end

elseif length(unstableEigenvalue_index) > 2
    warning('This shouldn''t be happening...')
    return
end

%%% Use index to determine the unstable eigenvalue
unstableEigenvalue_real = real(diag_eigenValues_new(unstableEigenvalue_index));

%%% Create a temporary vector to help identify the stable eigenvalue, since
%%% the two have the relationship of real(EV1_unstable) = 1/real(EV_stable)
tempVec = abs(real(unstableEigenvalue_real) - 1./real(diag_eigenValues_new));

%%% Find index of stable eigenvalue
stableEigenvalue_index = find(tempVec == min(tempVec));

%%% Output the corresponding eigenvectors
eVec_stable   = eigenVectors(:, stableEigenvalue_index);
eVec_unstable = eigenVectors(:, unstableEigenvalue_index);
end % function