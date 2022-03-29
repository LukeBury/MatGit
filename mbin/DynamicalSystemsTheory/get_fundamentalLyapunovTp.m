function [Tp_fundamentalLyap] = get_fundamentalLyapunovTp(A_mat_Lp)
%%% Description
% Returns the fundamental (minimum) time period for Lyapunov orbits at
% either L1 or L2 (whichever point A is evaluated at)
%       
% ------------------------------------------------------------------------
%%% Inputs
% A_mat - [6x6] Dynamics/State jacobian evaluated at L1 or L2
% ------------------------------------------------------------------------
%%% Outputs
% Tp_fundamentalLyap - [scalar] Fundamental (minimum) time period of
%                       Lyapunov family at either L1 or L2
% ------------------------------------------------------------------------
% Created: 4/14/21
% Author : Luke Bury, luke.bury@colorado.edu
% ========================================================================
%%% Get eigenvalues and eigenvectors of A matrix
[eigenvectors, eigenvalues] = eig(A_mat_Lp);

%%% Find the eigenvector that only bumps 'x' and 'y-dot', as this
%%% corresponds to Lyapunovs
for kk = 1:6
    absRealEigenvector = abs(real(eigenvectors(:,kk)));
    absRealEigenvector = absRealEigenvector ./ norm(absRealEigenvector);
    
    if (absRealEigenvector(1) >= 1e-3) && (absRealEigenvector(5) >= 1e-3)
        if sum(absRealEigenvector([2,3,4,6])) < 1e-12
            Tp_fundamentalLyap = 2*pi / abs(imag(eigenvalues(kk,kk)));
            return
        else
            continue
        end
    else
        continue
    end
end

warning('Fundamental Lyapunov time period not found')

end % function