function [ rowVec ] = colNorm( mat )
%%% Use:
%   Takes the norm of each column of matrix and condenses to row vector
%
%%% Inputs:
%   1) Matrix [nxm]
%
%%% Outputs:
%   1) Row vector [1xm]
% ========================================================================
%%% Preallocate output vector
rowVec = zeros(1,size(mat,2));

%%% Loop through columns, take norm, store in output vector
for col = 1:size(mat,2)
    n = norm(mat(:,col));
    rowVec(col) = n;
end

end

