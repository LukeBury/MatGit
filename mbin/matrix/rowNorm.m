function [ colVec ] = rowNorm( mat )
%%% Use:
%   - Takes the norm of each row of matrix and condenses to column vector
%
%%% Inputs:
%   1) Matrix [nxm]
%
%%% Outputs:
%   1) Column vector [nx1]
% ========================================================================
%%% Preallocate output vector
colVec = zeros(size(mat,1),1);

%%% Loop through rows, take norm, store in output vector
for row = 1:size(mat,1)
    n = norm(mat(row,:));
    colVec(row) = n;
end

end

