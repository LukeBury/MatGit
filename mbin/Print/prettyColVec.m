function prettyColVec(vector)
%%% Description
%       Prints an input vector to the command window as a column vector in
%       a convenient format for copy-pasting-using
%       
% ------------------------------------------------------------------------
%%% Inputs
%       vector - [1xn] or [nx1] input vector
% ------------------------------------------------------------------------
%%% Outputs
%
% ------------------------------------------------------------------------
% Created: 07/10/19
% Author : Luke Bury, luke.bury@colorado.edu
% ========================================================================
% -------------------------------------------------
%%% Print column vector
% -------------------------------------------------
fprintf('[')
for kk = 1:length(vector)-1
    if kk == 1
        fprintf('%1.15f;\n',vector(kk))
    else
        fprintf(' %1.15f;\n',vector(kk))
    end

end
fprintf(' %1.15f];\n\n',vector(end))

end % function 