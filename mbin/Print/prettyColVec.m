function prettyColVec(vector, varargin)
%%% Description
%       Prints an input vector to the command window as a column vector in
%       a convenient format for copy-pasting-using
%       
% ------------------------------------------------------------------------
%%% Inputs
%       vector - [1xn] or [nx1] input vector
%       varargin - {nx1} option for switching to row vector
%              ex: varargin(1) = 'row'
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
%%% Loop through the variable-length input argument list and look for axis
%%% shadow specfications
for arg_index = 1:length(varargin)
    %%% Setting current argument
    arg = varargin{arg_index};
    
    %%% If the current argument is 'row', plot the x-axis shadow
    if isequal(lower(arg), 'row')
        
        fprintf('[')
        for kk = 1:length(vector)-1
            if kk == 1
                fprintf('%1.16f,',vector(kk))
            else
                fprintf(' %1.16f,',vector(kk))
            end

        end
        fprintf(' %1.16f];\n\n',vector(end))

        return
    end
end


fprintf('[')
for kk = 1:length(vector)-1
    if kk == 1
        fprintf('%1.16f;\n',vector(kk))
    else
        fprintf(' %1.16f;\n',vector(kk))
    end

end
fprintf(' %1.16f];\n\n',vector(end))

end % function 