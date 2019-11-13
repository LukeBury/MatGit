function [output] = isNegative(input)
%%% Description
%       Returns a zero if the input is a negative number
%       
% ------------------------------------------------------------------------
%%% Inputs
%       input - [1x1] Input number
% ------------------------------------------------------------------------
%%% Outputs
%       output - [1x1] Output number
% ------------------------------------------------------------------------
% Created: 09/13/10
% Author : Luke Bury, luke.bury@colorado.edu
% ========================================================================
if input < 0
    output = 1;
else
    output = 0;
end
end % function