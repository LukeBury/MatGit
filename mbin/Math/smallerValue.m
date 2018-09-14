function [ valSmall ] = smallerValue( val1, val2 )
%%% Calculates smaller of two values. If equal, returns value
%%% Inputs:
%       1) val1 - first value for comparison
%       2) val2 - second value for comparison
%%% Outputs: 
%       1) valSmall - smaller of the two values
%=========================================================================
diff = val1 - val2;

if diff < 0
    valSmall = val1;
elseif diff > 0
    valSmall = val2;
elseif diff == 0
    valSmall = val1;
end

end

