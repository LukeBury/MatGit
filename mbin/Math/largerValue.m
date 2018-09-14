function [ valLarge ] = largerValue( val1, val2 )
%%% Calculates larger of two values. If equal, returns value
%%% Inputs:
%       1) val1 - first value for comparison
%       2) val2 - second value for comparison
%%% Outputs: 
%       1) valSmall - larger of the two values
%=========================================================================
diff = val1 - val2;

if diff < 0
    valLarge = val2;
elseif diff > 0
    valLarge = val1;
elseif diff == 0
    valLarge = val1;
end

end

