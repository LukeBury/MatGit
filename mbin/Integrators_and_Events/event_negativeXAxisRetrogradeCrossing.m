function [value, isterminal, direction] = event_negativeXAxisRetrogradeCrossing(t,X,prms)
%%% Event function watching for when "value" = 0
%%% Inputs:
%          t - normalized time vector
%          X - initial state [6x1]
%          prms - (u, R2_n)

value = X(2);
isterminal = 0; 
direction  = 1; 
end




