function [value, isterminal, direction] = event_xyRegionEscape_CR3Bn(t,X,prms)
%%% Designed for standard normalized CR3BP
%%% Event function watching for when "value" = 0
%%% Inputs:
%          t - normalized time vector
%          X - initial state [6x1]
%          prms - (xRegion, yRegion)


[out1] = isNegative(abs(X(1)) - prms.xRegion);

[out2] = isNegative(abs(X(2)) - prms.yRegion);

value = out1 + out2 - 2;
isterminal = 1;
direction = 0;
end