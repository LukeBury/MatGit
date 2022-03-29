function [ contourPointsOut ] = getContourPoints( X, Y, Z, contourLevel )
%%% Returns points of certain level of a given contour map
%%% Inputs:
%       1) X - X components of meshgrid
%       2) Y - Y components of meshgrid
%       3) Z - Z components of meshgrid
%       4) contourLevel - value pertaining to level of contour desired
%
%%% Outputs: 
%       1) contourPoints - [2xn] matrix of points defining contour shape
%=========================================================================

%%% Create contour points
figure(18373)
[contourPointsOut,handle] = contour(X,Y,Z,[contourLevel, contourLevel],'Visible','off');

%%% Clean up weirdness
contourPointsOut(contourPointsOut>1) = NaN; % fixing weird columns that contour returns
contourPointsOut = contourPointsOut(:,all(~isnan(contourPointsOut))); % Cutting out 'NaN' columns

%%% Closing plot that inevitably appears
close(figure(18373))
end

