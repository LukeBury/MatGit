function [ colorMatrix ] = colorScale(colors,n )
%%% colorScale
%   Description: Takes in two colors and outputs a matrix of intermediate
%   colors useful for showing progression in plots
%
%%% Inputs
%   1) colors [2x3], two rows, each a rgb color code
%   2) n [1x1], desired number of color steps for output
%
%%% Outputs
%   1) colorMatrix [nx3], matrix where each row is a color
% =========================================================================

%%% Assigning intial and final colors by index
c11 = colors(1,1); c12 = colors(1,2); c13 = colors(1,3);
c21 = colors(2,1); c22 = colors(2,2); c23 = colors(2,3);

for cc = 1:n-1
    colorMatrix(cc,1) = c11 + (cc-1)*(c21-c11)/n;
    colorMatrix(cc,2) = c12 + (cc-1)*(c22-c12)/n;
    colorMatrix(cc,3) = c13 + (cc-1)*(c23-c13)/n;
end
colorMatrix(n,:) = [c21, c22, c23];

end

