function PlotBoi3Norm()
%%% Description
%       Make the 3D plots you've always dreamed about ... labels the three
%       axes in the normalized system
%       
% --------------------------------------------------------------
%%% Inputs
%       xAxisLabel     - [string] x-axis label 
%       yAxisLabel     - [string] y-axis label 
%       zAxisLabel     - [string] z-axis label 
%       fs             - [1x1] font size
%       interpLanguage - [string] string interpretation language .. use
%                        'LaTex'
% --------------------------------------------------------------
%%% Outputs
%       N/A
% --------------------------------------------------------------
%%% Author
%       Luke Bury, luke.bury@colorado.edu
% ===============================================================
%%% Auto-set parameters 
xAxisLabel = '$x_n$';
yAxisLabel = '$y_n$';
zAxisLabel = '$z_n$';
fs = 20;
interpLanguage = 'LaTex';

%%% Set labels, font size, and interpreter
xlabel(xAxisLabel,'FontName','Times New Roman','Fontsize',fs,'Interpreter',interpLanguage)
ylabel(yAxisLabel,'FontName','Times New Roman','Fontsize',fs,'Interpreter',interpLanguage)
zlabel(zAxisLabel,'FontName','Times New Roman','Fontsize',fs,'Interpreter',interpLanguage)

%%% Make background of plot white
set(gcf,'color','white')

%%% Turn on the grid
grid on

end