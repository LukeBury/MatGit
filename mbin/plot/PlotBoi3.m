function PlotBoi3(xAxisLabel,yAxisLabel,zAxisLabel,fs,interpLanguage)
%%% Description
%       Make the 3D plots you've always dreamed about
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
%%% Set labels, font size, and interpreter
if nargin == 4
    xlabel(xAxisLabel,'FontName','Times New Roman','Fontsize',fs)
    ylabel(yAxisLabel,'FontName','Times New Roman','Fontsize',fs)
    zlabel(zAxisLabel,'FontName','Times New Roman','Fontsize',fs)
elseif nargin == 5
    xlabel(xAxisLabel,'FontName','Times New Roman','Fontsize',fs,'Interpreter',interpLanguage)
    ylabel(yAxisLabel,'FontName','Times New Roman','Fontsize',fs,'Interpreter',interpLanguage)
    zlabel(zAxisLabel,'FontName','Times New Roman','Fontsize',fs,'Interpreter',interpLanguage)
end

%%% Make background of plot white
set(gcf,'color','white')

%%% Turn on the grid
grid on

end