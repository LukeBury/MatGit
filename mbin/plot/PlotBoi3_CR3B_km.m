function PlotBoi3_CR3B_km(fs)
%%% Description
%       Make the 3D plots you've always dreamed about
%       
% --------------------------------------------------------------
%%% Inputs
%       fs - font size
% --------------------------------------------------------------
%%% Outputs
%       N/A
% --------------------------------------------------------------
%%% Author
%       Luke Bury, luke.bury@colorado.edu
% ===============================================================
%%% Set labels, font size, and interpreter  
xlabel('$x$ ($km$)','FontName','Times New Roman','Fontsize',fs,'Interpreter','LaTex')
ylabel('$y$ ($km$)','FontName','Times New Roman','Fontsize',fs,'Interpreter','LaTex')
zlabel('$z$ ($km$)','FontName','Times New Roman','Fontsize',fs,'Interpreter','LaTex')


%%% Make background of plot white
set(gcf,'color','white')

%%% Turn on the grid
grid on

end