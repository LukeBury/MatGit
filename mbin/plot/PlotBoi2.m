function PlotBoi2(A,B,s,interpLanguage)
%%%%    
%%%%    Use: PlotBoi('x-axis label','y-axis label',FontSize)
%%%%
%%%%    Created by: The One and Only Luke Bury
%%%%
if nargin == 3
    xlabel(A,'FontName','Times New Roman','Fontsize',s)
    ylabel(B,'FontName','Times New Roman','Fontsize',s)
elseif nargin == 4
    xlabel(A,'FontName','Times New Roman','Fontsize',s,'Interpreter',interpLanguage)
    ylabel(B,'FontName','Times New Roman','Fontsize',s,'Interpreter',interpLanguage)
end
set(gcf,'color','white')
grid on
end