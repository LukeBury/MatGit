function PlotBoi3(A,B,C,s,interpLanguage)
%%%%    
%%%%    Use: PlotBoi('x-axis label','y-axis label','z-axis label',fontsize)
%%%%
%%%%    Created by: The One and Only Luke Bury
%%%%
if nargin == 4
    xlabel(A,'FontName','Times New Roman','Fontsize',s)
    ylabel(B,'FontName','Times New Roman','Fontsize',s)
    zlabel(C,'FontName','Times New Roman','Fontsize',s)
    set(gcf,'color','white')
    grid on
elseif nargin == 5
    xlabel(A,'FontName','Times New Roman','Fontsize',s,'Interpreter',interpLanguage)
    ylabel(B,'FontName','Times New Roman','Fontsize',s,'Interpreter',interpLanguage)
    zlabel(C,'FontName','Times New Roman','Fontsize',s,'Interpreter',interpLanguage)
    set(gcf,'color','white')
    grid on
end
end