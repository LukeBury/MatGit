clear
clc
close all

%% ------------------------------------------------------------------------
%%%                           Acquiring Data
%%%------------------------------------------------------------------------

Data  = xlsread('Data');

%%%                              Parsing

T = wrev(Data(1:10,3)); % Deg Celsius
E = wrev(Data(1:10,4)); % Ft*lbs force
bo = Data(1:10,5);
ao = Data(1:10,6);
bc = Data(1:10,7);
ac = Data(1:10,8);

% perc = (bo.*bc)./(ao.*ac)
%% ------------------------------------------------------------------------
%%%                           Prob 1
%%%------------------------------------------------------------------------
close all
hold all
plot(T,E,'rx','markersize',8,'linewidth',2)


p1 = polyfit(T(1:6),E(1:6),1);
p3 = polyfit(T(7:10),E(7:10),1);


plot(T(1:6),p1(1)*T(1:6)+p1(2),'b');
plot(T(6:7),E(6:7),'b-');
plot(T(7:10),p3(1)*T(7:10)+p3(2),'b');

PlotBoi('Temperature, C','Energy, ft*lbs')
%(b) 37.5 deg = transition temperature








%% ------------------------------------------------------------------------
%%%                           Prob 3
%%%------------------------------------------------------------------------













