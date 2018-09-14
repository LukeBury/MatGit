clear
clc
close all
hold all

gross1 = [1:10:18450]/1000;
gross2 = [18451:10:74900]/1000;
gross3 = [74901:10:151200]/1000;
gross4 = [151201:10:230450]/1000;
gross5 = [230451:10:411500]/1000;
gross6 = [411501:10:464850]/1000;
gross7 = [464851:10:600000]/1000;


%%% Gross
plot(gross1,gross1,'--r','linewidth',2)
plot(gross2,gross2,'--r','linewidth',2)
plot(gross3,gross3,'--r','linewidth',2)
plot(gross4,gross4,'--r','linewidth',2)
plot(gross5,gross5,'--r','linewidth',2)
plot(gross6,gross6,'--r','linewidth',2)
plot(gross7,gross7,'--r','linewidth',2)


%%% After Tax
plot(gross1,(1-0.1)*gross1,'b','linewidth',2)
plot(gross2,(1-0.15)*gross2,'b','linewidth',2)
plot(gross3,(1-0.25)*gross3,'b','linewidth',2)
plot(gross4,(1-0.28)*gross4,'b','linewidth',2)
plot(gross5,(1-0.33)*gross5,'b','linewidth',2)
plot(gross6,(1-0.35)*gross6,'b','linewidth',2)
plot(gross7,(1-0.396)*gross7,'b','linewidth',2)

%%% Bracket Lines
plot([18.45 18.45],get(gca,'ylim'),'--k')
plot([74.9 74.9],get(gca,'ylim'),'--k')
plot([151.2 151.2],get(gca,'ylim'),'--k')
plot([230.45 230.45],get(gca,'ylim'),'--k')
plot([411.5 411.5],get(gca,'ylim'),'--k')
plot([464.85 464.85],get(gca,'ylim'),'--k')

PlotBoi('Gross Income (k)','Income After Tax (k)')