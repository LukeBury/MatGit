clc
close all

%% ------------------------------------------------------------------------
%%%                                    #1
%%%------------------------------------------------------------------------
close all
N_mill = [145601;7785;7365;4911;3861;2421];
N_sand = [516210;13008;9736;6023;4285;3837];
stress = [29700.90;47521.45;50491.54;53461.63;56431.72;59401.81];

sA = (max(stress)-min(stress))/2 % = 1.4850e+04 psi




%%% 1c
pm = polyfit(log(N_mill),log(stress),1);
ps = polyfit(log(N_sand),log(stress),1);

m_m = -pm(1); % 0.1721
m_c = pm(2);
m_c = exp(m_c) % 2.2962e+05

sig_m = m_c./(N_mill.^m_m);

s_m = -ps(1); % 0.1362
s_c = ps(2);
s_c = exp(s_c) % 1.7683e+05

figure
hold all
plot(log(N_mill),log(stress),'o')
plot(log(N_mill),pm(1)*log(N_mill)+pm(2))
legend('Not Blasted')
PlotBoi('log(N)','log(Stress)')

figure
hold all
plot(log(N_sand),log(stress),'x')
plot(log(N_mill),ps(1)*log(N_mill)+ps(2))
legend('Blasted')
PlotBoi('log(N)','log(Stress)')


% figure
% hold all
% plot(N_mill,stress,'o','linewidth',2,'markersize',9)
% plot(N_mill,sig_m-(sig_m(end)-stress(end)))
% % legend('Not Blasted','Blasted')
% PlotBoi('Number of Cycles','Stress, psi')
x = [2000:100:6e05];

figure
hold all
plot(N_mill,stress,'o','linewidth',2,'markersize',9)
plot(N_sand,stress,'x','linewidth',2,'markersize',9)
plot(x, m_c./(x.^m_m))
plot(x, s_c./(x.^s_m))
legend('Not Blasted','Blasted')
PlotBoi('Number of Cycles','Stress, psi')


%% ------------------------------------------------------------------------
%%%                                    #3
%%%------------------------------------------------------------------------

N = [0;1540;2740;3440;3900;4190;4300;4330];
a = [1.25;1.375;1.5;1.625;1.75;1.875;2;2.125];

aww = a./3.25;
qaw = [1.875;1.789;1.712;1.634;1.567;1.5;1.443;1.401];
%%% x axis = N
figure
hold all
plot(N,a,'-o','linewidth',2,'markersize',7)
PlotBoi('Number of Cycles','Crack Length, inches')

p3 = polyfit(N,a,2)
x = [0:4330];
% plot(x,p3(1)*x.^2+p3(2).*x+p3(3),'linewidth',2)
% legend('Data','Polynomial Fit')
% %%% x axis = a
% figure
% hold all
% plot(a,N,'-o','linewidth',2,'markersize',7)
% PlotBoi('Crack Length, inches','Number of Cycles')
% 
% p3 = polyfit(a,N,2)
% x = [1.25:.01:2.125];
% plot(x,p3(1)*x.^2+p3(2).*x+p3(3))


%%%%%%%%%%%  
%%%%%%%%%%%                       (c)
%%%%%%%%%%%
%%

aw = [.2;.3;.4;.5;.6;.7;.75;.8];
Qaw= [2.07;2.09;1.836;1.634;1.46;1.351;1.3325;1.314];

dP = (1-.1)*2222;
B  = 0.5;
W  = 3.25;

dK = (dP/(B*sqrt(pi*W)))*(2+aw).*Qaw./((1-aw).^1.5);
dkk= (dP/(B*sqrt(pi*W)))*(2+aww).*qaw./((1-aww).^1.5);

figure
plot(a,dK,'-','linewidth',2,'markersize',7)
PlotBoi('Crack Length, inches','Stress Intensity Factor')


dadN = 2*p3(1)*N + p3(2);


%%% e
figure
plot(log(dK),log(dadN),'-o')
PlotBoi('log(dK)','log(da/dN)')




