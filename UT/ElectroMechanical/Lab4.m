clc


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%                    Data Organization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% hold all
% plot(calV,calX)
% 
% %confirm accuracy of fit
% dist = poly(1)*calV+poly(2)
% plot(calV,dist)

d = @(V) poly(1)*V + poly(2);

V0 = [0.691992105	0.68887958	0.6963989	0.700756925	0.7071288	0.70554187	0.70065921	0.71580799	0.7048827	0.70871581	0.694396855	0.679150425	0.672058165];
V1 = [0.831408815	0.90174553	0.950280645	1.008569325	1.060241625	1.09628894	1.1096312	1.10969225	1.073999055	1.01494142	0.948266715	0.86395256	0.80717759];
V2 = [0.96715082	1.072070285	1.15725115	1.27359618	1.371118095	1.42789314	1.420702945	1.408752325	1.352331465	1.25168447	1.13861078	1.053796165	0.964221005];

X0 = [0 d(V0) 0];
X1 = [0 d(V1) 0];
X2 = [0 d(V2) 0];

D1 = X1 - X0;
D2 = X2 - X0;

% X locations of measurement points (+- .002)
x = [0 .106 .136 .157 .186 .216 .244 .261 .274 .308 .341 .374 .401 .429 .521];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%                    Graphing Deflections
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure % Measured Deflection
hold all
plot(x,-(X0-X0),'-o','linewidth',2,'markersize',6)
plot(x,-D1,'-o','linewidth',2,'markersize',6)
plot(x,-D2,'-o','linewidth',2,'markersize',6)
ylabel('Deflection Distance (m)')
xlabel('Distance along bar (m)')

legend('No Loading','Measured Light Loading','Measured Heavy Loading')


% Put in a vertical grid
plot([.136 .136],[0,.005],'--b')
plot([.106 .106],[0,.005],'--b')
plot([.157 .157],[0,.005],'--b')
plot([.186 .186],[0,.005],'--b')
plot([.216 .216],[0,.005],'--b')
plot([.244 .244],[0,.005],'--b')
plot([.261 .261],[0,.005],'--b')
plot([.274 .274],[0,.005],'--b')
plot([.308 .308],[0,.005],'--b')
plot([.341 .341],[0,.005],'--b')
plot([.374 .374],[0,.005],'--b')
plot([.401 .401],[0,.005],'--b')
plot([.429 .429],[0,.005],'--b')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%                    Loading Data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure %%% Loading Plot
hold all
plot(Loading_delta,Loading_force,'linewidth',2)
xlabel('Delta L (mm)')
ylabel('Loading Force (N)')
grid on

%% Marking the theoretical Pcr
size(Loading_delta)
vec = zeros(528,1);
vec(1:end) = 244.4849;
plot(Loading_delta,vec,'linewidth',2)
legend('Measured Loading vs Shortening','Theoretical Pcr')
% total end shortening
Loading_delta(end);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%                    Equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A1 = .5*D1(8) %m
A2 = .5*D2(8); %m

L = .521; %m

figure %%% Theoretical Deflection
hold all
plot(x,-(X0-X0),'-o','linewidth',2,'markersize',6)
plot(x,-D1,'-o','linewidth',2,'markersize',6)
plot(x,-D2,'-o','linewidth',2,'markersize',6)
ylabel('Deflection Distance (m)')
xlabel('Distance along bar (m)')


y_x1 = A1*(cos(2*pi*x/L)-1);
y_x2 = A2*(cos(2*pi*x/L)-1);


plot(x,y_x1,'linewidth',2)
plot(x,y_x2,'linewidth',2)
legend('No Loading','Measured Light Loading','Measured Heavy Loading','Theoretical Light Loading','Theoretical Heavy Loading')
%%% Moment of Intertia
II = (1/12)*(.02545)*(.00157^3);

E = 2.05*(10^11); %Pa

%%% Theoretical Critical Loading
Pcr = 4*(pi^2)*E*I/(L^2)
dPcr = .002*8*(pi^2)*E*I/(L^3)

%%% Measured Critical Loading
250; %N

%%% Working on shitty integral

fun1 = @(k) sqrt(1+((2*pi*A1/L)*sin(2*pi.*k./L)).^2);
fun2 = @(k) sqrt(1+((2*pi*A2/L)*sin(2*pi.*k./L)).^2);

q1 = integral(fun1,0,L);
q2 = integral(fun2,0,L);

%%% Theoretical End Shortening
short1 = norm((L - q1) * 1000) %mm
short2 = norm((L - q2) * 1000) %mm

%%% Real End Shortening
Loading_delta(end) %mm, stopping at P = 295.8 N


%% 2nd measured one was .1096 mm w/ P = 1509 N




















