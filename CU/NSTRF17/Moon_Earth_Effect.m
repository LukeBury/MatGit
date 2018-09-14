clear
clc
close all
% ------------------------------------------------------------------------
%%% Setting Initial Values
% ------------------------------------------------------------------------
%%% Time Constraints
t0 = 0;
tf = 10000; % sec

%%% Moon Parameters
M_radius = 1737; % km
M_a = 379700; % km
uM = 4904.8695; % km^3 / s^2
nM = sqrt(uM / M_a^3); % rad/s

%%% Earth Parameters
uE = 398600.4418; % km^3 / s^2
E_pos = [M_a, 0, 0]; % km

%%% Initial Hopper State
r0 = [0, M_radius, 0]; % km
% v0 = [0, .00657, 0]; % km/s 
v0 = [0, 1, 0]; % km/s (big jump)
% v0 = [0, 1.25, 0]; % km/s (craziness)


% %%% Wild State
%r0 = [M_radius, 0, 0]; % km
%v0 = [1.06334, .3, 0]; % km/s
% ------------------------------------------------------------------------
%%% Numerically Integrating the State
% ------------------------------------------------------------------------
%%% Setting integrator accuracy
options = odeset('RelTol',1E-10,'AbsTol',1E-10);

%%% Moon-only case
[Tm,Vm] = ode45(@M_EOMIntegrator,[t0:.1:tf],[r0,v0],options,M_radius,uM);

%%% Moon-and-Earth case
[Te,Ve] = ode45(@ME_EOMIntegrator,[t0:.1:tf],[r0,v0],options,M_radius,uM,uE,E_pos,nM);

% ------------------------------------------------------------------------
%%% Plotting the Motions
% ------------------------------------------------------------------------
hold all

ice_linewidth = 6;
%%% Plotting Moon surface
th = 0:.001:2*pi;
x = M_radius * cos(th);
y = M_radius * sin(th);
p1 = plot(x, y,'k','linewidth',ice_linewidth);

%%% Coloring in Moon
n = 5000;
THETA=linspace(0,2*pi,n);
RHO=ones(1,n)*(M_radius);
[X,Y] = pol2cart(THETA,RHO);
h=fill(X,Y,[0.9 0.9 0.9]);

%%% Plotting Moon-only motion
track_linewidth = 3;
p2 = plot(Vm(:,1),Vm(:,2),'color',[0.8500 0.3250 0.0980],'linewidth',track_linewidth);

%%% Plotting Moon-Earth motion
p3 = plot(Ve(:,1),Ve(:,2),'color',[0.9290 0.6940 0.1250],'linewidth',track_linewidth);
PlotBoi('X, km','Y, km')
% xlabel('X (km)'),ylabel('Y (km)'),zlabel('Z (km)')
grid on
axis square
% disp_theta = acos(dot(r0,[Ve(end,1),Ve(end,2),0]) / M_radius^2);
% disp = disp_theta * M_radius;
disp = norm([Vm(end,1),Vm(end,2),0] - [Ve(end,1),Ve(end,2),0]);
% title(['Initial Y Velocity = ',num2str(v0(2)),' km/s | \Delta = ',...
%     num2str(disp),' km'])
% legend([p2 p3],['No Earth Effects'],['Earth Effects'],...
%     'Location','southeast')
% dim = [.65 0 1 .26];
% str = 'To Earth ---->';
% annotation('textbox',dim,'String',str,'FitBoxToText','on');

%%% Case 1
title(['Moon, Vi = [',num2str(v0(1)),', ',num2str(v0(2)),'] km/s, '...
   '\Delta = ',num2str(disp),' km'])
legend([p2 p3],['No Earth Effects'],['Earth Effects'],...
    'Location','northwest')
dim = [.65 0 1 .9];
str = 'To Earth ---->';
annotation('textbox',dim,'String',str,'FitBoxToText','on');
xlim([-30 30])
ylim([1650 2200])
%%% Limits of plot
%xlim([-2e-3 16e-3])
%ylim([1560.795 1560.82])

% ------------------------------------------------------------------------
%%% Moon-only Numerical Integrator
% ------------------------------------------------------------------------
function [ dY ] = M_EOMIntegrator(t,Y,M_radius,uM)
dY = zeros(6,1);

% Unpack the Moon state vector
y = Y(1:3);
dy = Y(4:6);

%%% Ensure we haven't crashed into the Moon
if sqrt(y(1)^2 + y(2)^2) < M_radius
    return
end

% Dynamics of the system
ddy = (-uM/(norm(y)^3))*y;

% Output the derivative of the state
dY(1:3) = dy;
dY(4:6) = ddy;
end

% ------------------------------------------------------------------------
%%% Moon-and-Earth Numerical Integrator (w/ centripetal)
% ------------------------------------------------------------------------
function [ dY ] = ME_EOMIntegrator(t,Y,M_radius,uM,uE,M_pos,nM)
dY = zeros(6,1);

%%% Unpack the Moon-Earth state vector
y = Y(1:3);
dy = Y(4:6);

%%% Ensure we haven't crashed into the Moon
if sqrt(y(1)^2 + y(2)^2) < M_radius
    return
end

%%% Vector to Earth's center
toEarth = M_pos' - y;

%%% Dynamics of the system
ddy = (-uM/(norm(y)^3))*y + (uE/(norm(toEarth)^3))*toEarth - (nM^2)*M_pos';
%%% Centripetal acceleration magnitude
cent = norm((nM^2)*M_pos');
%%% Earth acceleration magnitude
E3 = norm((uE/(norm(toEarth)^3))*toEarth);

%%% Output the derivative of the state
dY(1:3) = dy;
dY(4:6) = ddy;
end
