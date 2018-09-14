clear
clc
close all
% ------------------------------------------------------------------------
%%% Setting Initial Values
% ------------------------------------------------------------------------
%%% Time Constraints
t0 = 0;
% tf = 2374; % sec
 tf = 37194; % sec

%%% Europa Parameters
E_radius = 1560.8; % km
E_a = 671100; % km
uE = 3203.413216; % km^3 / s^2
nE = sqrt(uE / E_a^3); % rad/s

%%% Jupiter Parameters
uJ = 126672520; % km^3 / s^2
J_pos = [E_a, 0, 0]; % km

%%% Initial Hopper State
% r0 = [0, E_radius, 0]; % km
% v0 = [0, .00657, 0]; % km/s (for 10 second unperturbed flight)
% v0 = [0, 1, 0]; % km/s (big jump)
% v0 = [0, 1.25, 0]; % km/s (craziness)
% v0 = [-1.5, .3, 0]; % km/s (energy loss)

%%% Wild State
r0 = [E_radius, 0, 0]; % km
v0 = [1.06334, .3, 0]; % km/s
% ------------------------------------------------------------------------
%%% Numerically Integrating the State
% ------------------------------------------------------------------------
%%% Setting integrator accuracy
options = odeset('RelTol',1E-10,'AbsTol',1E-10);

%%% Europa-only case
[Te,Ve] = ode45(@E_EOMIntegrator,[t0:.1:tf],[r0,v0],options,E_radius,uE);

%%% Europa-and-Jupiter case
[Tj,Vj] = ode45(@EJ_EOMIntegrator,[t0:.1:tf],[r0,v0],options,E_radius,uE,uJ,J_pos,nE);

% ------------------------------------------------------------------------
%%% Plotting the Motions
% ------------------------------------------------------------------------
hold all
ice_linewidth = 6;
%%% Plotting Europa surface
th = 0:.001:2*pi;
x = E_radius * cos(th);
y = E_radius * sin(th);
p1 = plot(x, y,'linewidth',ice_linewidth);

%%% Coloring in Europa
n = 5000;
THETA=linspace(0,2*pi,n);
RHO=ones(1,n)*(E_radius);
[X,Y] = pol2cart(THETA,RHO);
h=fill(X,Y,'c');

%%% Plotting Europa-only motion
track_linewidth = 3;
p2 = plot(Ve(:,1),Ve(:,2),'linewidth',track_linewidth);

%%% Plotting Europa-Jupiter motion
p3 = plot(Vj(:,1),Vj(:,2),'linewidth',track_linewidth);
PlotBoi('X, km','Y, km')
% xlabel('X (km)'),ylabel('Y (km)'),zlabel('Z (km)')
% set(gcf,'color','white')
grid on
axis square
% disp_theta = acos(dot(r0,[Vj(end,1),Vj(end,2),0]) / E_radius^2) * 180/pi
% disp = disp_theta * E_radius
disp = norm([Vj(end,1),Vj(end,2),0] - [Ve(end,1),Ve(end,2),0]);
Vf = sqrt(Vj(end,4)^2+Vj(end,5)^2);
% title(['Vi = [',num2str(v0(1)),', ',num2str(v0(2)),'] km/s | \Delta = ',...
%     num2str(disp),' | Vf = ',num2str(Vf)...
%     ,' km/s'])

% title(['Vi = ',num2str(sqrt(v0(1)^2+v0(2)^2)),' km/s | \Delta = ',...
%     num2str(disp),' | Vf = ',num2str(Vf)...
%     ,' km/s'])
% legend([p2 p3],['No Jupiter Effects'],['Jupiter Effects'],...
%     'Location','southeast')
% dim = [.65 0 1 .26];
% str = 'To Jupiter ---->';
% annotation('textbox',dim,'String',str,'FitBoxToText','on');

% %%% Case 1
% title(['Europa, Vi = [',num2str(v0(1)),', ',num2str(v0(2)),'] km/s'...
%     ', \Delta = ',num2str(disp),' km'])
% legend([p2 p3],'No Jupiter Effects','Jupiter Effects',...
%     'Location','northwest')
% dim = [.65 0 1 .9];
% str = 'To Jupiter ---->';
% annotation('textbox',dim,'String',str,'FitBoxToText','on');
% xlim([-200 800])
% ylim([1250 2300])


%%% Case Wild
title(['Europa, Vi = [',num2str(v0(1)),', ',num2str(v0(2)),'] km/s'])
legend([p2 p3],'No Jupiter Effects','Jupiter Effects',...
    'Location','northwest')
dim = [.65 0 1 .9];
str = 'To Jupiter ---->';
annotation('textbox',dim,'String',str,'FitBoxToText','on');
xlim([-.3*E_radius 2.3*E_radius])
ylim([-1.1*E_radius 1.5*E_radius])


%%% Limits of plot
% xlim([-2e-3 16e-3])
% ylim([1560.795 1560.82])

% ------------------------------------------------------------------------
%%% Stating Assumptions
% ------------------------------------------------------------------------
fprintf('ASSUMPTIONS:\n')
fprintf('-No atmosphere\n')
fprintf('-No rotation (other than basic centripetal force)\n')
fprintf('-Planet is point mass\n')
fprintf('-Moon is point mass at center\n')
fprintf('-Moon is circular\n')
fprintf('-All 2-D\n')

% ------------------------------------------------------------------------
%%% Europa-only Numerical Integrator
% ------------------------------------------------------------------------
function [ dY ] = E_EOMIntegrator(t,Y,E_radius,uE)
dY = zeros(6,1);

% Unpack the Europa state vector
y = Y(1:3);
dy = Y(4:6);

%%% Ensure we haven't crashed into Europa
if sqrt(y(1)^2 + y(2)^2) < E_radius
    return
end

% Dynamics of the system
ddy = (-uE/(norm(y)^3))*y;

% Output the derivative of the state
dY(1:3) = dy;
dY(4:6) = ddy;
end

% ------------------------------------------------------------------------
%%% Europa-and-Jupiter Numerical Integrator (w/ centripetal)
% ------------------------------------------------------------------------
function [ dY ] = EJ_EOMIntegrator(t,Y,E_radius,uE,uJ,J_pos,nE)
dY = zeros(6,1);

%%% Unpack the Europa-Jupiter state vector
y = Y(1:3);
dy = Y(4:6);

%%% Ensure we haven't crashed into Europa
if sqrt(y(1)^2 + y(2)^2) < E_radius
    return
end

%%% Vector to Jupiter's center
toJu = J_pos' - y;

%%% Dynamics of the system
ddy = (-uE/(norm(y)^3))*y + (uJ/(norm(toJu)^3))*toJu - (nE^2)*J_pos';
%%% Centripetal acceleration magnitude
cent = norm((nE^2)*J_pos');
%%% Jupiter acceleration magnitude
J3 = norm((uJ/(norm(toJu)^3))*toJu);

%%% Output the derivative of the state
dY(1:3) = dy;
dY(4:6) = ddy;
end
