clear
clc
close all
% ------------------------------------------------------------------------
%%% Setting Initial Values
% ------------------------------------------------------------------------
%%% Time Constraints
t0 = 0;
tf = 8000; % sec

%%% Enceladus Parameters
E_radius = 252; % km
E_a = 237948; % km
uE = 7.2095; % km^3 / s^2
nE = sqrt(uE / E_a^3); % rad/s

%%% Saturn Parameters
uS = 37931187; % km^3 / s^2
S_pos = [E_a, 0, 0]; % km

%%% Initial Hopper State
r0 = [0, E_radius, 0]; % km
% v0 = [0, .00657, 0]; % km/s 
v0 = [0, .02, 0]; % km/s 
%v0 = [0, .0111, 0]; % km/s


% ------------------------------------------------------------------------
%%% Numerically Integrating the State
% ------------------------------------------------------------------------
%%% Setting integrator accuracy
options = odeset('RelTol',1E-10,'AbsTol',1E-10);

%%% Enceladus-only case
[Te,Ve] = ode45(@E_EOMIntegrator,[t0:.1:tf],[r0,v0],options,E_radius,uE);

%%% Enceladus-and-Saturn case
[Ts,Vs] = ode45(@ES_EOMIntegrator,[t0:.1:tf],[r0,v0],options,E_radius,uE,uS,S_pos,nE);

% ------------------------------------------------------------------------
%%% Plotting the Motions
% ------------------------------------------------------------------------
hold all
ice_linewidth = 6;
%%% Plotting Enceladus surface
th = 0:.001:2*pi;
x = E_radius * cos(th);
y = E_radius * sin(th);
p1 = plot(x, y,'linewidth',ice_linewidth);

%%% Coloring in Enceladus
n = 5000;
THETA=linspace(0,2*pi,n);
RHO=ones(1,n)*(E_radius);
[X,Y] = pol2cart(THETA,RHO);
h=fill(X,Y,[.5 .85 1]);

%%% Plotting Enceladus-only motion
track_linewidth = 3;
p2 = plot(Ve(:,1),Ve(:,2),'color',[0.8500 0.3250 0.0980],'linewidth',track_linewidth);

%%% Plotting Enceladus-Saturn motion
p3 = plot(Vs(:,1),Vs(:,2),'color',[0.9290 0.6940 0.1250],'linewidth',track_linewidth);
PlotBoi('X, km', 'Y, km')
% xlabel('X (km)'),ylabel('Y (km)'),zlabel('Z (km)')
grid on
axis square
%disp_theta = acos(dot(r0,[Vs(end,1),Vs(end,2),0]) / E_radius^2);
%disp = disp_theta * E_radius;
disp = norm([Vs(end,1),Vs(end,2),0] - [Ve(end,1),Ve(end,2),0]);
% title(['Initial Y Velocity = ',num2str(v0(2)),' km/s | \Delta = ',...
%     num2str(disp),' km'])
% legend([p2 p3],['No Saturn Effects'],['Saturn Effects'],...
%     'Location','southeast')
% dim = [.65 0 1 .26];
% str = 'To Saturn ---->';
% annotation('textbox',dim,'String',str,'FitBoxToText','on');

title(['Enceladus, Vi = [',num2str(v0(1)),', ',num2str(v0(2)),'] km/s'...
    ', \Delta = escape'])
legend([p2 p3],['No Saturn Effects'],['Saturn Effects'],...
    'Location','northwest')
dim = [.65 0 1 .9];
str = 'To Saturn ---->';
annotation('textbox',dim,'String',str,'FitBoxToText','on');
% xlim([-23 23])
% ylim([244 260])
xlim([-50  50])
ylim([220 280])
% ------------------------------------------------------------------------
%%% Enceladus-only Numerical Integrator
% ------------------------------------------------------------------------
function [ dY ] = E_EOMIntegrator(t,Y,E_radius,uE)
dY = zeros(6,1);

% Unpack the Enceladus state vector
y = Y(1:3);
dy = Y(4:6);

%%% Ensure we haven't crashed into Enceladus
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
%%% Enceladus-and-Saturn Numerical Integrator (w/ centripetal)
% ------------------------------------------------------------------------
function [ dY ] = ES_EOMIntegrator(t,Y,E_radius,uE,uS,S_pos,nE)
dY = zeros(6,1);

%%% Unpack the Enceladus-Saturn state vector
y = Y(1:3);
dy = Y(4:6);

%%% Ensure we haven't crashed into the Moon
if sqrt(y(1)^2 + y(2)^2) < E_radius
    return
end

%%% Vector to Saturn's center
toSaturn = S_pos' - y;

%%% Dynamics of the system
ddy = (-uE/(norm(y)^3))*y + (uS/(norm(toSaturn)^3))*toSaturn - (nE^2)*S_pos';
fprintf('------------------')

%%% Centripetal acceleration magnitude
cent = norm((nE^2)*S_pos');
%%% Saturn acceleration magnitude
S3 = norm((uS/(norm(toSaturn)^3))*toSaturn);

%%% Output the derivative of the state
dY(1:3) = dy;
dY(4:6) = ddy;
end
