function HW4_HumericalIntegration
%%% Demonstrate numerical integration %%%
clear
clc
clf

u=398600.4415;
Time=11*86164/167;
a=(((Time/(2*pi))^2)*u)^(1/3)
vel=((-u/a+2*u/a)^.5);
%a=6876.921342; %km
%vel=7.6133285761;%km/s
vy=vel*cosd(28.5);
vz=vel*sind(28.5);
% An initial condition:
r= [a,0,0];    %km   
v=[0,vy,vz];      %km/s

% Time bounds (sec)
t0      = 0;            
tf      = 11*86164;   

options = odeset('RelTol',1E-10,'AbsTol',1E-10);
[T,V]=ode45(@EOMIntergrator,[t0,tf],[r,v],options,u);


R=[V(:,1) V(:,2) V(:,3)];
vel=[V(:,4) V(:,5) V(:,6)];
plot3(V(:,1),V(:,2),V(:,3))
xlabel('X (km)'),ylabel('Y (km)'),zlabel('Z (km)')
grid on
axis square

f=[V(end,1) V(end,2) V(end,3)];
m=r-f;
Mis=norm(m) %km


function [ dY ] = EOMIntergrator(t,Y,u)
dY = zeros(6,1);

% Unpack the state vector
y = Y(1:3);
dy = Y(4:6);

% Dynamics of the system
ddy = (-u/(norm(y)^3))*y;

% Output the DERIVATIVE of the state
dY(1:3) = dy;
dY(4:6) = ddy;
