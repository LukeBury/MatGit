function HW5_NumericalIntegration
clear
clc

u=398600.4415;
radius=.0037/2;
A=pi*(radius^2); %km^2
m=400; %kg
Cd=2;
gamma=.006; %kg/km^3
Time=11*86164/167; %sec
ai=(((Time/(2*pi))^2)*u)^(1/3); %km
vel=((-u/ai+2*u/ai)^.5); %km/s
vy=vel*cosd(28.5); %km/s
vz=vel*sind(28.5); %km/s
% An initial condition:
r= [ai,0,0];    %km   
v=[0,vy,vz];      %km/s

% Time bounds (sec)
t0      = 0;            
tf      = 11*86164;   

options = odeset('RelTol',1E-10,'AbsTol',1E-10,'Events',@events);
[T,V,te,Ye,ie]=ode45(@EOMIntergrator,[t0,tf],[r,v],options,u,Cd,A,m,gamma);


R=[V(:,1) V(:,2) V(:,3)]; %km
vel=[V(:,4) V(:,5) V(:,6)]; %km/s
misclosure=norm(r-Ye(end,1:3)) %km

%find groundtrack misclosure
we=2*pi/86164;
mislong=(we*(T(end)-te(end)))*180/pi %misclosure in longitude (degrees)
for i=1:1:length(T)
    a(i,1)=norm(R(i,:));
end
poly=polyfit(T,a,1);
kmperday=poly(1,1)*60*60*24 %taking slope from km/sec to km/day



function [ dY ] = EOMIntergrator(t,Y,u,Cd,A,m,gamma)
dY = zeros(6,1);

% Unpack the state vector
y = Y(1:3);
dy = Y(4:6);

% Dynamics of the system
ddy = (-u/(norm(y)^3))*y-.5*(Cd*A/m)*gamma*norm(dy)*dy;

% Output the DERIVATIVE of the state
dY(1:3) = dy;
dY(4:6) = ddy;

function [value,isterminal,direction] = events(t,Y,u,Cd,A,m,gamma)
%find the points where the z comp = zero
value = Y(3);    
isterminal = 0;   
direction = 0;   
