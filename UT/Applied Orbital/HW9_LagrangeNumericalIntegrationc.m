function HW9_LagrangeNumericalIntegrationc

clear
clc
% clf

GMe=398600.44;
GMm=GMe/81.30056;

% Constant Earth-Moon distance
Rem=[384400,0,0]; %km

%initial conditions of satellite

Rsm=[0.150934288618019*384400,0, 0]; %km - vector from satellite to moon

m=384400-64185.8;
Res=[(1-0.150934288618019)*384400, 0, 0]; %km - vector from earth to satellite
Ves=[0,0,0]; %km/s

w=sqrt((GMe+GMm)/(norm(Rem)^3));

% Time bounds (sec)
t0      = 0;            
tf      = 20*86164;   

options = odeset('RelTol',1E-10,'AbsTol',1E-10);
[T,V]=ode45(@EOMIntegratorc,[t0:10:tf],[Res,Rsm,Ves],options,GMe,GMm,Rem,w);

Re=[V(:,1) V(:,2) V(:,3)];
Rm=[V(:,4) V(:,5) V(:,6)];
Ve=[V(:,7) V(:,8) V(:,9)];

dPosition=[sqrt((V(:,1)-V(1,1)).^2+(V(:,2)-V(1,2)).^2+(V(:,3)-V(1,3)).^2)];
dPosition(1)
dPosition(end)
plot(T,dPosition)
xlabel('Time')
ylabel('Change in Position')




function [ dY ] = EOMIntegratorc(t,Y,GMe,GMm,Rem,w)
dY = zeros(9,1);

% Unpack the state vector
y = Y(1:3);
m = Y(4:6);
dy = Y(7:9);

% Dynamics of the system

ddx=-GMe*y(1)/(norm(y)^3)+GMm*(m(1)/(norm(m)^3)-1/(norm(Rem)^2))+2*w*dy(2)+y(1)*(w^2);
ddy=-GMe*y(2)/(norm(y)^3)+GMm*(m(2)/(norm(m)^3))-2*w*dy(1)+y(2)*(w^2);
ddz=-GMe*y(3)/(norm(y)^3)+GMm*(m(3)/(norm(m)^3));

% Output the DERIVATIVE of the state
dY(1:3) = dy;
dY(4:6) = dy;
dY(7) = ddx;
dY(8) = ddy;
dY(9) = ddz;
