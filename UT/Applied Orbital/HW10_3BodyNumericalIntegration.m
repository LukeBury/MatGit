function HW10_3BodyNumericalIntegration

clear
clc
clf

GMe=398600.4415; %km^3/s^2
GMm=GMe/81.30056;
GMs=1.327124e11;

%Moon
Rm=[-291608,-274979,36271]; %km
Vm=[.643535,-.730984,-.011506]; %km/s

%Sun
Rs=[26499000,-144697300,0]; %km
Vs=[29.794,5.469,0] %km/s

%Moon to Sun
Rms=Rs-Rm;
Vms=Vs-Vm;


% Time bounds (sec)
t0      = 0;            
tf      = 400*86164;   

options = odeset('RelTol',1E-12,'AbsTol',1E-12);
[T,V]=ode45(@EOMIntegrator,[t0:100:tf],[Rm,Rs,Vm,Vs],options,GMe,GMm,GMs);

Rm=[V(:,1) V(:,2) V(:,3)];
Rs=[V(:,4) V(:,5) V(:,6)];
Vm=[V(:,7) V(:,8) V(:,9)];
Vs=[V(:,10) V(:,11) V(:,12)];

m=[0 0 1];
for k=1:length(Rm)
    h(k,:)=cross(Rm(k,:),Vm(k,:));
    i(k,:)=acosd(h(k,3)/(sqrt(h(k,1)^2+h(k,2)^2+h(k,3)^2)));
    n(k,:)=cross(m,h(k,:));
    omega(k,:)=acosd(n(k,1)/(sqrt(n(k,1)^2+n(k,2)^2+n(k,3)^2)));
end


plot(T,i)
xlabel('Time')
ylabel('Inclination')
figure
plot(T,omega)
xlabel('Time')
ylabel('RAAN')
fit=polyfit(T,omega,1)
rate=fit(1)



function [ dY ] = EOMIntegrator(t,Y,GMe,GMm,GMs)
dY = zeros(12,1);

% Unpack the state vector
m = Y(1:3);
s = Y(4:6);
vm = Y(7:9);
vs = Y(10:12);

% Dynamics of the system

ddx=-(GMe+GMm)*m/(norm(m)^3)+GMs*((s-m)/(norm(s-m)^3)-s/(norm(s)^3));
ddy=-(GMe+GMs)*s/(norm(s)^3);

% Output the DERIVATIVE of the state
dY(1:3) = vm;     %vel moon
dY(4:6) = vs;      %vel sun
dY(7:9) = ddx;     %acc moon
dY(10:12) = ddy;    %acc sun

