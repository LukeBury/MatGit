clear
clc

u=398600;t1=1200;t2=2700;rp=7000;ra=9000;

a=(ra+rp)/2;
n=sqrt(u/(a^3));
e=(ra-rp)/(ra+rp);
T=2*pi*sqrt((a^3)/u);

%t1
E=n*t1;
for i=0:100
    M=t1*n;
    dE=(M-E+e*sin(E))/(1-e*cos(E));
    E=E+dE;
end
theta1=2*atan(tan(E/2)*sqrt((1-e)/(1+e)));

%t2
E=n*t2;
for i=0:100
    M=t2*n;
    dE=(M-E+e*sin(E))/(1-e*cos(E));
    E=E+dE;
end
theta2=2*atan(tan(E/2)*sqrt((1+e)/(1-e)));

d_theta=theta2-theta1;
%in degrees
Degrees=d_theta*180/pi