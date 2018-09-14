clear
clc

a=2*pi*3/86164;b=2*pi*86400/86164; c=2*pi*86403/86164;
rx=[cos(a) -sin(a) 0;sin(a) cos(a) 0;0 0 1];
ry=[cos(b) -sin(b) 0;sin(b) cos(b) 0;0 0 1];
rz=[cos(c) -sin(c) 0;sin(c) cos(c) 0;0 0 1];
k=[0;0;1];


r1=[3325.396441;5472.597483;-2057.129050];
r2=[3309.747175;5485.240159;-2048.664333];
ra=[4389.882255;-4444.406953;-2508.462520];
rb=[4402.505030;-4428.002728;-2515.303456];

r22=rx*r2
dr1=r22-r1
v1=dr1/3
h1=cross(r22,v1)
n1=cross(k,h1)
AN1=acos(n1(1)/norm(n1));
AN1d=AN1*180/pi
i1=acos(h1(3)/norm(h1))
i1d=i1*180/pi

raa=ry*ra
rbb=rz*rb
dr2=rbb-raa
v2=dr2/3
h2=cross(rbb,v2)
n2=cross(k,h2)
AN2=acos(n2(1)/norm(n2));
AN2d=AN2*180/pi
i2=acos(h2(3)/norm(h2))
i2d=i2*180/pi

dAN=AN2-AN1;
dANd=dAN*180/pi
di=i2-i1;
did=di*180/pi

u=398600.4415; a=norm(r1); j=.00108265; ae=6378.1363;n=sqrt(u/(a^3));
nr=(-3*n*j/2)*((ae/a)^2)*cos(i1);
nrf=(nr*180/pi)*60*60*24

% %b
% sigg=(360/365.25)
% sig=(360/365.25)*(pi/180)/(24*60*60)
% aa=6378.1363+705;nn=sqrt(u/(aa^3));
% x=(-sig*2/(3*nn*j))*((aa/ae)^2)
% iii=acos(x)
% iiii=iii*180/pi