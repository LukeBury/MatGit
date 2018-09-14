clear
clc
clf
x=-1:.1:1;
y=1./(1+25*x.^2);
n=11;
hold on

%Piecewise linear interpolation
xa=linspace(-1,1,n+1);
ya=interp1(x,y,xa,'linear');
scatter(xa,ya,'r');
pa=polyfit(xa,ya,n);
Fa=polyval(pa,x);
plot(x,Fa,'r')

%Piecewise cubic spline interpolation
xb=linspace(-1,1,n+1);
yb=interp1(x,y,xb,'spline');
scatter(xb,yb,'g')
pb=polyfit(xb,yb,n);
Fb=polyval(pb,x);
plot(x,Fb,'g')

%Piecewise cubic hermite interpolation
xc=linspace(-1,1,n+1);
yc=interp1(x,y,xc,'pchip');
scatter(xc,yc,n);
pc=polyfit(xc,yc,n);
Fc=polyval(pc,x);
plot(x,Fc,'k')

u=@(v)1./(1+25*v.^2);
v=-1:.01:1;
plot(v,u(v))