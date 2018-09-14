clear
clc

x=-1:.02:1;
y=1./(1+25*x.^2);
plot(x,y)

n=10;
xp=linspace(-1,1,n+1);
yp=1./(1+25*xp.^2);
hold on
scatter(xp,yp,'r')
p=polyfit(xp,yp,n);
F=polyval(p,x);
plot(x,F,'r')

for i=1:n+1
    xb(i)=cos((2*i-1)/11*pi/2);
end
nb=size(xb,2);
yb=1./(1+25*xp.^2);
scatter(xb,yb,'g')
pb=polyfit(xb,yb,nb);
Fb=polyval(pb,xb,'g');
plot(xb,Fb,'g')

legend('F(x)','','a','','b')