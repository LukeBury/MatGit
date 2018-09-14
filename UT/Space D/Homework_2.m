clear
clc
clf
 
a = 7000; e = 0; 
theta = [0:pi/32:2*pi];
p = a*(1-e^2);
r = p./(1+e*cos(theta));
x=r.*cos(theta);
y=r.*sin(theta);
 
hold on 
plot(x,y,'r')
 
a = 7000; e = .5; 
p = a*(1-e^2);
r = p./(1+e*cos(theta));
x=r.*cos(theta);
y=r.*sin(theta);
 
hold on 
plot(x,y,'b')

a = 12000; e = .1; 
p = a*(1-e^2);
r = p./(1+e*cos(theta));
x=r.*cos(theta);
y=r.*sin(theta);
 
hold on 
plot(x,y,'g')
 
  
a = 17000; e = .9; 
p = a*(1-e^2);
r = p./(1+e*cos(theta));
x=r.*cos(theta);
y=r.*sin(theta);
 
hold on 
plot(x,y,'k')
 
legend('a=7000,e=0','a=7000, e=.5','a=12000, e=.1','a=17000, e=.9')
