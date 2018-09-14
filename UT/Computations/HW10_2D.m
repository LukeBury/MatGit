clear
clc
clf
x=0:.01:3;
L2=-4/9*x.*(x-3);
L3=(1/9)*x.*(2*x-3);

hold on
plot(x,L2)
plot(x,L3)
grid on
hold off
figure

s=[0:.01:3];
u=@(s) (s.*(2.*s-3))/200-(3*s.*(s-3))/200;
plot(s,u(s))