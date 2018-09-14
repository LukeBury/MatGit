clear
clc
x=[0:.01:3];
L2=-4/9 * x.*(x-3);
L3=1/4 *x.*(2*x-3);

hold on
plot(x,L2)
plot(x,L3)
grid on
