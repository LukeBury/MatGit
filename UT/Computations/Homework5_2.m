clear
clc
clf

t=0:.01:1;



hold on
plot3(t,t,0*t)
plot3(0*t,0*t,t)
plot3(t,-t,0*t)

grid on
axis square