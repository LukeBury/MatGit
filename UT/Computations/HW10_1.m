clear
clc
clf
x=[0:.01:3];
u=@(x) (x.*(2*3-x))/200

plot(x,u(x))