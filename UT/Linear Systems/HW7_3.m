clear all
close all

num     = 1.6^2;
den     = [1 2*.75*1.6 1.6^2];
sys     = tf(num,den);
f       = 1; 

dt      = 0.001;
t       = [0:dt:50];
hold all
%Sinusoidal input
u       = ones(length(t),1)';
y       = lsim(sys,u,t);
plot(t,y)

k=.75; w=1.6;K=sqrt(1-k^2);
f=(1-exp(-k*w*t).*(cos(w*K*t)+k/K*sin(w*K*t)));
plot(t,f)
xlabel('time')
ylabel('y')
