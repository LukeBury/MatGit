clear
clc

num     = 1.6^2;
den     = [1 2*.75*1.6 1.6^2];
sys     = tf(num,den);
f       = 1; 

dt=.001;
t=[0:dt:40];

g=impulse(sys,t);

u=sin(t*f);
y=lsim(sys,u,t);

G=fft(g)*dt;
bode(sys)

func=polyval(num,i*f)/polyval(den,i*f);
A=abs(func)
phi=angle(func)