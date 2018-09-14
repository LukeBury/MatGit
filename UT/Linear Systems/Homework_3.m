clear;clc;clf
n=5;
dt=1/10000;
t1=[0:dt:pi]; 
t2=[dt+pi:dt:2*pi]; 
t=[t1,t2]';T=2*pi;
for ii=1:length(t1)
    f(ii)=sin(t1(ii));
end
for iii=length(t1):length(t)
    f(iii)=0;
end
f=f';
w0=2*pi/T;
kk=0;
for k=-n:n
    kk=kk+1;
    Psi(:,kk)=(1/sqrt(T))*exp(-w0*k*t*i);
    if abs(round(k/2)-k/2)<1e-5
        ca(kk,1)=-(sqrt(2)/(sqrt(pi)*(k^2-1)));
    elseif k==1
        ca(kk,1)=sqrt(pi)*i/(2*sqrt(2));
    elseif k==-1
        ca(kk,1)=-sqrt(pi)*i/(2*sqrt(2));
    else
        ca(kk,1)=0;
    end
end
c=Psi'*f*dt;
fn=Psi*c;
fa=Psi*ca;
plot(t,f,'k',t,fa,'r',t,fn,'b--'); 
grid on;
disp([c ca])
boldify