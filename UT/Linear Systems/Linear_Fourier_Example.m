clear all
close all

%Fourier series expansion for the square wave example (lecture notes)

N   = 20;

dt  = 1/10000;        %time step
t1 =[0:dt:pi];
t2 =[pi+dt:dt:2*pi];
t =[t1(:);t2(:)];
T   = 2*pi;           %period

f = [ones(length(t1),1); -ones(length(t2),1)]; %square wave function

% form basis functions
w0  = 2*pi/T;

kk  = 0;
for k = -N : N,
    kk = kk + 1;
    Psi(:,kk) = (1/sqrt(T))*exp(-w0*k*t*1i);
    if abs(round(k/2)-k/2)> 1e-5
        ca(kk,1) = -2*sqrt(2)/(sqrt(pi)*i*k);   %analytically computed coefficient
    else
        ca(kk,1) = 0;
    end
end
c   = Psi'*f*dt;

fn  = Psi*c;    %numerically computed function
fa  = Psi*ca;   %analytically computed function

plot(t,f,'k',t,fa,'r',t,fn,'b--');grid on;  boldify

disp('Comparison of numerically and analytically computed Coefficients') 
disp([c ca])