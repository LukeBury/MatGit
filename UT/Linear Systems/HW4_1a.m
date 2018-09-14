clear
clc
clf
N=101;
dt=2*pi/(N-1);
t1=[0:dt:pi]';
t2=[pi+dt:dt:2*pi]';
t=[t1(:);t2(:)];
ff=[t1;zeros(length(t2),1)];
w=exp(i*2*pi/N);


for k=0:N-1;
    v(k+1,1)=w^(-k);
end

for k=1:100;
    T(:,1)=v.^0;
    T(:,k+1)=v.^(k);
end

a=fft(ff);
aa=T*ff;

f   = 0:1:N-1;
f   = f/N;
hold all
axis square
plot(f,abs(a))
plot(f,abs(aa),'r.')
title('Frequency')

fi  = ifft(a);
ffi = conj(T)*aa/N;

figure
hold all
axis square
plot(f,abs(fi))
plot(f,abs(ffi),'r.')
title('IDFT')