clear
clc

N=201;
dt=2*pi/(N-1);
t=[0:dt:2*pi]';
ff=[.4*cos(40*pi*t)+.3*sin(10*pi*t)+2*sin(pi*t/5)];
w=exp(i*2*pi/N);

for k=0:N-1;
    v(k+1,1)=w^(-k);
end

for k=1:200;
    T(:,1)=v.^0;
    T(:,k+1)=v.^(k);
end

a=fft(ff);
aa=T*ff;

f=0:1:N-1;
f=f/N;
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

%filtering
ff_filt=aa;
k=dt*N;
for i=1:N/2
    if i>=k
        ff_filt(i)=0;
    end
end
for i=round ((N/2)):N
    if i<=N-k
        ff_filt(i)=0;
    end
end

figure
plot(f,abs(ff_filt))
title('Frequency Filtered')

ffi_filt=conj(T)/N*ff_filt;
figure
plot(t,abs(ffi_filt))
title('IDFT Filtered')
