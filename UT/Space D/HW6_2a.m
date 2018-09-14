clear;clc;
u=398600; R1=7000;
Vc1=sqrt(u/R1);
for i=1:94
    R2=7000+1000*(i-1);
    a=(R1+R2)/2;
    Vtp=sqrt(2*u/R1-u/a);
    Vc2=sqrt(u/R2);
    Vta=sqrt(2*u/R2-u/a);
    DVh(i)=(Vtp-Vc1)+(Vc2-Vta);
    R(i)=R2/R1;
end
axis square
hold all
boldify
plot(R,DVh)

for i=1:94
    R2=7000+1000*(i-1);
    Ra=4*R2;
    a1=(R1+Ra)/2;
    a2=(R2+Ra)/2;
    Vc2=sqrt(u/R2);
    Vap=sqrt(2*u/R1-u/a1);
    Vba=sqrt(2*u/Ra-u/a2);
    Vaa=sqrt(2*u/Ra-u/a1);
    Vbp=sqrt(2*u/R2-u/a2);
    DVb(i)=(Vap-Vc1)+(Vba-Vaa)+(Vbp-Vc2);
    R(i)=R2/R1;
end
plot(R,DVb)

for i=1:94
    if DVb(i) < DVh(i)
        display('Bi-elliptic transfers are more efficient at values of R greater than:')
        R(i-1)
        break
    end
end
