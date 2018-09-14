clear;clc;
u=398600; R1=7000;
Vc1=sqrt(u/R1);
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