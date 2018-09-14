clear
clc

r2=[14420.98524 39621.3313 0];
r1=[0 -6978 0];
t=12600;
u=398600.4415;
R1=norm(r1);
R2=norm(r2);
c=r2-r1;
C=norm(c);

a=(R1+R2)/2-1;
da=1;

while da > 1E-4
    a=a+da;
    alpha=2*asind(((R1+R2+C)/(4*a))^.5)*pi/180;
    beta=2*asind(((R1+R2-C)/(4*a))^.5)*pi/180;

    n=sqrt(u/(a^3));
    Dt=(1/n)*(alpha-beta-(sin(alpha)-sin(beta)));
    dt=12600-Dt;

    dDt=1.5*sqrt(a/u)*(alpha-beta-(sin(alpha)-sin(beta)))+(1/(2*sqrt(u)))*(-((1-cos(alpha))/(cos(alpha/2)))*sqrt(R1+R2+C)+(1-cos(beta))*sqrt(R1+R2-C)/cos(beta/2));

    da=dt/dDt;
end
a

psi=(alpha-beta)/2;
sin(psi);
x=(R1-R2)^2-C^2;
e=sqrt(1+x/(4*a^2*(sin(psi)^2)))