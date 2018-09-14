clear
clc

r1=[5039.04 4228.257 0];
r2=[-34538.727 24184.277 0];
c=r2-r1;
R1=6578;
R2=42164;
C=44324.286;
a=24371;

% alpha=155.4178*pi/180;
% beta=24.5822*pi/180;
u=398600.4415;
% n=1.6594297E-4;
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