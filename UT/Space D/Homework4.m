clear
clc

b=pi/180; %for degrees to radians
a=7016;e=.05;i=30*b;m=0;w=30*b;c=10*b;u=398600;

r=a*(1-e^2)/(1+e*cos(c));
p=a*(1-e^2);

r_vec=[r*cos(c);r*sin(c);0];
v_vec=sqrt(u/p)*[-sin(c);e+cos(c);0];

R11=cos(m)*cos(w)-sin(m)*sin(w)*cos(i);
R12=-cos(m)*sin(w)-sin(m)*cos(w)*cos(i);
R13=sin(m)*sin(i);
R21=sin(m)*cos(w)+cos(m)*sin(w)*cos(i);
R22=-sin(m)*sin(w)+cos(m)*cos(w)*cos(i);
R23=-cos(m)*sin(i);
R31=sin(w)*sin(i);
R32=cos(w)*sin(i);
R33=cos(i);
R_trans=[R11,R12,R13;R21,R22,R23;R31,R32,R33];

%position
R_IJK=(R_trans*r_vec)
%velocity
V_IJK=(R_trans*v_vec)

