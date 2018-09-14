clear
clc

b=pi/180; %for degrees to radians
a=7016;e=.7;i=30*b;m=0;w=30*b;u=398600;
c=[0:5*b:2*pi];
n=length(c);

for ii=1:n
    r=a*(1-e^2)/(1+e*cos(c(ii)));
    p=a*(1-e^2);
    
    r_vec=[r*cos(c(ii));r*sin(c(ii));0];
    
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
    
    hold on
    grid on
    %position
    R_IJK=(R_trans*r_vec);
    plot3(R_IJK(1),R_IJK(2),R_IJK(3),'+')
    plot3(0,0,0,'x')
    axis square
end