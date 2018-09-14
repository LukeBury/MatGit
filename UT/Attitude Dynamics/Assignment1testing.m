clear 
clc
% clf

ib=[3/sqrt(29),4/sqrt(29),2/sqrt(29)];
jb=[-7*sqrt(2/24447),107/sqrt(48894),-193/sqrt(48894)];
kb=[-17*sqrt(2/843),19/sqrt(1686),13/sqrt(1686)];

% 
% hold on
% xlabel('x')
% ylabel('y')
% zlabel('z')
% 
% a=[0,0,0];
% plot(a)
% 
% quiver3(1,2,3,ib(1),ib(2),ib(3),'r')
% quiver3(1,2,3,jb(1),jb(2),jb(3),'b')
% quiver3(1,2,3,kb(1),kb(2),kb(3),'g')
% 
% 
% 
% quiver3(1,2,3,3,4,2,'m')
% quiver3(1,2,3,2,7,-5,'m')
% quiver3(4,6,5,-1,3,-7,'m')
% 

ia=[1,0,0];ja=[0,1,0];ka=[0,0,1];

% Cba=zeros(3,3);
% 
% Cba(1,1)=ib(1);
% Cba(1,2)=ib(2);
% Cba(1,3)=ib(3);
% Cba(2,1)=jb(1);
% Cba(2,2)=jb(2);
% Cba(2,3)=jb(3);
% Cba(3,1)=kb(1);
% Cba(3,2)=kb(2);
% Cba(3,3)=kb(3);
% 
% Cab=Cba';
% Cab*[2;-1;3];

% syms c1 c2 c3 s1 s2 s3;
% 
% R1=[c3 s3 0;-s3 c3 0;0 0 1];
% R2=[c2 0 -s2;0 1 0;s2 0 c2];
% R3=[1 0 0;0 c1 s1;0 -s1 c1];
syms c1 c2 s1 s2;c3=1;s3=0;

R1=[c3 s3 0;-s3 c3 0;0 0 1];
R2=[c2 0 -s2;0 1 0;s2 0 c2];
R3=[1 0 0;0 c1 s1;0 -s1 c1];
R1*(R2*R3)













