clear
clc
clf

l=3; EA=100; f=1; K=zeros(2,2); F=zeros(2,1);n=2;u=0;
J=l/2;
ksi=[-1/sqrt(3),1/sqrt(3)];
x =l/2+l/2*ksi;
w =[1,1];

for k=1:n  %goes through gauss points 
    for i=1:n
        for j=1:n
            K(i,j)=K(i,j)+dL(i,x(k))*dL(j,x(k))*w(k)*J*EA;
        end
        F(i)=F(i)+f*L(i,x(k))*J*w(k);
    end
end
K
F
beta=K\F

for e=1:2
    syms X
    u=u+beta(e)*L(e,X);
end
u

s=(0:.01:3);

u=@(s) (s.*(2.*s-3))/200-(3*s.*(s-3))/200; %sum i=1 to 2 L_i(x)*beta_i
plot(s,u(s))