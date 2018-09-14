clear
clc
clf
%n= number of point;
n=100; L=3; EA=100; f=1; dx=L/n;
x=[0:dx:L];

K=diag(ones(n-1,1),1)+diag(ones(n-1,1),-1);
for i=1:n
    for j=1:n
        if i==j
            K(i,j)=-2;
        end
    end
end
K(n,n)=-1;
F=-f/EA*dx^2*ones(n,1);
F(n)=0;
u=K\F;

for o=1:n
    U(o+1)=u(o);
end
plot(x,U)