clear
clc
x=[-2 0 2];
y=[4 2 8];
syms X
L(1)=(X-x(2))/(x(1)-x(2))*(X-x(3))/(x(1)-x(3));
L(2)=(X-x(1))/(x(2)-x(1))*(X-x(3))/(x(2)-x(3));
L(3)=(X-x(1))/(x(3)-x(1))*(X-x(2))/(x(3)-x(2));
%showing the polynomials
L1=simplify(L(1))
L2=simplify(L(2))
L3=simplify(L(3))
p=2;sum=0;
for i=1:p+1
    sum=sum+y(i)*L(i);
    f=sum;
end
Lagrange_Polynomial=simplify(f)

hold on
scatter(x,y)
v=@(u)u.^2+u+2;
u=-5:.1:5;
plot(u,v(u))