clear
clc

x=[1,1]; tol=10^-4;i=0;i_max=40;error=100;
J(1,1)=2*x(1);
J(1,2)=2*x(2);
J(2,1)=2*x(1);
J(2,2)=-1;
%F(x1,x2)=[x^2+y^2-5;x^2-y+1];
while i<i_max && tol<error
    F(1)=x(1)^2+x(2)^2-5;
    F(2)=x(1)^2-x(2)+1;
    a=J;
    b=-F; %Ta and books use this b value
            %Instead of b=J*x-F
    x_old = x;
    x=gauss(a,b);
    x = x_old + x;
    error=norm(F,2);
    i=i+1;
end
x
Iterations=i