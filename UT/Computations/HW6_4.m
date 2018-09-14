clear
clc

%f(x)=2sin(sqrt(x))-x
x=.5;x_old=0;
tol=.0001;
error=inf;
i=0;
while(error>tol)
    x=2*sin(sqrt(x));
    error=abs(x-x_old);
    x_old=x;
    i=i+1;
end
x
i