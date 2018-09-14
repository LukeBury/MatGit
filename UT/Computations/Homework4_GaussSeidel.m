clear
clc

a=[10 2 -1;-3 6 2;1 1 5]; b=[27;-61.5;-21.5];
n=size(a,1);
y=0.8;

e=100; tol=10^-12;
x=zeros(n,1);

while e > tol
    xold=x;
    for k = 1:n
        sum1 = 0; sum2 = 0;
        for j = 1:k-1
            sum1 = sum1 + a(k,j)*x(j);
        end
        for j=k+1:n
            sum2=sum2+a(k,j)*x(j);
        end
        x(k)=(b(k)-sum1-sum2)/a(k,k);
        x = y*x + (1-y)*xold;
    end
    e=norm(b-a*x,2)/norm(b,2);
end
x
a\b