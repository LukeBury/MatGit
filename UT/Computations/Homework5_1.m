clear
clc

%Pivoted matrix to be more Gauss-Siedel suitable
a=[-8 1 -2;2 -6 -1;-3 -1 7]; b=[-20;-38;-34];
n=size(a,1);
y=.8;

i=0;
e=100; tol=.0001;
x=zeros(n,1);

while e >= tol
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
        x(k) = y*x(k) + (1-y)*xold(k);
    end
    e=norm(b-a*x,2)/norm(b,2);
    i=i+1;
end
iterations=i
x
a\b
