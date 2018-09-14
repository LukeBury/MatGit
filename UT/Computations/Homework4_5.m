clear
clc

a=[10 2 -1;-3 6 2;1 1 5];b=[27;-61.5;-21.5]
%a=[8 3 1;2 4 -1;-6 0 7];b=[12;5;1];
%a=[3 1 -1;1 4 -1;1 1 5];b=[3;4;7];
%a=[-1 3 5; -2 4 5; 0 2 -1];b=[7;-3;1]

n=size(a,1);

e=100; tol=10^-12;
x=zeros(n,1);

while e > tol
    for k = 1:n
        sum1 = 0; sum2 = 0;
        for j = 1:k-1
            sum1 = sum1 + a(k,j)*x(j);
        end
        for j=k+1:n
            sum2=sum2+a(k,j)*x(j);
        end
        x(k)=(b(k)-sum1-sum2)/a(k,k);
    end
    e=norm(b-a*x,2)/norm(b,2);
end
x
a\b