clear
clc

t=10^-8;
a=[t^2,1;1,1];
b=[1;0];
cond(a)

n=size(a,1);

l=0;
k=1;

%Forward Elimination
for k=1:n-1
    for i=k+1:n
        l(i,k)=a(i,k)/a(k,k);
        a(i,k)=0;
        for j=k+1:n
            a(i,j)=a(i,j)-l(i,k)*a(k,j);
        end
        b(i)=b(i)-l(i,k)*b(k);
    end
end


%Backwards Substitution
x(n)=b(n)/a(n,n); %defines the function
for i=n-1:-1:1
    sum=0;
    for j=i+1:n
        sum=sum+a(i,j)*x(j);
    end
    x(i)=(b(i)-sum)/a(i,i);
end
x

x=a\b