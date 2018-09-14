clear
clc

a=[1 2 0;-2 1 2;1 3 1]
e=100; tol=10^-12;

n=size(a,1);
x=ones(n,1);
x_old=x;

while e>=tol
    y=a*x;
    L=max(y);
    x=y/L;
    e=norm(x-x_old,2)/norm(x_old,2);
    x_old=x;
end
My_eign=L
eig_vector=x'/norm(x)
%normalized eig_vector to match matlab format

%Matlab Answer
[V,d]=eig(a)
Matlab_eign=d(3,3)
Matlab_vector=V(1:3,3)'
