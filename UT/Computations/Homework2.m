clear
clc

A=[4 7;1 2;5 6]';
B=[1 0 1;0 2 3;4 0 0];

[N,M]=size(A);
C=0;

if size(A,2)~=size(B,1)
    fprintf('Multiplication does not work \n')
else
    fprintf('Multiplication works \n')
    for j=1:N
        for k=1:M
            sum=0;
            for l=1:M
               sum = sum + A(j,l)*B(l,k);
            end
            C(j,k)=sum;
        end
    end
end
A*B;
C
