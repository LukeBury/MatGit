clc; clear

syms w M1 M2 k

a = -3*w^2;
b = 0;
c = 0;
d = (7/8)*(M2/M1)*w^2;

A = [0 0 1 0; 
    0 0 0 1; 
    a b 0 2*w; 
    c d -2*w 0]; 

simplify(eig(A))