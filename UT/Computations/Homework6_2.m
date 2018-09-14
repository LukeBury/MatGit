clear
clc
f_x=1; i=0;
%function = x^2-x;
a=0; b=3;
tol=10^-12;
iter=50;

while (abs(f_x)>tol) && (i<iter)
    x=(a+b)/2;
    f_x=x^2-2;
    f_a=a^2-2;
    f_b=b^2-2;
    if (f_a*f_x)<0
        b=x;
        a=a;
    elseif (f_b*f_x)<0
        a=x;
        b=b;
    end
    i=i+1;
end
i
root=x
error = abs(b-a)/2^i