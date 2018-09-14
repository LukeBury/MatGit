clear
clc

c=1; x_old=1; i=0;
%function = x^2-x;
a=0; b=3;
tol=10^-12;
iter=50;
error=100;
while (error>tol) && (i<iter)
    f_a=a^2-2;
    f_b=b^2-2;
    c=b-(f_b*(a-b))/(f_a-f_b);
    f_c=c^2-2;
    if (f_a*f_c)<0
        b=c;
    elseif (f_b*f_c)<0
        a=c;
    end
    error=abs(x_old-c);
    x_old=c;
    i=i+1;
end
c
error
i