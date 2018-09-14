function[I,h,e]=simpson(f,a,b,n)
%n=subintervals
h=(b-a)/n;
I=f(a);
for i=1:n-1
    if mod(i,2)==0
        I=I+2*f(i*h);
    else
        I=I+4*f(i*h);
    end
end

I=I+f(b);
I=I*h/3;
dfourf=@(x) 2400*(-9+20*x);
K=(b-a)/2;
e=abs(-h^4/180 *(b-a)*dfourf(K));
end
