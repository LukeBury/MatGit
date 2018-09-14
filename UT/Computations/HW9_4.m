clear
clc

%f(x)=.2+25*x-200*x^2+675*x^3-900*x^4+400*x^5;
a=0; b=.8;
n_base=[3 5 7 9 11 13 15]-1;
f=@(x).2+25*x-200*x^2+675*x^3-900*x^4+400*x^5;
for i=1:length(n_base)
    n=n_base(i);
    [I,h,e]=simpson(f,a,b,n);
    I_v(i)=I;
    h_v(i)=h;
    e_v(i)=e;
end
out=[n_base' h_v' I_v' e_v']

fprintf('%20s %20s %20s 25%s \n', 'Sub-intervals,n','Width,h','Error','Approx of Integral')
for ii=1:7
fprintf('%20d %20f %20.10f %25.10f \n',out(ii,:))
end