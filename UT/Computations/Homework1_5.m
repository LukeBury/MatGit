clear
clc

x0 = 30*pi/180;
x1 = 31*pi/180;
y=0;

syms x;
f=sin(x);
for n = 1:10
    t=taylor(f,x,n,x0);
    y(n)=subs(t,x,x1)
end

e = abs(sin(x1)-y)/sin(x1)*100;
a = [y' e'];

for i = 1:10
    fprintf('Taylor value of sign(31), %10.7f. \tThe percent error, %5.15f%% \n',a(i,:))
end
