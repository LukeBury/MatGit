clear
clc

GMe=398600.44;
GMm=GMe/81.30056;

% Constant Earth-Moon distance
Rem=384400; %km

psi=GMm/(GMm+GMe);

A=psi/(1-psi);
syms x
f=  @(x) (x^3)*(x^2-3*x+3)/((1+x+x^2)*((1-x)^3))-A;
df= @(x) (3*x^3*(x^2 - 3*x + 3))/((x - 1)^4*(x^2 + x + 1)) - (3*x^2*(x^2 - 3*x + 3))/((x - 1)^3*(x^2 + x + 1)) - (x^3*(2*x - 3))/((x - 1)^3*(x^2 + x + 1)) + (x^3*(2*x + 1)*(x^2 - 3*x + 3))/((x - 1)^3*(x^2 + x + 1)^2);


x0=.1; %guess

e=1e-15;

for i=1:25
    y=f(x0);
    yprime=df(x0);
    x1=x0-y/yprime;
    if (abs(x1-x0)/abs(x1)<e)
        break
    end
    x0=x1;
end
x0