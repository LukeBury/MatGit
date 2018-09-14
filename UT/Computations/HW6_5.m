clear
clc
g=9.81;cD=.33; v=40; t=5;
error=inf;tol=10^-12;

m=.5; m_old=1;
%y=sqrt(g*m/cD)*tanh(sqrt(g*cD/m)*t)-v;
%Dy=.5*sqrt(g/(m*cD))*tanh(sqrt(g*cD/m)*t)-g/(2*m)*t*sech(sqrt(g*cD/m)*t)^2;

while error>tol
    y_m=sqrt(g*m/cD)*tanh(sqrt(g*cD/m)*t)-v;
    Dy_m=.5*sqrt(g/(m*cD))*tanh(sqrt(g*cD/m)*t)-g/(2*m)*t*sech(sqrt(g*cD/m)*t)^2;
    m=m_old-y_m/Dy_m;
    error=abs(m-m_old)/abs(m);
    m_old=m;
end
m

fun=@(x) (g*x/cD)^.5 * tanh(sqrt(g*cD/x)*t)-v;
fzero(fun,m)