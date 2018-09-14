clear
clc
clf

speed=[30 45 60 75 90 120];
thinking=[5.6 8.5 11.1 14.5 16.7 22.4];
braking=[5.0 12.3 21.0 32.9 47.6 84.7];
hold on
scatter(speed, thinking,'g')
scatter(speed, braking,'b')
xlabel('Speed (km/h)')
ylabel('Meters')
legend('Thinking','Braking')

%linear for thinking
x=speed';
y=thinking';
n=length(x);
sx=sum(x);
sy=sum(y);
sx2=sum(x.*x);
sxy=sum(x.*y);
sy2=sum(y.*y);
a(1)=(n*sxy-sx*sy)/(n*sx2-sx^2);
a(2)=sy/n-a(1)*sx/n;

x_t=linspace(min(x),max(x),2);
y_t=a(1)*x_t+a(2);
plot(x_t,y_t,'g')

%Quadratic for breaking
x=speed';
y=braking';
for i=1:n
    zeta(:,i)=[1,x(i),x(i)^2];
end
zeta=zeta';
A=zeta'*zeta;
B=zeta'*y;
X=(A\B)'; %produces correct values but backwards

u=30:.02:120;
v=@(u).0059*u.^2+.0009*u-.095;
plot(u,v(u))
d=30:.02:120;
s=@(d).1865*d+.08;
distance=s(65)+v(65)

