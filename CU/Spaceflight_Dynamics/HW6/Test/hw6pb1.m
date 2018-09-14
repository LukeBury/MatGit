clear


mu=398600.441; %%km^3/s^2
i=deg2rad(51.6441);
O=deg2rad(225.6859);
e=.0006691;
w=deg2rad(42.3347);
M0=deg2rad(63.3332);
n=15.54051996*2*pi/24/60/60;%%rev per day to rad per sec
a=(mu/n^2)^(1/3);%km
JD=JulianD(10,11,2016,1,0);
JDi=JulianD(10,5,2016,22,6.1144384);
tf=(JD-JDi)*60*24*60;
M=M0+n*tf;


Tu=(JD-2451545)/36525;
we=7.2921158553*10^(-5);
%thetaG0=100.460618+36000.7705361*Tu+.00038793*Tu^2-2.6*10^(-8)*Tu^3;
thetaG0=(67310.54841+(876600*60*60+8640184.812866)*Tu+.093104*Tu^2-6.2*10^(-6)*Tu^3)/240;
thetaG0=deg2rad(thetaG0);
M0=M;
tf=3*60*60;
dt=60;
j=1;
while dt<tf
    M=M0+n*dt;
    [Reci,Veci]=COE2RVM(a,e,i,O,w,M,mu);
    thetaG=thetaG0+we*dt;
    Recef=eci2ecef(Reci,thetaG);
    long=atan2(Recef(2),Recef(1));
    lat=atan2(Recef(3),Recef(2)/sin(long));
    lattitude(j)=rad2deg(lat);
    longitude(j)=rad2deg(long);
    dt=dt+60;
    j=j+1;
end
figure
load worldmap2384.dat;
x= worldmap2384(:,1);
y= worldmap2384(:,2);
hold on
plot (x,y)
scatter(longitude,lattitude)
title('Homework 6 Problem 1 ISS Ground Track')
hold off