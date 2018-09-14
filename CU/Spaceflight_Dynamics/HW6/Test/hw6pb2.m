clear

mu=398600.4415; %%km^3/s^2
i=deg2rad(51.6441);
O=deg2rad(225.6859);
e=.0006691;
w=deg2rad(42.3347);
M0=deg2rad(63.3332);
n=15.54051996*2*pi/24/60/60;%%rev per day to rad per sec
a=(mu/n^2)^(1/3);%km

latB=deg2rad(40.01);%given latitude of topocentric frame
lonB=deg2rad(254.83);%given longitude of topocentric frame
altB=1.615; %given altidue of topocentric frame

JD=JulianD(10,11,2016,1,0);
JDi=JulianD(10,5,2016,22,6.144384);
tf=(JD-JDi)*60*24*60;
M=M0+n*tf;

Tu=(JD-2451545)/36525;
we=7.2921158553*10^(-5);
%thetaG0=100.460618+36000.7705361*Tu+.00038793*Tu^2-2.6*10^(-8)*Tu^3;
thetaG0=(67310.54841+(876600*60*60+8640184.812866)*Tu+.093104*Tu^2-6.2*10^(-6)*Tu^3)/240;
thetaG0=deg2rad(thetaG0);

M0=M;
tf=3*60*60;
dt=1;
k=1;
while dt<tf
    M=M0+n*dt;
    [Reci,Veci]=COE2RVM(a,e,i,O,w,M,mu);
    thetaG=thetaG0+we*dt;
    Recef=eci2ecef(Reci,thetaG);
    [Rtopo]=ecef2topo(Recef,latB,lonB,altB);
   
    el=atan2(Rtopo(3),sqrt(Rtopo(1)^2+Rtopo(2)^2));%solving for the elevation
    el=rad2deg(el);%putting it in degrees
    if el<=0&&k>1
        break
    end
    if el>0
        B=atan2(Rtopo(2),-Rtopo(1));%solving for the azimuth
        B=rad2deg(B);%putting it in degrees
        if B<0
            B=B+360;%making sure matlab gives possitive coordinates
        end  
        azimuth(k)=B;
        zenithel(k)=90-el;
        t(k)=dt;
        k=k+1;
    end
    
    if dt == 50 * 60
        el
    end
    
    dt=dt+1;
end
%els'
figure
azimuth=deg2rad(azimuth);
polar(azimuth,zenithel)
title('Homework 6 Problem 2 ISS Pass')