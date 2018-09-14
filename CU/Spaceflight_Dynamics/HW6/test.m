clear
close all
clc
rad2deg = 180/pi;
r = 1;
i = 1;
for theta = 0:.001:2*pi
    x(i) = r*cos(theta);
    y(i) = r*sin(theta);
    i = i +1;
end
plot(x,y)
figure
plot(x,atan2(y,x)*rad2deg)

atan2(1,1)*rad2deg
atan2(1,-1)*rad2deg
atan2(-1,-1)*rad2deg
atan2(-1,1)*rad2deg