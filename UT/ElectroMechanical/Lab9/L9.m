clc
close all

t=[0:1/1000:1-1/1000];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% A
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wA1=2*pi/.07;
wA2=2*pi/.04;
wA3=2*pi/.032;
wA4=2*pi/.027;
wA5=2*pi/.024;
wA6=2*pi/.023;
wA7=2*pi/.021;
wA8=2*pi/.02;

wA=[wA1 wA2 wA3 wA4 wA5 wA6 wA7 wA8];

vA1=1.007;
vA2=2.036;
vA3=3.003;
vA4=4.023;
vA5=5.070;
vA6=6.07;
vA7=7.00;
vA8=7.99;

vA=[vA1 vA2 vA3 vA4 vA5 vA6 vA7 vA8];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% B
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wB1=2*pi/.077;
rtB1=.005;%%%%%%%%%%%%%%%%%%%%%% precision tuncates at .001
wB2=2*pi/.042;
rtB2=.002;
wB3=2*pi/.033;
rtB3=.002;
wB4=2*pi/.028;
rtB4=.002;
wB5=2*pi/.025;
rtB5=.001;
wB6=2*pi/.023;
rtB6=.001;
wB7=2*pi/.022;
rtB7=.001;
wB8=2*pi/.021;
rtB8=.001;

wB=[wB1 wB2 wB3 wB4 wB5 wB6 wB7 wB8];
rtB=[rtB1 rtB2 rtB3 rtB4 rtB5 rtB6 rtB7 rtB8];

vB1=1.038;
vB2=2.098;
vB3=3.024;
vB4=3.970;
vB5=4.99;
vB6=5.99;
vB7=7.08;
vB8=8.07;

vB=[vB1 vB2 vB3 vB4 vB5 vB6 vB7 vB8];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% C
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vC=2.144;
wC=2*pi/.039;
rtC=.001;%almost .002

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(t,D8)

wD1=2*pi/.07;
wD2=2*pi/.041;
wD3=2*pi/.031;
wD4=2*pi/.027;
wD5=2*pi/.025;
wD6=2*pi/.023;
wD7=2*pi/.021;
wD8=2*pi/.021;

wD=[wD1 wD2 wD3 wD4 wD5 wD6 wD7 wD8];

rtD1=.001;
rtD2=.001;
rtD3=.001;
rtD4=.001;
rtD5=.001;
rtD6=.001;
rtD7=.001;
rtD8=.001;

rtD=[rtD1 rtD2 rtD3 rtD4 rtD5 rtD6 rtD7 rtD8];

vD1=1.026;
vD2=2.036;
vD3=2.996;
vD4=4.065;
vD5=5.025;
vD6=5.97;
vD7=7.06;
vD8=8.05;

vD=[vD1 vD2 vD3 vD4 vD5 vD6 vD7 vD8];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plot(t,A1)
figure
plot(t,A8)
figure
plot(t,B1)
figure
plot(t,B8)
figure
plot(t,D1)
figure
plot(t,D8)
% figure(1)
% plot(vA,wA,'-o','markersize',10,'linewidth',2)
% 
% figure(2)
% plot(vB,wB,'-x','markersize',10,'linewidth',2)
% 
% figure(3)
% plot(vD,wD,'-x','markersize',10,'linewidth',2)
% 
% figure(4)
% plot(vB,rtB,'-x','markersize',10,'linewidth',2)
% 
% figure(5)
% plot(vD,rtD,'-x','markersize',10,'linewidth',2)









