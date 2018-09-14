clear 
clc

% Sun's gravitational parameter in km^3/s^2
mu_Sun=1.32712428E11;

% Planet orbital radii in km
r_Earth=149597870;
r_Venus=108208600.379;
r_Mars=227939183.827;
r_Jupiter=778298355.829;

% Hohmann transfer for Earth to Venus
atrans2Venus=(r_Earth+r_Venus)/2;
vtransa1=sqrt(((2*mu_Sun)/r_Earth)-(mu_Sun/atrans2Venus));
vi_Earth=sqrt(mu_Sun/r_Earth);
vtransb1=sqrt(((2*mu_Sun)/r_Venus)-(mu_Sun/atrans2Venus));
vf_Venus=sqrt(mu_Sun/r_Venus);
deltaVa1=vtransa1-vi_Earth;
deltaVb1=vf_Venus-vtransb1;
deltaVtot1=abs(deltaVa1)+abs(deltaVb1);
% Transfer duration
Ttrans1=pi*sqrt((atrans2Venus^3)/mu_Sun);
Tt1=Ttrans1/(60*60*24*365.25);


%Hohmann transfer for Earth to Mars
atrans2Mars=(r_Earth+r_Mars)/2;
vtransa2=sqrt(((2*mu_Sun)/r_Earth)-(mu_Sun/atrans2Mars));
vi_Earth=sqrt(mu_Sun/r_Earth);
vtransb2=sqrt(((2*mu_Sun)/r_Mars)-(mu_Sun/atrans2Mars));
vf_Mars=sqrt(mu_Sun/r_Mars);
deltaVa2=vtransa2-vi_Earth;
deltaVb2=vf_Mars-vtransb2;
deltaVtot2=abs(deltaVa2)+abs(deltaVb2);
% Transfer duration
Ttrans2=pi*sqrt((atrans2Mars^3)/mu_Sun);
Tt2=Ttrans2/(60*60*24*365.25);
% Phase angle
omega_tgt=sqrt(mu_Sun/((r_Mars)^3));
omega_int=sqrt(mu_Sun/(r_Earth)^3);
alpha_L=omega_tgt*Ttrans2;
phase_angle=alpha_L-pi;
synod=2*pi/(omega_int-omega_tgt);
synod_yrs=synod/(60*60*24*365.25);


%Hohmann transfer from Earth to Jupiter
atrans2Jupiter=(r_Earth+r_Jupiter)/2;
vtransa3=sqrt(((2*mu_Sun)/r_Earth)-(mu_Sun/atrans2Jupiter));
vi_Earth=sqrt(mu_Sun/r_Earth);
vtransb3=sqrt(((2*mu_Sun)/r_Jupiter)-(mu_Sun/atrans2Jupiter));
vf_Jupiter=sqrt(mu_Sun/r_Jupiter);
deltaVa3=vtransa3-vi_Earth;
deltaVb3=vf_Jupiter-vtransb3;
deltaVtot3=abs(deltaVa3)+abs(deltaVb3);
% Transfer duration
Ttrans3=pi*sqrt((atrans2Jupiter^3)/mu_Sun);
Tt3=Ttrans3/(60*60*24*365.25);

Earth2Venus=table(r_Venus,deltaVa1,deltaVb1,deltaVtot1,Tt1);
Earth2Mars=table(r_Mars,deltaVa2,deltaVb2,deltaVtot2,Tt2,phase_angle,synod_yrs);
Earth2Jupiter=table(r_Venus,deltaVa3,deltaVb3,deltaVtot3,Tt3);

% Problem #3
mu_Earth=3.986E5; %km^3/s^2
h=6000; %km
a_tgt=h+6378.1363; %km
omega_trgt=sqrt(mu_Earth/(a_tgt^3)); %rad/s
k_tgt=1; %number of orbits
k_int=1;
theta_phase=30*(pi/180); %radians
tao_phase=((2*pi*k_tgt)+theta_phase)/omega_trgt; %seconds
a_phase=(mu_Earth*(tao_phase/(2*pi*k_int))^2)^(1/3); %km
DeltaVtot=2*abs(sqrt(((2*mu_Earth)/a_tgt)-(mu_Earth/a_phase))-sqrt(mu_Earth/a_tgt)); %km/s
tao_trans=2*pi*sqrt((a_phase^3)/mu_Earth); 
transtime_hours=tao_trans/(60*60);
ans=table(a_tgt, omega_trgt, a_phase, tao_trans) 

clear
clc

% Rendevous between ISS and Progress supply vessel.
% Givens:
h_trgt=400; %km
r_trgt=6378.1363+h_trgt; %km
mu_Earth=3.986E5; %km^3/s^2
w=sqrt(mu_Earth/(r_trgt^3)); %rad/s


% Initial position and velocity vectors.
r_init=[200 300 50]; %m
x0=r_init(1);
y0=r_init(2);
z0=r_init(3);
v_init=[0.1 0.0 -0.1]; %m/s
vx0=v_init(1);
vy0=v_init(2)
vz0=v_init(3);

% Position and velocity to begin 5 minute transfer.
t1=5*60; %seconds

% This computation shows the velocity necessary to get the vehicle to meet
% the target. 
vy0=((6*x0*(w*t1-sin(w*t1))-y0)*w*sin(w*t1)-2*w*x0*(4-3*cos(w*t1))*(1-cos(w*t1)))/...
    ((4*sin(w*t1)-3*w*t1)*sin(w*t1)+4*(1-cos(w*t1))^2);
vx0=-((w*x0)*(4-3*cos(w*t1))+2*(1-cos(w*t1))*vy0)/sin(w*t1)
vz0=-z0*w*cot(w*t1);
ans = table(w,x0,t1,vy0)
% This is the deltaV
Vmeet=[vx0 vy0 vz0]
DeltaV_1=v_init-Vmeet;
Delta_V_1=norm(DeltaV_1)

% Now compute the second burn to get r and v to zero.
xt=(vx0*cos(w*t1))+(2*vy0+3*w*x0)*sin(w*t1);
yt=(4*vy0+6*w*x0)*cos(w*t1)-(2*vx0*sin(w*t1))-(6*w*x0+3*vy0);
zt=(-z0*w*sin(w*t1))+((vz0*cos(w*t1)));
ans = table(vx0,w,t1,x0,vy0)
% This is the velocity change necessary to get the vehicle to stop.
DeltaStop=[xt yt zt]
Delta_V_2=norm(DeltaStop)

% This is the total deltaV for the entire maneuver.
Delta_V_tot=Delta_V_1+Delta_V_2
