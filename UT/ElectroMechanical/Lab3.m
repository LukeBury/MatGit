clc
clf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Problem 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Convert voltage to displacement
d = @(V) 11.032*V;

% How much higher right side is naturally
done = d(aaR)-d(aaL);

% How much each side dropped with weight
dL = d(aaL)-d(abL);
dR = d(aaR)-d(abR);

% How much more the right side dropped than the left
dtwo = dR-dL;

% Initial, weightless angle
angle1 = atand(done/38.1);

% Final angle w/ weight
angle2 = atand((done-dtwo)/38.1);

% How much the weight at shear center twisted the setup (want 0)
deltaangle = angle1-angle2;

% Dammit


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Problem 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Constants
% Applied Force
W = .250*9.8; % Newtons

C = 38.1; % mm

% Chordwise Coordinate of weight
yf = 0; % mm

% Shear Center
e = 17.07; % +- .005 mm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Build Twist Matrix
Twist = zeros(5,5);


% Column 1 (c)
t1 = d(bccR)-d(bccL);
Twist(1,1) = atand(t1/C);


% Column 2 (e)
t2 = d(becR)-d(becL);
Twist(1,2) = atand(t2/C);

t3 = d(beeR)-d(beeL);
Twist(2,2) = atand(t3/C);


% Column 3 (g)
t4 = d(bgcR)-d(bgcL);
Twist(1,3) = atand(t4/C);

t5 = d(bgeR)-d(bgeL);
Twist(2,3) = atand(t5/C);

t6 = d(bggR)-d(bggL);
Twist(3,3) = atand(t6/C);


% Column 4 (i)
t7 = d(bicR)-d(bicL);
Twist(1,4) = atand(t7/C);

t8 = d(bieR)-d(bieL);
Twist(2,4) = atand(t8/C);

t9 = d(bigR)-d(bigL);
Twist(3,4) = atand(t9/C);

t10 = d(biiR)-d(biiL);
Twist(4,4) = atand(t10/C);

% Column 5 (k)

t11 = d(bkcR)-d(bkcL);
Twist(1,5) = atand(t11/C);

t12 = d(bkeR)-d(bkeL);
Twist(2,5) = atand(t12/C);

t13 = d(bkgR)-d(bkgL);
Twist(3,5) = atand(t13/C);

t14 = d(bkiR)-d(bkiL);
Twist(4,5) = atand(t14/C);

t15 = d(bkkR)-d(bkkL);
Twist(5,5) = atand(t15/C);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Build Deflection Matrix

%%% Find averaged positions of L&R data
defc = zeros(1,1);
defc(1) = (d(bccL)+d(bccR))/2000; %bcc

defe = zeros(2,1);
defe(1) = (d(becL)+d(becR))/2000; %bec
defe(2) = (d(beeL)+d(beeR))/2000; %bee

defg = zeros(3,1);
defg(1) = (d(bgcL)+d(bgcR))/2000; %bgc
defg(2) = (d(bgeL)+d(bgeR))/2000; %bge
defg(3) = (d(bggL)+d(bggR))/2000; %bgg

defi = zeros(4,1);
defi(1) = (d(bicL)+d(bicR))/2000; %bic
defi(2) = (d(bieL)+d(bieR))/2000; %bie
defi(3) = (d(bigL)+d(bigR))/2000; %big
defi(4) = (d(biiL)+d(biiR))/2000; %bii

defk = zeros(5,1);
defk(1) = (d(bkcL)+d(bkcR))/2000; %bkc
defk(2) = (d(bkeL)+d(bkeR))/2000; %bke
defk(3) = (d(bkgL)+d(bkgR))/2000; %bkg
defk(4) = (d(bkiL)+d(bkiR))/2000; %bki
defk(5) = (d(bkkL)+d(bkkR))/2000; %bkk




%%% Make spanwise coordinates
sw = [139.7; 279.4; 419.1; 558.8; 698.5]./1000; %m
% 
% hold all % Plotting measured deflection
% plot(sw(1),(defc-.0024)/(defc-.0024),'-o','linewidth',2)
% plot(sw(1:2),(defe-.0024)/(defc-.0024),'-o','linewidth',2)
% plot(sw(1:3),(defg-.0024)/(defc-.0024),'-o','linewidth',2)
% plot(sw(1:4),(defi-.0024)/(defc-.0024),'-o','linewidth',2)
% plot(sw(1:5),(defk-.0024)/(defc-.0024),'-o','linewidth',2)
% xlabel('Spanwise Coordinate (m)')
% ylabel('Normalized Bending Deflection')
% title('Measured Bending Deflection Trend Based on Spanwise Shear Force Location')
% legend('@ .1397 m','@ .2794 m', '@ .4191 m','@ .5588 m','@ .6985 m')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Theoretical Calculations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%% Bending Moment
Mc = W*sw(1);
Me = W*sw(2);
Mg = W*sw(3);
Mi = W*sw(4);
Mk = W*sw(5);

%%%%%%%%%%%%%%%% Torque
T = W*e;

%%%%%%%%%%%%%%%% Bending Deflection

wc = zeros(1,1);
wc(1,1) = Mc * (sw(1)*(sw(1)^2)/2 - (sw(1)^3)/6);

we = zeros(2,1);
we(1,1) = Me * (sw(2)*(sw(1)^2)/2 - (sw(1)^3)/6);
we(2,1) = Me * (sw(2)*(sw(2)^2)/2 - (sw(2)^3)/6);

wg = zeros(3,1);
wg(1,1) = Mg * (sw(3)*(sw(1)^2)/2 - (sw(1)^3)/6);
wg(2,1) = Mg * (sw(3)*(sw(2)^2)/2 - (sw(2)^3)/6);
wg(3,1) = Mg * (sw(3)*(sw(3)^2)/2 - (sw(3)^3)/6);

wi = zeros(4,1);
wi(1,1) = Mi * (sw(4)*(sw(1)^2)/2 - (sw(1)^3)/6);
wi(2,1) = Mi * (sw(4)*(sw(2)^2)/2 - (sw(2)^3)/6);
wi(3,1) = Mi * (sw(4)*(sw(3)^2)/2 - (sw(3)^3)/6);
wi(4,1) = Mi * (sw(4)*(sw(4)^2)/2 - (sw(4)^3)/6);

wk = zeros(5,1);
wk(1,1) = Mk * (sw(5)*(sw(1)^2)/2 - (sw(1)^3)/6);
wk(2,1) = Mk * (sw(5)*(sw(2)^2)/2 - (sw(2)^3)/6);
wk(3,1) = Mk * (sw(5)*(sw(3)^2)/2 - (sw(3)^3)/6);
wk(4,1) = Mk * (sw(5)*(sw(4)^2)/2 - (sw(4)^3)/6);
wk(5,1) = Mk * (sw(5)*(sw(5)^2)/2 - (sw(5)^3)/6);


% 
% figure % Theoretical deflection
% hold all
% plot(sw(1),(defc - wc + .1499)/(defc - wc + .1499),'-x','linewidth',2)
% plot(sw(1:2),(defc - we + .1499)/(defc - wc + .1499),'-x','linewidth',2)
% plot(sw(1:3),(defc - wg + .1499)/(defc - wc + .1499),'-x','linewidth',2)
% plot(sw(1:4),(defc - wi + .1499)/(defc - wc + .1499),'-x','linewidth',2)
% plot(sw(1:5),(defc - wk + .1499)/(defc - wc + .1499),'-x','linewidth',2)
% xlabel('Spanwise Coordinate (m)')
% ylabel('Normalized Bending Deflection')
% title('Theoretical Bending Deflection Trend Based on Spanwise Shear Force Location')
% legend('@ .1397 m','@ .2794 m', '@ .4191 m','@ .5588 m','@ .6985 m')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Do Tha Twist

TheoTwist = zeros(5,1);

TheoTwist(1) = T*sw(1);
TheoTwist(2) = T*sw(2);
TheoTwist(3) = T*sw(3);
TheoTwist(4) = T*sw(4);
TheoTwist(5) = T*sw(5);
% 
% figure %% Theoretical Twist - Normalized
% plot(sw,TheoTwist/TheoTwist(5),'linewidth',2)
% xlabel('Spanwise Coordinate (m)')
% ylabel('Normalized Sectional Twist')
% title('Theoretical Twist Trend')
% 
% figure %% Measured Twist Based on Spanwise Shear Force Location
% hold all
% plot(sw(1),Twist(1,1),'x','linewidth',2)
% plot(sw(1:2),Twist(1:2,2),'linewidth',2)
% plot(sw(1:3),Twist(1:3,3),'linewidth',2)
% plot(sw(1:4),Twist(1:4,4),'linewidth',2)
% plot(sw(1:5),Twist(1:5,5),'linewidth',2)
% xlabel('Spanwise Coordinate (m)')
% ylabel('Clockwise Sectional Twist (degrees)')
% title('Measured Twist Based on Spanwise Shear Force Location')
% legend('@ .1397 m','@ .2794 m', '@ .4191 m','@ .5588 m','@ .6985 m')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Calculate EI and GJ
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EI
%guessing initial height at .047 m

dispc = .047 - defc;
dispe = .047 - defe;
dispg = .047 - defg;
dispi = .047 - defi;
dispk = .047 - defk;


EIc = zeros(1,1);
EIc(1,1) = (Mc * (sw(1)*(sw(1)^2)/2 - (sw(1)^3)/6))/dispc;

EIe = zeros(2,1);
EIe(1,1) = (Me * (sw(2)*(sw(1)^2)/2 - (sw(1)^3)/6))/dispe(1);
EIe(2,1) = (Me * (sw(2)*(sw(2)^2)/2 - (sw(2)^3)/6))/dispe(2);

EIg = zeros(3,1);
EIg(1,1) = (Mg * (sw(3)*(sw(1)^2)/2 - (sw(1)^3)/6))/dispg(1);
EIg(2,1) = (Mg * (sw(3)*(sw(2)^2)/2 - (sw(2)^3)/6))/dispg(2);
EIg(3,1) = (Mg * (sw(3)*(sw(3)^2)/2 - (sw(3)^3)/6))/dispg(3);

EIi = zeros(4,1);
EIi(1,1) = (Mi * (sw(4)*(sw(1)^2)/2 - (sw(1)^3)/6))/dispi(1);
EIi(2,1) = (Mi * (sw(4)*(sw(2)^2)/2 - (sw(2)^3)/6))/dispi(2);
EIi(3,1) = (Mi * (sw(4)*(sw(3)^2)/2 - (sw(3)^3)/6))/dispi(3);
EIi(4,1) = (Mi * (sw(4)*(sw(4)^2)/2 - (sw(4)^3)/6))/dispi(4);

EIk = zeros(5,1);
EIk(1,1) = (Mk * (sw(5)*(sw(1)^2)/2 - (sw(1)^3)/6))/dispk(1);
EIk(2,1) = (Mk * (sw(5)*(sw(2)^2)/2 - (sw(2)^3)/6))/dispk(2);
EIk(3,1) = (Mk * (sw(5)*(sw(3)^2)/2 - (sw(3)^3)/6))/dispk(3);
EIk(4,1) = (Mk * (sw(5)*(sw(4)^2)/2 - (sw(4)^3)/6))/dispk(4);
EIk(5,1) = (Mk * (sw(5)*(sw(5)^2)/2 - (sw(5)^3)/6))/dispk(5);

EI = [EIc;EIe;EIg;EIi;EIk];

% histfit(EI)

%%% w/ std dev = 2.3401 +- 1.405824595

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GJ

%Switching to radians
Twistc = Twist(1,1)*pi/180;
Twiste = Twist(1:2,2)*pi/180;
Twistg = Twist(1:3,3)*pi/180;
Twisti = Twist(1:4,4)*pi/180;
Twistk = Twist(1:5,5)*pi/180;



GJc = zeros(1,1);
GJc(1) = T * sw(1) / Twistc(1);

GJe = zeros(2,1);
GJe(1) = T * sw(1) / Twiste(1);
GJe(2) = T * sw(2) / Twiste(2);

GJg = zeros(3,1);
GJg(1) = T * sw(1) / Twistg(1);
GJg(2) = T * sw(2) / Twistg(2);
GJg(3) = T * sw(3) / Twistg(3);

GJi = zeros(4,1);
GJi(1) = T * sw(1) / Twisti(1);
GJi(2) = T * sw(2) / Twisti(2);
GJi(3) = T * sw(3) / Twisti(3);
GJi(4) = T * sw(4) / Twisti(4);

GJk = zeros(5,1);
GJk(1) = T * sw(1) / Twistk(1);
GJk(2) = T * sw(2) / Twistk(2);
GJk(3) = T * sw(3) / Twistk(3);
GJk(4) = T * sw(4) / Twistk(4);
GJk(5) = T * sw(5) / Twistk(5);

GJ = [GJc;GJe;GJg;GJi;GJk];

histfit(GJ)

% w/ std dev =  291.03966 +- 205.6391532









% plot(sw(1),Twist(1,1),'x','linewidth',2)
% plot(sw(1:2),Twist(1:2,2),'linewidth',2)
% plot(sw(1:3),Twist(1:3,3),'linewidth',2)
% plot(sw(1:4),Twist(1:4,4),'linewidth',2)
% plot(sw(1:5),Twist(1:5,5),'linewidth',2)


















