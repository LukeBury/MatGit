clear
clc


Mrod = .011;    %kg
r1   = .002375; %m
r2   = .003175; %m
Lrod = .409575; %m
 
Mex  = .0005;   %kg
Lex  = .01905;  %m

Dex  = .22384;  %m

Msq  = .0002;   %kg
b    = .010922; %m

Dsq  = .24517;  %m






%Moment of inertia for rod centered about fulcrum (kg*m^2)
Irod = (1/12)*(Mrod)*(3*(r2^2+r1^2)+Lrod^2);

%Calculate volume ratio of bar for hollow part to solid part of excess rod
v1 = pi*(r1^2)*Lex;
v2 = pi*(r2^2)*Lex;
Vratio = v1/v2;
SolidPortion = 1-Vratio;

%If excess rod weren't hollow, its mass would be...
MexTotal = Mex/SolidPortion;

%Weight of filled in hollow part of excess rod...
Mhollow = Vratio * MexTotal;


%Moment of inertia for excess rod length as defined by parallel axis
%theorem (kg*m^2)

I_exrod = (1/3)*MexTotal*(Lex^2)-(1/3)*Mhollow*(Lex^2);
PA_exrod = Mex*(Dex^2);

%Moment of inertia for garolite square (kg*m^2)
I_sq = (1/12)*Msq*b^2;
PA_sq = Msq*(Dsq^2);

%Total moment of inertia (kg*m^2)
I = Irod + I_exrod + PA_exrod + I_sq + PA_sq;


fprintf('The moment of inertia is %d kg*m^2\n', I);

