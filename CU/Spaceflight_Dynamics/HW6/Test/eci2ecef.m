function [Recef] = eci2ecef(Reci,thetaGST)
th=thetaGST;
R3th=[cos(th),sin(th),0;-sin(th),cos(th),0;0,0,1];
%%rotating about the 3rd axis to location Greenwich is at the time with respect
%%to the inertial frame.
Recef=R3th*Reci; %rotating the inertial postion