
clc
clear all
close all

AU=149597870;   %(km) provides generally accepted AU measurement for convenience 

mu(1) = 2.2032e+004;        % all MU: km3/s2
mu(2) = 3.24858596e+005;      
mu(3) = 398600.4418;  
mu(4) = 4.28283758e+004;     
mu(5) = 1.26712764e+008;  
mu(6) = 3.7940585e+007;   
mu(7) = 5.794549e+006;    
mu(8) = 6.836527e+006;    
mu(9) = 9.7178e+002;        
mu(10) = 1.32712442099E11; 



r(1)=0.38709831;   % all mean distances of planets from sun: AU
r(2)=0.72332982;
r(3)=1.000001018;
r(4)=1.523679342;
r(5)=5.202603191;
r(6)=9.554909596;
r(7)=19.21844606;
r(8)=30.11038687;
r(9)=39.544674;
r(10)=696000/AU;

for i=1:9
rsoi(i)=r(i)*(mu(i)/mu(10))^(2/5);
end


