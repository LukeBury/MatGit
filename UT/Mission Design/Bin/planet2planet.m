function [rPd,vPd,rPa,vPa,vd,va,exitflag]=planet2planet(type,eventd,dplanet,JDd,rnbd,...
    coepd,eventa,aplanet,JDa,rnba,coepa,rev,grade)
% Determines the departure and arrival states of the spacecraft and planets
% assuming patched conics method

% eventd: Defines if departure node is a planet or not
% eventa: Defines if arrival node is a planet or not
% dplanet: Departure planet number, found in "pheph.M"
% aplanet: Arrival planet number, found in "pheph.M"
% JDd: Departure julian date
% JDa: Arrival julian date
% rnbd: Departure non-body position vector (typically non-zero if used)
% rnba: Arrival non-body position vector (typically non-zero if used)
% grade: prograde (0) or retrograde orbit (1)
% vd: Velocity required for spacecraft departure [km/s]
% va: Velocity at spacecraft arrival [km/s]
AU=149597870.691;   %km
mu = 1.32712440017987E11;      %km3/s2 
eps= 23.439291111; % Obliquity of ecliptic at J2000
dxtol=1E-7;
tof = (JDa-JDd)*86400;   %time of flight [sec]

if dplanet==0
    if (eventd==40)||(eventd==41)
    rPd(1:3,1)=rnbd.*AU;
    vPd(1:3,1)=0;
    elseif (eventd==17)||(eventd==27)
%     coepd(6)=kepdt(coepd(1:7),JDd,mu);   
    [rPd, vPd]=orbel2rv(mu,coepd(1),coepd(2),coepd(3),coepd(4),coepd(5),coepd(6));
    [rPd,vPd]=kepuv(mu,(JDd-coepd(7))*86400,rPd,vPd,0,dxtol);
    end
else
    
    switch type
        case 1
    [rPd,vPd]=pleph(JDd,dplanet,11,1);
    Rot1 = [1 0 0;
       0 cosd(eps) sind(eps);
       0 -sind(eps) cosd(eps)];
    %Convert from Mean Equator and Mean Equinox of J2000 to Mean Ecliptic and
    %Mean Equinox of J2000 
    rPd=Rot1*rPd;
    vPd=Rot1*vPd;
   
        case 2
    [rPd,vPd]=pleph2o(dplanet,JDd);
     	
	    case 3
    [rPd,vPd]=pleph3o(dplanet,JDd);
	
        otherwise
    disp('Ephemeris type not available')
    end
end

if aplanet==0
    if (eventa==40)||(eventa==41)
    rPa(1:3,1)=rnba.*AU;    
    vPa(1:3,1)=0;
    elseif (eventa==17)||(eventa==27)
%     coepa(6)=kepdt(coepa(1:7),JDa,mu);
    [rPa, vPa]=orbel2rv(mu,coepa(1),coepa(2),coepa(3),coepa(4),coepa(5),coepa(6));    
    [rPa,vPa]=kepuv(mu,(JDa-coepa(7))*86400,rPa,vPa,0,dxtol);
    end
else
    switch type
        case 1
    [rPa,vPa]=pleph(JDa,aplanet,11,1);
    Rot1 = [1 0 0;
       0 cosd(eps) sind(eps);
       0 -sind(eps) cosd(eps)];
    %Convert from Mean Equator and Mean Equinox of J2000 to Mean Ecliptic and
    %Mean Equinox of J2000 
    rPa=Rot1*rPa;
    vPa=Rot1*vPa;
        case 2
    [rPa,vPa]=pleph2o(aplanet,JDa);
	    case 3
    [rPa,vPa]=pleph3o(aplanet,JDa);
        otherwise
    disp('!! Ephemeris type not available !!')
    end

end

%Hybrid lambert targeter with Dr. Izzo (ESA) fast convergence solution and 
%Lancaster & Blancard's robust solution 
 [vd,va,exitflag]=lambert(rPd,rPa,tof,rev,grade,mu);
 if (exitflag ~=1)||any(isnan(vd))||any(isnan(va))
     if exitflag==1 %Catch any NaN cases that exitflag didn't
     exitflag=-1;
     end
     disp(['!! Lambert Targeter failed between planet ',num2str(dplanet),' and planet ',num2str(aplanet),'!!'])
 end
