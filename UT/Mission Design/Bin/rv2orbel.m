function [a,e,i,o,w,ta,Case] = rv2orbel(r,v,mu)

%Determines Classical orbital elements from state vectors of position and 
%velocity (XYZ-inertial)
%****NOTE****
% If Circular Equatorial: set w = 0, o = 0, and ta = lamdaTrue
% If Circular Inclined: set w = 0, and ta = u
% If Elliptical Equatorial: set RAAN = 0 and w = whatTrue

% r(3): position column vector
% v(3): velocity column vector
% mu: Gravitational parameter [km3/s2]
% a: Semimajor axis  [km]
% e: Eccentricity  [unitless]
% i: Inclination [deg]
% o: Right ascention of ascending node [deg], in equitorial plane
%       between vernal equinox (Ihat) and ascending node
% w: Argument of perigee [deg], between ascending node and periapsis
% ta: True Anomaly [deg]
% lamdaTrue: True longitude [deg], eastward from Ihat to the position of
%       satellite
% u: Argument of latitude [deg], between ascending node and the satellite's
%       position vector
% whatTrue: True longitude of periapsis [deg], between vernal equinox
% (Ihat) and eccentricity vector
% Case: if 0 then normal classical orbital elements have been determined
%       if 1 then special case of Circular Inclined orbit (ta=u)
%       if 2 then special case of Circular Equatorial orbit (ta=lamdaTrue)
%       if 3 then speical case of Elliptical Equatorial orbit (w=whatTrue)

%Code created by Martin Brennan
%Last edited Aug 1, 2014


%identify cut-off point for approximating zero (accounts for numerical
%error)
zero=10^-10;
rmag = sqrt(r'*r);
vmag = sqrt(v'*v);
h = [r(2)*v(3) - r(3)*v(2);... % cross product
     r(3)*v(1) - r(1)*v(3);...
     r(1)*v(2) - r(2)*v(1)];

hmag = sqrt(h'*h);
n= [-h(2);h(1);0];% cross product (zhat,h)
nmag = sqrt(n'*n);
evec = ((vmag*vmag-mu/rmag)*r-(r'*v)*v)/mu;
e = sqrt(evec'*evec);
E = vmag*vmag/2-mu/rmag;

if e==1.0       %parabolic orbit
%     p = hmag*hmag/mu;
    a = Inf;
else            %all other orbits
    a = -mu/(2*E);
%     p = a*(1-e*e);
end

i = acosd(h(3)/hmag);

if e>zero
    if (abs(i)<zero)||(180-abs(i)<zero) 
        % If Elliptical Equatorial:
        Case = 3;
        o=0;         % set o = 0 and w = whatTrue 
        whatTrue = acosd(evec(1)/e);
        
        if ((evec(2)<0)&&(abs(i)<zero))||((-evec(2)<0)&&(180-abs(i)<zero))
            whatTrue = 360 -whatTrue;
        end
        w=whatTrue;
        evecr=evec'*r/(e*rmag);

        if evecr>=1.0
            ta=0;
        elseif evecr<=-1.0
            ta=180;    
        else
        ta = acosd(evecr);
        end
        if (r'*v)<0
            ta = 360-ta;
        end
        ta=mod(ta,360);
    else            % Normal cases
        Case = 0;
        o = acosd(n(1)/nmag);
        if n(2)<0
            o = 360-o;
        end
        w = acosd((n'*evec)/(nmag*e));
        if evec(3)<0
            w = 360-w;
        end
        evecr=evec'*r/(e*rmag);
        if evecr>1.0
            ta=0;
        elseif evecr<-1.0
            ta=180;
        else
        ta = acosd(evecr);
        end
        if (r'*v)<0
            ta = 360-ta;
        end   
        ta=mod(ta,360);
    end 
else
     e=0;
    if ((abs(i)>zero))&&(180-abs(i)>zero)    % If Circular Inclined: 
        Case = 1;
        w=0;        % set w = 0, and TA = u
        if (n'*r)/(nmag*rmag)>1.0
            u=0;
        else
        u = acosd((n'*r)/(nmag*rmag));
        end
        if r(3)<0
            u = 360-u;
        end
        ta = u;
        ta=mod(ta,360);
        o = acosd(n(1)/nmag);
        if n(2)<0
            o = 360-o;
        end
    else
        % If Circular Equatorial: 
        Case = 2;
        w=0;        % set w = 0, o = 0, and TA = lamdaTrue
        o=0;
        lamdaTrue = acosd(r(1)/rmag);

        if (180-abs(i)<zero) 
            lamdaTrue = -lamdaTrue;
        end
        
        if (r(2)<0)
            lamdaTrue = 360-lamdaTrue;
        end
        ta = lamdaTrue;
        ta=mod(ta,360);
    end
end

