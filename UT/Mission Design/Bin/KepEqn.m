function [anomaly]=KepEqn(mu,M,e,p,dt)

%Determines respective anomaly (eccentric, parabolic, hyperbolic) by solving 
%kepler's equation

% INPUTS
% mu: gravitational parameter [km3/s2 or ER3/TU2]
% M: mean anomaly: [rad]
% e: eccentricity
% p: semiparameter [km]
% dt: change in time [sec or TU]
% tol:  tolerance of convergence
tol=1e-8;

% OUTPUTS
% anomaly (rad)

%Coded by Marty Brennan Jan. 2013

%If not parabolic trajectory, set p=0 and dt=0 
%If parabolic trajectory use p and t

M=mod(M,2*pi);
if e < 1.0  %for elliptic orbits

    if (((M>-pi)&&(M<0))||(M>pi))
        E=M-e; %(rad)
    else
        E=M+e; %(rad)
    end
    
    E0 = 0; %initialize E0
    E1 = E;
    while abs(E1-E0)>tol
        E0 = E1; %(rad)
        E1 = E0+(M-E0+e*sin(E0))/(1-e*cos(E0)); %(rad)
    end
   anomaly = E1; %(rad)
   
elseif e >1.0 %for hyperbolic orbits
      if e<1.6
        if (((M>-pi)&&(M<0))||(M>pi))
            H=M-e; %(rad)
        else
            H=M+e; %(rad)
        end 
    else
        if (e<3.6)&&(abs(M)>pi)
            H=M-sign(M)*e; %(rad)
        else
            H01=M/(e-1); %(rad)
            H02=M^2/(e*(e-1)*sinh(M/(e-1))-M); %(rad)
            H=(H01+H02)/2; %(rad)
        end 
    end
    
    H0 = 0;
    H1 = H; %(rad)
    while abs(H1-H0)>tol
        H0 = H1; %(rad)
        H1 = H0+(M+H0-e*sinh(H0))/(e*cosh(H0)-1); %(rad)
    end
   anomaly = H1; %(rad)
   
elseif e==1.0 %for parabolic orbits
    np=2*sqrt(mu/p^3);
    s=acot(3/2*np*dt)/2;
    w=atan((tan(s))^(1/3));
    B=2*cot(2*w);
    anomaly = B; %(rad)
end

















