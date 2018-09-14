function [r2,v2]=kepuv(mu,dt,r0,v0,disp,dxtol)

%INPUTS:
%mu: gravitational parameter of body (DU^3/TU^2)
%dt: time of flight (TU)
%r0: initial radius col vec (DU)
%v0: initial velocity col vec (DU/TU)
%disp:  display iterations (0:no, 1:yes)
%dxtol: convergence tolerance for universal kepler eqn

%OUTPUTS:
%r: final radius col vec (DU)
%v: final velocity col vec (DU/TU)

% AU=149597870; 
ztol= 0.01; % z tolerance for defining parabolic-case buffer zone
% alphatol = 1/(1000*AU); %alpha tolerance for defining parabolic-case buffer zone
etol = 0.01;
dx=Inf;
r0mag=sqrt(r0'*r0);
v0mag=sqrt(v0'*v0);


alpha=-v0mag^2/mu+2/r0mag; % =1/a from E
% a=r0mag/(2-r0mag*v0mag^2/mu);
% e=sqrt(1-norm(cross(r0,v0))^2*alpha/mu)
e=((v0mag^2*r0-(r0'*v0)*v0)/mu-r0/r0mag);
e=sqrt(e'*e);
sqrtmu=sqrt(mu);

% if alpha >alphatol  %**replaced alpha check with e, for nondimensionality
if any(dt) %check if dt=0
if e < (1-etol)
    x0=sqrtmu*dt*alpha;
% elseif alpha<-alphatol
elseif e > (1+etol)
    x0=sign(dt)/sqrt(-alpha)*log(-2*mu*alpha*dt/((r0'*v0)+sign(dt)*...
        sqrt(-mu/alpha)*(1-r0mag*alpha)));
% elseif abs(alpha)<alphatol
elseif abs(e-1) < etol
    h=cross(r0,v0);
    p=norm(h)^2/mu;
    s=0.5*(pi/2-atan(3*dt*sqrt(mu/p^3)));
%     s=0.5*acot(3*dt*sqrt(mu/p^3)); %omit due do acot() 
    w=atan((tan(s))^(1/3));
    x0=2*sqrt(p)*cot(2*w);
%     x0=0;
end
% i=0;
x02=x0*x0;
z=x02*alpha;

% if disp==1
% fprintf('  Iter           x                      z\n')
% fprintf('%4.f %15.16e %15.16e \n', i,x0,z)
% end

while abs(dx)>dxtol
    z=x02*alpha;
% i=i+1 ;
%choose correct form of C and S for conic type
if z>ztol
%     0
    sqrtz=sqrt(z);
    C=(1-cos(sqrtz))/z;
    S=(sqrtz-sin(sqrtz))/sqrt(z^3);
elseif z<-ztol
%     2
    sqrtnz=sqrt(-z);
    C=(1-cosh(sqrtnz))/z;
    S=(sinh(sqrtnz)-sqrtnz)/sqrt((-z)^3);
elseif abs(z)<=ztol
%     1
    C=1/2-z/24+z^2/720-z^3/40320+z^4/3628800;       %1/2!-z/4!+z^2/6!-z^3/8!+z^4/10!
    S=1/6-z/120+z^2/5040-z^3/362880+z^4/39916800;   %1/3!-z/5!+z^2/7!-z^3/9!+z^4/11!
end
    rx=x02*C+x0*(r0'*v0)*(1-z*S)/sqrtmu+r0mag*(1-z*C);
    dx=(sqrtmu*dt-x02*x0*S-x02*(r0'*v0)*C/sqrtmu-r0mag*x0*(1-z*S))/rx;
    
    x0=x0+dx;
    x02=x0*x0;
% if disp==1
% fprintf('%4.i %15.16e %15.16e \n', i,x0,z)
% end

end

f=1-x02/r0mag*C;
g=dt-x02*x0/sqrtmu*S;
r2 = f*r0 + g*v0;
r2mag=norm(r2);

fdot=sqrtmu/(r0mag*r2mag)*x0*(z*S-1);
gdot=1-x02/r2mag*C;
v2 = fdot*r0 + gdot*v0;

% % fgcheck=f*gdot-fdot*g
% if disp==1
% %check that f*gdot-fdot*g=1
% fgcheck=f*gdot-fdot*g
% a0=a
% alpha=-norm(v2)^2/mu+2/r2mag; % =1/a from E
% af=1/alpha
% end

else
    r2=r0;
    v2=v0;
end

end


