function [ a,e,i,omega,w, nu  ] = co2el(R,V,u)

r=norm(R);
v=norm(V);

energy=v^2.0/2.0-u/r;

a=-u/(2.0*energy);

e_vec=(1.0/u)*[(v^2.0-u/r).*R-dot(R,V).*V];

e=norm(e_vec);

h_vec=cross(R,V);

h=norm(h_vec);

%i=acos(h_vec(3)/h)*180/pi;
i=atan2(norm(cross(h_vec,[0 0 1])),dot(h_vec,[0 0 1]))*180/pi;

if abs(i)>1e-12

n_vec=cross([0 0 1],h_vec);

n=norm(n_vec);

%omega=acos(n_vec(1)/n)*180/pi;
omega=atan2(norm(cross(n_vec,[1 0 0])),dot(n_vec,[1 0 0]))*180/pi;

if (n_vec(2)<0.0)
    omega=360.0-omega;
end

if abs(e)<1e-10
    e_vec=[1 0 0]';
end

%w=acos(dot(n_vec,e_vec)/(e*n))*180/pi;
w=atan2(norm(cross(n_vec,e_vec)),dot(n_vec,e_vec))*180/pi;

if (e_vec(3)<0)
    w=360.0-w;
end

else
   
    n_vec=[1 0 0];

n=norm(n_vec);

omega=0;

if abs(e)<1e-10
    e_vec=[1 0 0]';
end

w=atan2(e_vec(2),e_vec(1))*180/pi;

if (h_vec(3)<0)
    w=360.0-w;
end

if w>360
    w=w-360;
end

if w<0
    w=w+360;
end

end

%nu=acos(dot(e_vec,R)/(e*r))*180/pi;
nu=atan2(norm(cross(e_vec,R)),dot(e_vec,R))*180/pi;

if (dot(R,V)<0)
   nu=360-nu;
end

end