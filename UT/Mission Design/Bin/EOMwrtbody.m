function dy =EOMwrtbody(t,y,mu)

r(1:3,1) = y(1:3);
v(1:3,1) = y(4:6);

rmag = norm(r);
%    Xdot

dy(1:3,1) = y(4:6);
dy(4,1) = -mu/rmag^3*r(1,1);
dy(5,1) = -mu/rmag^3*r(2,1);
dy(6,1) = -mu/rmag^3*r(3,1);

