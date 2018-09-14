clc; clear

syms x y z v

r1 = sqrt((x+v)^2 + y^2 + z^2);
r2 = sqrt((x-1+v)^2 + y^2 + z^2);

% Newton Raphson
x_eqn = x - (1-v)*(x+v)/r1^3 - v*(x-1+v)/r2^3;
y_eqn = y - (1-v)*y/r1^3 - v*y/r2^3;

x_diff = diff(x_eqn, x);
y_diff = diff(y_eqn, y);

% STM Partials
O = (x^2 + y^2)/2 + (1-v)/r1 + v/r2;

Oxx = diff(diff(O,x),x)

Oyy = diff(diff(O,y),y)

Ozz = diff(diff(O,z),z)

Oxy = diff(diff(O,x),y)

Oyx = diff(diff(O,y),x)

Oxz = diff(diff(O,x),z)

Ozx = diff(diff(O,z),x)

Oyz = diff(diff(O,y),z)

Ozy = diff(diff(O,z),y)


