clear
clc

% ------------------------------------------------------------------------
%%% 2b
% ------------------------------------------------------------------------
syms x y z u a J2 J3

acc = -(3*u*a*a*J2*z*z)/(2*((x^2+y^2+z^2)^(5/2))) + (u*a*a*J2)/(2*((x^2+y^2+z^2)^(3/2))) - (5*u*a*a*a*J3*z*z*z)/(2*((x^2+y^2+z^2)^(7/2))) + (3*u*a*a*a*J3*z)/(2*((x^2+y^2+z^2)^(5/2)));

dx = diff(acc,x);

dy = diff(acc,y);

dz = diff(acc,z);

a = [dx; dy; dz];

dadx = simplify(diff(a,x));
dady = simplify(diff(a,y));
dadz = simplify(diff(a,z));
dadu = simplify(diff(a,u));
dadJ2 = simplify(diff(a,J2));
dadJ3 = simplify(diff(a,J3));

fprintf('dadx \nx-hat: %s\ny-hat: %s\nz-hat: %s\n\n',dadx(1),dadx(2),dadx(3))
fprintf('dady \nx-hat: %s\ny-hat: %s\nz-hat: %s\n\n',dady(1),dady(2),dady(3))
fprintf('dadz \nx-hat: %s\ny-hat: %s\nz-hat: %s\n\n',dadz(1),dadz(2),dadz(3))
fprintf('dadu \nx-hat: %s\ny-hat: %s\nz-hat: %s\n\n',dadu(1),dadu(2),dadu(3))
fprintf('dadJ2 \nx-hat: %s\ny-hat: %s\nz-hat: %s\n\n',dadJ2(1),dadJ2(2),dadJ2(3))
fprintf('dadJ3 \nx-hat: %s\ny-hat: %s\nz-hat: %s\n\n',dadJ3(1),dadJ3(2),dadJ3(3))
clear x y z u a J2 J3







