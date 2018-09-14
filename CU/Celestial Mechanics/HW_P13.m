clc

% syms m1 m2 m3 r12 r23 r13 G real
% assume(m1>0); assume(m1==m2); assume(m1==m3);
% M = m1+m2+m3;
% 
% C = (G*G/M)*((m1*m2/r12 + m2*m3/r23 + m1*m3/r13)^2)*(m1*m2*r12*r12 + m2*m3*r23*r23 + m1*m3*r13*r13);
% dr12 = simplify(diff(C,r12))
% 
% dr23 = simplify(diff(C,r23))
% 
% dr13 = simplify(diff(C,r13))
% 
% r13 = r12+r23;
% a = subs(dr13);
% a = simplify(a)
% solve(a,r23)


% syms m r12 r23 G real
% r13 = r12+r23;
% assume(r12>0); assume(r23>0); assume(G>0); assume(m>0);
% C = (1/3)*G*G*(m^5)*((1/r12 + 1/r23 + 1/r13)^2)*(r12^2 + r23^2 + r13^2)
% 
% C = subs(C);
% C = simplify(C)
% a = diff(C,r12);
% a = simplify(a)
% solve(a,r12)


syms m r12 r23 G real
r13 = r12+r23;
assume(r12>0); assume(r23>0); assume(G>0); assume(m>0);
C = (1/3)*G*G*(m^5)*((1/r12 + 1/r23 + 1/r13)^2)*(r12^2 + r23^2 + r13^2);

eqn1 = diff(C,r12)


eqn2 = diff(C,r23)
(G^2*m^5*(4*r12 + 2*r23)*(1/(r12 + r23) + 1/r12 + 1/r23)^2)/3 - (2*G^2*m^5*(1/(r12 + r23)^2 + 1/r12^2)*(1/(r12 + r23) + 1/r12 + 1/r23)*((r12 + r23)^2 + r12^2 + r23^2))/3
(G^2*m^5*(2*r12 + 4*r23)*(1/(r12 + r23) + 1/r12 + 1/r23)^2)/3 - (2*G^2*m^5*(1/(r12 + r23)^2 + 1/r23^2)*(1/(r12 + r23) + 1/r12 + 1/r23)*((r12 + r23)^2 + r12^2 + r23^2))/3


solve(eqn1,r12)
solve(eqn2,r23)

















