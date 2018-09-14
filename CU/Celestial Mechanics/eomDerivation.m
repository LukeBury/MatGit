clear
clc
tic %
syms m1 m2 m3 m4 m5 real
syms R1 R2 R3 R4 real
syms r12 r23 r34 r45
syms G real

toc %%
U = -G*(m1*m2/R1 + m1*m3/(R2 + m2*R1/(m1+m2)) + m1*m4/(R3 + m3*R2/(m1+m2+m3) + m2*R1/(m1+m2)) + m1*m5/(R4 + m4*R3/(m1+m2+m3+m4) + m3*R2/(m1+m2+m3) + m2*R1/(m1+m2)) + m2*m3/(R2 - m1*R1/(m1+m2)) + m2*m4/(R3 + m3*R2/(m1+m2+m3) - m1*R1/(m1+m2)) + m2*m5/(R4 + m4*R3/(m1+m2+m3+m4) + m3*R2/(m1+m2+m3) - m1*R1/(m1+m2)) + m3*m4/(R3 - (m1+m2)*R2/(m1+m2+m3)) + m3*m5/(R4 + m4*R3/(m1+m2+m3+m4) - (m1+m2)*R2/(m1+m2+m3)) + m4*m5/(R4 - (m1+m2+m3)*R3/(m1+m2+m3+m4)));

toc %%
dR1 = -diff(U,R1);
dR2 = -diff(U,R2);
dR3 = -diff(U,R3);
dR4 = -diff(U,R4);

toc %%
R1 = r12;
R2 = r23 + R1*m1/(m1+m2);
R3 = r34 + (m1+m2)*R2/(m1+m2+m3);
R4 = r45 + (m1+m2+m3)*R3/(m1+m2+m3+m4);

toc %%
dR1_sub = subs(dR1);
dR2_sub = subs(dR2);
dR3_sub = subs(dR3);
dR4_sub = subs(dR4);

toc %%
sdR1 = simplify(dR1_sub)
sdR2 = simplify(dR2_sub)
sdR3 = simplify(dR3_sub)
sdR4 = dR4_sub
% sdR4 = simplify(dR4_sub)    this gets pretty nasty

toc %%
ddR1 = sdR1*(m1+m2)/(m2*m1)
%-(G*(m1 + m2)*((m1*m2)/r12^2 - (m1*m2*m3)/(r23^2*(m1 + m2)) + (m1*m2*m4)/((m1 + m2)*(r12 + r23 + r34)^2) - (m1*m2*m5)/((m1 + m2)*(r23 + r34 + r45)^2) + (m1*m2*m5)/((m1 + m2)*(r12 + r23 + r34 + r45)^2) + (m1*m2*m3)/((m1 + m2)*(r12 + r23)^2) - (m1*m2*m4)/((m1 + m2)*(r23 + r34)^2)))/(m1*m2)
ddR2 = sdR2*(m1+m2+m3)/(m3*(m1+m2))
%-(G*(m1 + m2 + m3)*((m2*m3)/r23^2 + (m1*m3)/(r12 + r23)^2 - (m3*m5*(m1 + m2))/((r34 + r45)^2*(m1 + m2 + m3)) + (m1*m3*m5)/((m1 + m2 + m3)*(r12 + r23 + r34 + r45)^2) - (m3*m4*(m1 + m2))/(r34^2*(m1 + m2 + m3)) + (m2*m3*m4)/((r23 + r34)^2*(m1 + m2 + m3)) + (m1*m3*m4)/((m1 + m2 + m3)*(r12 + r23 + r34)^2) + (m2*m3*m5)/((m1 + m2 + m3)*(r23 + r34 + r45)^2)))/(m3*(m1 + m2))
ddR3 = sdR3*(m1+m2+m3+m4)/(m4*(m1+m2+m3))
%-(G*(m1 + m2 + m3 + m4)*((m3*m4)/r34^2 + (m1*m4)/(r12 + r23 + r34)^2 + (m2*m4)/(r23 + r34)^2 - (m4*m5*(m1 + m2 + m3))/(r45^2*(m1 + m2 + m3 + m4)) + (m2*m4*m5)/((r23 + r34 + r45)^2*(m1 + m2 + m3 + m4)) + (m1*m4*m5)/((m1 + m2 + m3 + m4)*(r12 + r23 + r34 + r45)^2) + (m3*m4*m5)/((r34 + r45)^2*(m1 + m2 + m3 + m4))))/(m4*(m1 + m2 + m3))
ddR4 = sdR4*(m1+m2+m3+m4+m5)/(m5*(m1+m2+m3+m4))
%-(G*((m4*m5)/r45^2 + (m3*m5)/(r45 - ((m1 + m2)*(r23 + (m1*r12)/(m1 + m2)))/(m1 + m2 + m3) + (m4*(r34 + ((m1 + m2)*(r23 + (m1*r12)/(m1 + m2)))/(m1 + m2 + m3)))/(m1 + m2 + m3 + m4) + ((r34 + ((m1 + m2)*(r23 + (m1*r12)/(m1 + m2)))/(m1 + m2 + m3))*(m1 + m2 + m3))/(m1 + m2 + m3 + m4))^2 + (m1*m5)/(r45 + (m4*(r34 + ((m1 + m2)*(r23 + (m1*r12)/(m1 + m2)))/(m1 + m2 + m3)))/(m1 + m2 + m3 + m4) + (m3*(r23 + (m1*r12)/(m1 + m2)))/(m1 + m2 + m3) + ((r34 + ((m1 + m2)*(r23 + (m1*r12)/(m1 + m2)))/(m1 + m2 + m3))*(m1 + m2 + m3))/(m1 + m2 + m3 + m4) + (m2*r12)/(m1 + m2))^2 + (m2*m5)/(r45 + (m4*(r34 + ((m1 + m2)*(r23 + (m1*r12)/(m1 + m2)))/(m1 + m2 + m3)))/(m1 + m2 + m3 + m4) + (m3*(r23 + (m1*r12)/(m1 + m2)))/(m1 + m2 + m3) + ((r34 + ((m1 + m2)*(r23 + (m1*r12)/(m1 + m2)))/(m1 + m2 + m3))*(m1 + m2 + m3))/(m1 + m2 + m3 + m4) - (m1*r12)/(m1 + m2))^2)*(m1 + m2 + m3 + m4 + m5))/(m5*(m1 + m2 + m3 + m4))
r13 = r12 + r23;
r14 = r13 + r34;
r15 = r14 + r45;
r24 = r23 + r34;
r25 = r24 + r45;
r35 = r34 + r45;
tt = subs(ddR4);
tt
tt_s = simplify(tt)



% r12 = R1;
% r13 = R2 + m2*R1/(m1+m2);
% r14 = R3 + m3*R2/(m1+m2+m3) + m2*R1/(m1+m2);
% r15 = R4 + m4*R3/(m1+m2+m3+m4) + m3*R2/(m1+m2+m3) + m2*R1/(m1+m2);

% r23 = R2 - m1*R1/(m1+m2);
% r24 = R3 + m3*R2/(m1+m2+m3) - m1*R1/(m1+m2);
% r25 = R4 + m4*R3/(m1+m2+m3+m4) + m3*R2/(m1+m2+m3) - m1*R1/(m1+m2);

% r34 = R3 - (m1+m2)*R2/(m1+m2+m3);
% r35 = R4 + m4*R3/(m1+m2+m3+m4) - (m1+m2)*R2/(m1+m2+m3);

% r45 = R4 - (m1+m2+m3)*R3/(m1+m2+m3+m4);