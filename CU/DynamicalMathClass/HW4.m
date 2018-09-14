clear
clc

% % ------------------------------------------------------------------------
% %%% 4.15
% % ------------------------------------------------------------------------
% syms x y z lam sig b r real
% % dLdx = x*(1-lam*sig)
% % dLdy = y*(1-lam)
% % dLdz = 2*sig*z*(1-lam*b) + 2*sig*(lam*b*r-2*r)
% % dLdlam = -sig*r*x^2 - sig*y^2 - sig*b*z^2 + 2*sig*r*b*z
% 
% E = r*x^2 + sig*y^2 + sig*(z-2*r)^2;
% dLdx = 2*r*x - 2*lam*sig*r*x;
% dLdy = 2*sig*y - 2*lam*sig*y;
% dLdz = 2*sig*(z - 2*r) - 2*lam*sig*b*z + 2*lam*sig*r*b;
% dLdlam = -sig*r*x^2 - sig*y^2 - sig*b*z^2 + 2*sig*r*b*z;
% 
% sol = solve([dLdx == 0, dLdy == 0, dLdz == 0, dLdlam == 0],[x,y,z,lam]);
% xSol = simplify(sol.x);
% ySol = simplify(sol.y);
% zSol = simplify(sol.z);
% lSol = simplify(sol.lam);
% 
% % %%% ex from pg 145
% % E = x^2 + y^2 + (z - r - sig)^2;
% % dLdx   = x*(1-lam*sig)
% % dLdy   = y*(1-lam)
% % dLdz   = 2*z*(1-lam*b) - (2-lam*b)*(r+sig)
% % dLdlam = -sig*x^2 - y^2 - b*z^2 + (r+sig)*b*z
% % sol = solve([dLdx == 0, dLdy == 0, dLdz == 0, dLdlam == 0],[x,y,z,lam])
% % xSol = simplify(sol.x);
% % ySol = simplify(sol.y);
% % zSol = simplify(sol.z);
% % lSol = simplify(sol.lam);
% 
% for i = 1:length(xSol)
%     C(i) = simplify(subs(E,[x,y,z],[xSol(i), ySol(i), zSol(i)]));
%     C(i) = simplify(C(i));
% 
% end
% C(:)
% 
% sig = 10;
% b = 8/3; 
% % r = 28;
% r = .1;
% C1_eval = subs(C(1));
% C2_eval = double(subs(C(2)));
% C3_eval = double(subs(C(3)));
% C4_eval = double(subs(C(5)));
% 
% myC = double([C1_eval, C2_eval, C3_eval, C4_eval])
% 
% %%% comparing to txtbook
% R = ((r+sig)*.5*[0, 2, b/sqrt(b-1), b/sqrt(sig*(b-sig))]).^2


% % ------------------------------------------------------------------------
% %%% 4.17
% % ------------------------------------------------------------------------
% syms r r0 real
% A = [(r/r0)*exp(-4*pi) 0; 0 1];
% [v,e] = eig(A)



% % ------------------------------------------------------------------------
% %%% 4.18
% % ------------------------------------------------------------------------
% syms a b real
% A1 = [0 1 0; 1 -a 0; 0 0 -b];
% [v1,e1] = eig(A1)


% ------------------------------------------------------------------------
%%% 5.2
% ------------------------------------------------------------------------
syms a x y xs ys real

%%%%%%%%%%%%%%%%%%%%%%%%%% Part a

% yd = (x - 3*a*x^2 + .25*a*x^4*(.25*a^2*x^8 - x^2) + a*(.5*a^2*x^4)*x^3);
% q = -4 + 12*a*x + a^2*x^6;
% test = yd/q;
% simplify(test)

% xs = x;
% ys = .5*a*x^4;
% 
% yd = xs - 3*a*xs^2 + .5*ys*(ys^2-xs^2) + a*ys*xs^3;
% q = -4 + 12*a*x + a^2*x^6;
% test1 = subs(yd);
% test2 = test1/q

%%%%%%%%%%%%%%%%%%%%%%%%%% Part b
% dx = y + .5*x*(y^2-x^2) + a*x^4;
% dy = x - 3*a*x^2 + .5*y*(y^2-x^2) + a*y*x^3;
% 
% dxd_dx = diff(dx,x);
% dxd_dy = diff(dx,y);
% dyd_dx = diff(dy,x);
% dyc_dy = diff(dy,y);
% 
% Df = [dxd_dx, dxd_dy; dyd_dx, dyc_dy];
% x = 0; y = 0;
% Df_origin = subs(Df);
% 
% [v0, e0] = eig(Df_origin)

%%%%%%%%%%%%%%%%%%%%%%%%%% Part c

q = -4 + 12*x + x^6; % 1 0 0 0 0 12 -4
rts = roots([1 0 0 0 0 12 -4])
rt1 = real(rts(1))
rt2 = real(rts(end))

dx = y + .5*x*(y^2-x^2) + x^4;
dy = x - 3*x^2 + .5*y*(y^2-x^2) + y*x^3;

dxd_dx = diff(dx,x);
dxd_dy = diff(dx,y);
dyd_dx = diff(dy,x);
dyd_dy = diff(dy,y);

Df = [dxd_dx, dxd_dy; dyd_dx, dyd_dy];
x = rt1; y = .5*rt1^4;
Df_rt1 = vpa(subs(Df))

x = rt2; y = .5*rt2^4;
Df_rt2 = vpa(subs(Df))

[v_rt1, e_rt1] = eig(Df_rt1)

[v_rt2, e_rt2] = eig(Df_rt2)

% ------------------------------------------------------------------------
%%% 5.3
% ------------------------------------------------------------------------
% syms x y z d real
% xd = x*(1 - x - (1+d)*y - (1-d)*z);
% yd = y*(1 - y - (1+d)*z - (1-d)*x);
% zd = z*(1 - z - (1+d)*x - (1-d)*y);
% 
% Df = [diff(xd,x), diff(xd,y), diff(xd,z);...
%       diff(yd,x), diff(yd,y), diff(yd,z);...
%       diff(zd,x), diff(zd,y), diff(zd,z)];
% %%% Eq 1
% x = 0; y = 0; z = 0;
% Df_eq1 = subs(Df)
% [v1,e1] = eig(Df_eq1)
% 
% %%% Eq 2
% x = 0; y = 0; z = 1;
% Df_eq2 = subs(Df)
% [v2,e2] = eig(Df_eq2)
% 
% %%% Eq 3
% x = 0; y = 1; z = 0;
% Df_eq3 = subs(Df) 
% [v3,e3] = eig(Df_eq3)
% 
% %%% Eq 4
% x = 1; y = 0; z = 0;
% Df_eq4 = subs(Df)
% [v4,e4] = eig(Df_eq4)
% 
% %%% Eq 5
% x = 0; y = 1/d; z = -1/d;
% Df_eq5 = subs(Df)
% [v5,e5] = eig(Df_eq5)
% 
% %%% Eq 6
% x = -1/d; y = 0; z = 1/d;
% Df_eq6 = subs(Df)
% [v6,e6] = eig(Df_eq6)
% 
% %%% Eq 7
% x = 1/d; y = -1/d; z = 0;
% Df_eq7 = subs(Df)
% [v7,e7] = eig(Df_eq7)
% 
% %%% Eq 8
% x = 1/3; y = 1/3; z = 1/3;
% Df_eq8 = subs(Df)
% [v8,e8] = eig(Df_eq8)
% 
% H = d*x*y*(1-x-y)
% xd = diff(H,y)
% yd = -diff(H,x)

% ------------------------------------------------------------------------
%%% 5.4
% ------------------------------------------------------------------------
% A = [-1 0; 0 1];
% [vA,eA] = eig(A)
% 
% 
% syms s t real
% 
% eqn = expm(t*A-s*A)*[0; sin(s)]
% 
% int(eqn,s,t,inf)


% ------------------------------------------------------------------------
%%% 5.4
% ------------------------------------------------------------------------
% %%% Part A
% Df00 = [-1 0; 0 2];
% [v,e] = eig(Df00)
% 
% %%% Part B
% syms t sx s real
% A = [-1 0; 0 2];
% T0 = expm(t*A)*[sx;0]
% 
% g1 = [0; sx^2*exp(-2*s)];
% piS = [1 0; 0 0];
% piU = [0 0; 0 1];
% T1 = expm(t*A)*[sx;0] + int(expm((t-s)*A)*piS*g1,s,0,t) - int(expm((t-s)*A)*piU*g1,s,t,inf)
% 













