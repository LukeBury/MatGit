clear
clc
%%

% syms dt Ep y0 y1 y2 real
% 
% H0 = [1 0 0];
% H1 = [1 dt .5*dt^2];
% H2 = [1 2*dt 2*dt^2];
% 
% phi1 = [1 dt .5*dt^2; 0 1 dt; 0 0 1];
% phi2 = [1 2*dt 2*dt^2; 0 1 2*dt; 0 0 1];
% 
% H = [H0; H1*phi1; H2*phi2]
% 
% 
% % % This H is a stacked H*phi
% % H = [1 0 0; 1 2*dt 2*dt^2; 1 4*dt 6*dt^2]
% 
% R = eye(3)*Ep;
% 
% y = [y0;y1;y2];
% 
% P0 = inv(H'*inv(R)*H)
% 
% xhat0 = simplify(P0*H'*inv(R)*y)
% 
% syms x0 xdot0 a real
% y0 = x0;
% y1 = x0 + xdot0*dt + .5*a*dt^2;
% y2 = x0 + xdot0*2*dt + 2*a*dt^2;
% 
% subs(xhat0)

%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% trying without 'a' as a state
% syms dt Ep y0 y1 y2 real
% 
% H0 = [1 0];
% H1 = [1 dt];
% H2 = [1 2*dt];
% 
% phi10 = [1 dt; 0 1];
% phi20 = [1 2*dt; 0 1];
% 
% H = [H0; H1*phi10; H2*phi20]
% 
% R = eye(3)*Ep;
% 
% y = [y0;y1;y2];
% 
% P0 = inv(H'*inv(R)*H)
% 
% xhat0 = simplify(P0*H'*inv(R)*y)
% 
% syms x0 xdot0 a real
% y0 = x0;
% y1 = x0 + xdot0*dt + .5*a*dt^2;
% y2 = x0 + xdot0*2*dt + 2*a*dt^2;
% 
% subs(xhat0)

%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% trying w/ y = x, so H = [1 0]
% syms dt Ep y0 y1 y2 real
% 
% H0 = [1 0];
% H1 = [1 0];
% H2 = [1 0];
% 
% phi1 = [1 dt; 0 1];
% phi2 = [1 2*dt; 0 1];
% 
% H = [H0; H1*phi1; H2*phi2]
% 
% R = eye(3)*Ep;
% 
% y = [y0;y1;y2];
% 
% P0 = inv(H'*inv(R)*H)
% 
% xhat0 = simplify(P0*H'*inv(R)*y);
% 
% syms x0 x1 x2
% y0 = x0;
% y1 = x1;
% y2 = x2;
% 
% subs(xhat0)

%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% % syms dt sig y0 y1 y2 real
% % 
% % H0 = [1 0 0];
% % H1 = [1 0 0];
% % H2 = [1 0 0];
% % 
% % phi10 = [1 dt .5*dt^2; 0 1 dt; 0 0 1];
% % phi20 = [1 2*dt 2*dt^2; 0 1 2*dt; 0 0 1];
% % 
% % H = [H0; H1*phi10; H2*phi20]
% % 
% % R = eye(3)*sig;
% % 
% % y = [y0;y1;y2];
% % 
% % P0 = inv(H'*inv(R)*H)
% % 
% % xhat0 = simplify(P0*H'*inv(R)*y)
% % 
% % syms x0 xdot0 a real
% % y0 = x0;
% % y1 = x0 + xdot0*dt + .5*a*dt^2;
% % y2 = x0 + xdot0*2*dt + 2*a*dt^2;
% % 
% % subs(xhat0)

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
syms dt sig y0 y1 y2 real

H0 = [1 0];
H1 = [1 0];
H2 = [1 0];

phi10 = [1 dt; 0 1];
phi20 = [1 2*dt; 0 1];

H = [H0; H1*phi10; H2*phi20]

R = eye(3)*sig;

y = [y0;y1;y2];

P0 = inv(H'*inv(R)*H)

xhat0 = simplify(P0*H'*inv(R)*y)

syms x0 xdot0 a real
y0 = x0;
y1 = x0 + xdot0*dt + .5*a*dt^2;
y2 = x0 + xdot0*2*dt + 2*a*dt^2;

subs(xhat0)


