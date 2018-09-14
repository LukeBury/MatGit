clear
clc

fprintf('------------------------ 1ii ------------------------\n')
f = 0.8;
q = 10;

%%% From (i)
Pinf_i = q/(1-f^2)

%%% From (ii)
sig2 = 50;
Pinf_ii = sig2^2;
for i = 1:50
    Pinf_ii = f*Pinf_ii*f + q;
end
fprintf('With sig2 = 50, P50 = %f\n',Pinf_ii)

sig2 = 10;
Pinf_ii = sig2^2;
for i = 1:50
    Pinf_ii = f*Pinf_ii*f + q;
end
fprintf('With sig2 = 10, P10 = %f\n',Pinf_ii)

fprintf('------------------------ 1iv ------------------------\n')
P0 = 10*eye(2);
F = [0.99 0.2; 0 -0.76];
Q = [1 0.37; 0.37 2.5];
Pinf_iv = P0;
for i = 1:500
    Pinf_iv = F*Pinf_iv*F' + Q;
end
fprintf('Sufficiently long predicted covariance sim:\nP = \n')
fprintf('%1.4f %1.4f\n',Pinf_iv)

Pinf_dlyap = dlyap(F,Q);
fprintf('\nUsing dlyap.m:\nP = \n')
fprintf('%1.4f %1.4f\n',Pinf_dlyap)

f11 = F(1,1);
f12 = F(1,2);
f21 = F(2,1);
f22 = F(2,2);

Q = [Q(1,1); Q(1,2); Q(2,2)];
mat = [1-f11^2, -2*f11*f12, -f12^2;...
    -f11*f21, 1-f11*f22-f12*f21, -f12*f22;...
    -f21^2, -2*f21*f22, 1-f22^2];

P = mat\Q;
P = [P(1) P(2); P(2) P(3)];
fprintf('\nSolving analytically:\nP = \n')
fprintf('%1.4f %1.4f\n',P)









