clear
clc
close all
addpath(genpath('/Users/lukebury/Documents/MATLAB/mbin'))
tic
% ========================================================================
%%% Run/Plot Switches
% ========================================================================

% ========================================================================
%%% Importing Data
% ========================================================================
%%% General data on solar system bodies
bodies = getBodyData();

%%% Color options/schemes
colors = get_colors();

% ========================================================================
%%% Running
% ========================================================================
A = [5 0 3 0 1;...
    3 0 0 -2 0;...
    0 -2 4 1 0;...
    1 3 -4 1 3;...
    0 2 2 0 -1];

B = [0 1;...
    0 2;...
    0 0;...
    1 3;...
    1 1];

C = [0 0 1 0 0];
D = [0 0];

x0 = [1 1 1 1 1]';

eps = 1e-12;
gammas = [1,1,1];
as =     [1,eps,eps];
bs =     [1,eps,1];
strings = {'\gamma=1   a=1   b=1';...
           '\gamma=1   a=0   b=0';...
           '\gamma=1   a=0   b=1'};
for kk = 1:3
gamma = gammas(kk);
a = as(kk);
b = bs(kk);

R = [a 0; 0 b];

%  Q is so that x' Q x = gamma*x3^2
Q = [0 0 0 0 0;...
    0 0 0 0 0;...
    0 0 gamma 0 0;...
    0 0 0 0 0;
    0 0 0 0 0];

N = 0; 
% bc we don't have an N in our cost function, though in full 
% form, J = Integral {x'Qx + u'Ru + 2*x'Nu} dt



% [K,S,E] = lqr(A,B,Q,R,N) is an equivalent syntax for continuous-time
%     models with dynamics  dx/dt = Ax + Bu
[K,S,E] = lqr(A,B,Q,R,N);

Anew = A - B*K;
% u = -K * x;
sys = ss(Anew,B,C,D);

time = [0:.0001:4];
u = [0; 0]; % so Bu is [5x1]
u = repmat(u,1,length(time));
% xd = Ax + Bu
% xd = Ax + B(-Kx)
% xd = (A - BK)x
% so u is zero and A is (A-BK)


[y,t,x] = lsim(sys,u,time,x0);

us = -K * x';





figure(1)
subplot(3,2,[1 3 5]); hold all
plot(time,y,'linewidth',2)
legend(strings)
PlotBoi2('Time, s', 'x$_3$',14,'LaTex')

if kk == 1
    subplot(3,2,2); hold all
    plot(time,us(1,:),'linewidth',2,'color',colors.std.mag)
    plot(time,us(2,:),'linewidth',2,'color',colors.std.grn)
    legend('u1','u2')
    PlotBoi2('', 'Control',14,'LaTex')
    title('$\gamma$=1   a=1   b=1','Interpreter','LaTex')
elseif kk == 2
    subplot(3,2,4); hold all
    plot(time,us(1,:),'linewidth',2,'color',colors.std.mag)
    plot(time,us(2,:),'linewidth',2,'color',colors.std.grn)
    PlotBoi2('', 'Control',14,'LaTex')
    title('$\gamma$=1   a=0   b=0','Interpreter','LaTex')
elseif kk == 3
    subplot(3,2,6); hold all
    plot(time,us(1,:),'linewidth',2,'color',colors.std.mag)
    plot(time,us(2,:),'linewidth',2,'color',colors.std.grn)
    PlotBoi2('Time, s', 'Control',14,'LaTex')
    title('$\gamma$=1   a=0   b=1','Interpreter','LaTex')
end


end





matrix2latex(Q,'test');
type test


































